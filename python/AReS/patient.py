from collections import deque
from scipy.integrate import solve_ivp

import logging
import numpy as np
import scipy.signal as signal

from .disturbance import Disturbance
from .pharmacokinetics import Pharmacokinetic
from .pharmacodynamics import PharmacodynamicDoH, PharmacodynamicHemo, PharmacodynamicNMB
from .utils.enums import Interaction, DoHMeasure, Model, Drug, PatientPhase, BloodSampling


class Patient:
    def __init__(self, data: list,
                 internal_states: dict = None,
                 output_init: dict = None,
                 pk_models: dict = None,
                 pd_models: dict = None,
                 opiates: bool = True,
                 blood_sampling: BloodSampling = BloodSampling.ARTERIAL,
                 interaction: Interaction = Interaction.SURFACE,
                 doh_measure: DoHMeasure = DoHMeasure.BOTH,
                 disturbance_model: Disturbance = None,
                 volume_status: dict = None,
                 seed=None):

        # Patient demographic parameters
        self._height = data[0]
        self._weight = data[1]
        self._age = data[2]
        self._gender = data[3]  # 1 for female and 0 for male
        self._bmi = data[4]
        self._lbm = data[5]

        opiates = opiates if opiates is not None else True  # If opiates are used or not during the surgery
        blood_sampling = blood_sampling if blood_sampling is not None else BloodSampling.ARTERIAL  # The blood sampling site: arterial or venous

        # Initialize the internal states of the patient
        if internal_states is None:
            internal_states = {}
        self._x_prop = internal_states.get('x_prop', [0, 0, 0])             # Initial states (concentrations) of propofol pk model
        self._x_ce_prop = internal_states.get('x_ce_prop', 0)               # Initial state of propofol effect-site concentration
        self._x_delay = internal_states.get('x_delay',0)                    # Initial state of the delay part of propofol pd model

        self._x_wav_filtered =  [0, 0]                                      # Initial state of the filter representing wav monitor
        self._x_bis_lti = [0, 0, 0, 0]                                      # Initial state of the LTI part of the filter representing bis monitor
        self._x_bis_delay = 0                                                # Initial state of the delay part of the filter representing bis monitor

        self._x_remi = internal_states.get('x_remi', [0, 0, 0])             # Initial states (concentrations) of remifentanil pk model
        self._x_ce_remi = internal_states.get('x_ce_remi', 0)               # Initial state of remifentanil effect-site concentration

        self._x_nore = internal_states.get('x_nore', [0, 0, 0])             # Initial states (concentrations) of norepinephrine pk model
        self._x_nore_delayed = internal_states.get('x_nore_delayed', [0])   # Initial state of the delay part of norepinephrine pk model

        self._x_rocu = internal_states.get('x_rocu', [0, 0, 0])             # Initial states (concentrations) of rocuronium pk model
        self._x_ce_rocu = internal_states.get('x_ce_rocu', 0)               # Initial state of rocuronium effect-site concentration

        # Pk-models names for propofol, remifentanil, norepinephrine, and rocuronium
        if pk_models is None:
            pk_models = {}
        pk_model_prop = pk_models.get('prop', Model.ELEVELD)
        pk_model_remi = pk_models.get('remi', Model.ELEVELD)
        pk_model_nore = pk_models.get('nore', Model.JOACHIM)
        pk_model_rocu = pk_models.get('rocu', Model.DAHE)

        # Pd-models names for propofol and remifentanil
        if pd_models is None:
            pd_models = {}
        pd_model_prop = pd_models.get('prop', Model.ELEVELD)
        pd_model_remi = pd_models.get('remi', Model.ELEVELD)

        # Parameters to initialize the Pk and Pd model
        kwargs = {'age': self._age, 'weight': self._weight, 'height': self._height, 'gender': self._gender,
                  'bmi': self._bmi, 'lbm': self._lbm, 'opiates': opiates, 'blood_sampling': blood_sampling}

        if pd_model_prop == Model.PATIENT_SPECIFIC:
            # Parameters to define the Propofol patient-specific PD model
            if len(data) < 11:
                raise ValueError(
                    'The data list should contain 11 elements to initialize the patient with a patient-specific PD model.')
            pd_patient_specific_parameter = {
                'e0': data[6],
                'ke0': data[7],
                'delay': data[8],
                'ec50': 0.01 * data[6] * data[9],
                'gamma': data[10]
            }
            kwargs['pd_patient_specific_parameter'] = pd_patient_specific_parameter

        # Pharmacokinetic model objects
        self._pk_prop = Pharmacokinetic.create(Drug.PROPOFOL, pk_model_prop, **kwargs)
        self._pk_remi = Pharmacokinetic.create(Drug.REMIFENTANIL, pk_model_remi, **kwargs)
        self._pk_nore = Pharmacokinetic.create(Drug.NOREPINEPHRINE, pk_model_nore, **kwargs)
        self._pk_rocu = Pharmacokinetic.create(Drug.ROCURONIUM, pk_model_rocu, **kwargs)

        # Pharmacodynamic model objects
        self._pd_doh = PharmacodynamicDoH(model_prop=pd_model_prop, model_remi=pd_model_remi, **kwargs)
        if output_init is not None and 'doh' in output_init.keys():
            self._pd_doh.set_e0(output_init.get('doh'))

        # Initialize the hemodynamic model object
        self._pd_hemo = PharmacodynamicHemo(self._age)

        # Initialize the hemodynamic variables

        # If the user defines the initial values of hemodynamic variables, the base values in the hemodynamic model introduced
        # by Su et al.(2023) is calculated from the user defined hemoynamic variables; otherwise, base values are the ones
        # introduced in the study by Su et al. (2023).

        if output_init is not None:
            hr = output_init.get('hr') if 'hr' in output_init.keys() else self._pd_hemo.base_hr
            base_hr = hr - self._pd_hemo.base_hr * self._pd_hemo.ltde_hr

            sv = output_init.get('sv') if 'sv' in output_init.keys() else self._pd_hemo.base_sv
            base_sv = sv / (1 - self._pd_hemo.hr_sv * np.log(hr / (self._pd_hemo.base_hr * (
                    1 + self._pd_hemo.ltde_hr)))) - self._pd_hemo.base_sv * self._pd_hemo.ltde_sv

            co = output_init.get(
                'co') if 'co' in output_init.keys() else self._pd_hemo.base_hr * self._pd_hemo.base_sv / 1000

            if co != base_hr * base_sv / 1000:
                logging.warning(
                    f'Inconsistency in initial values: CO does not match the expected value derived from HR and SV. '
                    f'Computed error: {co - (base_hr * base_sv / 1000)}. This discrepancy may result in incorrect initialization of MAP.')

            map_ = output_init.get(
                'map') if 'map' in output_init.keys() else self._pd_hemo.base_hr * self._pd_hemo.base_sv * self._pd_hemo.base_tpr
            base_tpr = map_ / (co * 1000)

            self._hemodynamic_variables = [base_tpr, base_sv, base_hr, 0, 0]
        else:
            if seed is not None:
                # Add interpatient variability to the hemodynamic variables if seed is not None
                self._pd_hemo.interpatient_variability(seed)
            self._hemodynamic_variables = [self._pd_hemo.base_tpr, self._pd_hemo.base_sv, self._pd_hemo.base_hr, 0, 0]

        # Initialize the neuromuscular blockade model
        self._pd_nmb = PharmacodynamicNMB()

        # Patient simulation settings
        self._interaction = interaction
        self._doh_measure = doh_measure
        self._time = 0

        # Disturbance model and initial values of the disturbances lists
        self._disturbance_model = disturbance_model
        self._doh_dis = [0]
        self._hr_dis = [0]
        self._map_dis = [0]

        # Volume status model
        self._volume_status = volume_status
        self._volume_status_coeff = {'sv': 1, 'co': 1, 'map': 1, 'hr': 1}

        # Parameters to define the patient phase and the steady state (needed for the PID controller)
        self._patient_phase = PatientPhase.INDUCTION
        self._maintenance_time = 0
        self._steady_state = False
        self._steady_state_values = deque(maxlen=180)

        # Create the arrays to record the plasma and effect-site concentration, infusion rates, and the predicted stats
        self._cp_prop = np.array([], dtype=np.float64)  # Simulated propofol plasma concentraion
        self._ce_prop = np.array([], dtype=np.float64)  # Simulated propofol effect-site concentraion
        self._ce_del = np.array([], dtype=np.float64)   # Simulated variable that represent the effect of delay  in pd model on propofol effect-site concentraion
        self._ce_wav = np.array([], dtype=np.float64)   # Simulated variable that represent the effect of WAV filter on propofol effect-site concentraion
        self._ce_bis = np.array([], dtype=np.float64)   # Simulated variable that represent the effect of BIS filter on propofol effect-site concentraion

        self._cp_remi = np.array([], dtype=np.float64)  # Simulated remifentanil plasma concentraion
        self._ce_remi = np.array([], dtype=np.float64)  # Simulated remifentanil effect-site concentraion

        self._c_nore = np.array([], dtype=np.float64)   # Simulated norepinephrine blood concentraion
        self._cp_rocu = np.array([], dtype=np.float64)  # Simulated rocuronium plasma concentraion
        self._ce_rocu = np.array([], dtype=np.float64)  # Simulated rocuronium effect-site concentraion

        self._u_prop = []                                     # Propofol infusion rate
        self._u_remi = []                                     # Remifentanil infusion rate
        self._u_nore = []                                     # Norepinephrine infusion rate
        self._u_rocu = []                                     # Rocuronium infusion rate

        # Doh
        self._WAV = np.array([], dtype=np.float64)      # Simulated WAV Index
        self._BIS = np.array([], dtype=np.float64)      # Simulated BIS Index

        # Hemodynamic variables
        self._MAP = np.array([], dtype=np.float64)      # Simulated Mean Arterial Pressure
        self._CO = np.array([], dtype=np.float64)       # Simulated Cardiac Output
        self._HR = np.array([], dtype=np.float64)       # Simulated Heart Rate
        self._SV = np.array([], dtype=np.float64)       # Simulated Stroke Volume
        self._TPR = np.array([], dtype=np.float64)      # Simulated Total Peripheral Resistance

        # Neuromuscular blockade
        self._NMB_m0 = np.array([], dtype=np.float64)   # Simulated probability that m = 0 (NMB is moderate)
        self._NMB_m1 = np.array([], dtype=np.float64)   # Simulated probability that m = 1 (NMB is deep)
        self._NMB_m2 = np.array([], dtype=np.float64)   # Simulated probability that m = 2 (NMB is profound)
        self._NMB_m3 = np.array([], dtype=np.float64)   # Simulated probability that m = 3 (NMB is very profound)

    def get_patient_demographics(self):
        """
        :returns: The patient demographic parameters: age, height, weight, gender, bmi, lbm
        :rtype: dict
        """
        patient_data = {
            'age': self._age,
            'height': self._height,
            'weight': self._weight,
            'gender': self._gender,
            'bmi': self._bmi,
            'lbm': self._lbm
        }
        return patient_data

    def get_patient_state(self):
        """
        :returns: The current state of the patient.
            - WAV, BIS, MAP, CO, HR, SV, NMB_m0, NMB_m1, NMB_m2, NMB_m3, cp_prop, cp_remi, c_nore, cp_rocu, ce_prop,
                ce_remi, ce_del, ce_wav, ce_bis
        :rtype: dict
        """
        base_wav = self._pd_doh.get_e0()
        base_bis = self._pd_doh.get_e0()
        base_co = self._pd_hemo.base_sv * self._pd_hemo.base_hr / 1000
        base_map = self._pd_hemo.base_sv * self._pd_hemo.base_hr * self._pd_hemo.base_tpr
        base_hr = self._pd_hemo.base_hr
        base_sv = self._pd_hemo.base_sv

        state = {
            'WAV': self._WAV[-1] if len(self._WAV) > 0 else base_wav,
            'BIS': self._BIS[-1] if len(self._BIS) > 0 else base_bis,
            'MAP': self._MAP[-1] if len(self._MAP) > 0 else base_map,
            'CO': self._CO[-1] if len(self._CO) > 0 else base_co,
            'HR': self._HR[-1] if len(self._HR) > 0 else base_hr,
            'SV': self._SV[-1] if len(self._SV) > 0 else base_sv,
            'TPR': self._TPR[-1] if len(self._TPR) > 0 else self._pd_hemo.base_tpr,
            'NMB_m0': self._NMB_m0[-1] if len(self._NMB_m0) > 0 else 1,
            'NMB_m1': self._NMB_m1[-1] if len(self._NMB_m1) > 0 else 0,
            'NMB_m2': self._NMB_m2[-1] if len(self._NMB_m2) > 0 else 0,
            'NMB_m3': self._NMB_m3[-1] if len(self._NMB_m3) > 0 else 0,
            'cp_prop': self._cp_prop[-1] if len(self._cp_prop) > 0 else 0,
            'cp_remi': self._cp_remi[-1] if len(self._cp_remi) > 0 else 0,
            'c_nore': self._c_nore[-1] if len(self._c_nore) > 0 else 0,
            'cp_rocu': self._cp_rocu[-1] if len(self._cp_rocu) > 0 else 0,
            'ce_prop': self._ce_prop[-1] if len(self._ce_prop) > 0 else 0,
            'ce_remi': self._ce_remi[-1] if len(self._ce_remi) > 0 else 0,
            'ce_rocu': self._ce_rocu[-1] if len(self._ce_rocu) > 0 else 0,
            'ce_del': self._ce_del[-1] if len(self._ce_del) > 0 else 0,
            'ce_wav': self._ce_wav[-1] if len(self._ce_wav) > 0 else 0,
            'ce_bis': self._ce_bis[-1] if len(self._ce_bis) > 0 else 0
        }

        return state

    def get_patient_state_history(self):
        """
        :returns: The history of the patient state.
        :rtype: dict
        """
        return {'WAV': self._WAV, 'BIS': self._BIS, 'MAP': self._MAP, 'CO': self._CO, 'HR': self._HR, 'SV': self._SV,
                'TPR': self._TPR, 'NMB_m0': self._NMB_m0, 'NMB_m1': self._NMB_m1, 'NMB_m2': self._NMB_m2,
                'NMB_m3': self._NMB_m3, 'cp_prop': self._cp_prop, 'cp_remi': self._cp_remi, 'c_nore': self._c_nore,
                'cp_rocu': self._cp_rocu, 'ce_prop': self._ce_prop, 'ce_remi': self._ce_remi, 'ce_rocu': self._ce_rocu,
                'ce_del': self._ce_del, 'ce_wav': self._ce_wav, 'ce_bis': self._ce_bis}

    def get_patient_input_history(self):
        """
        :returns: The history of the patient inputs.
            - u_prop (np.ndarray): The propofol input from the beginning of the simulation to the current time in mg * sec^-1.
            - u_remi (np.ndarray): The remifentanil input from the beginning of the simulation to the current time in µg * sec^-1.
            - u_nore (np.ndarray): The norepinephrine input from the beginning of the simulation to the current time in nmol * sec^-1.
            - u_rocu (np.ndarray): The rocuronium input from the beginning of the simulation to the current time in µg * sec^-1.
        :rtype: dict
        """
        return {'u_prop': self._u_prop, 'u_remi': self._u_remi, 'u_nore': self._u_nore, 'u_rocu': self._u_rocu}

    def get_patient_disturbances(self):
        """
        :returns: The disturbances of the patient for the doh and co.
        :rtype: dict
        """
        return {'doh': self._doh_dis, 'hr': self._hr_dis, 'map': self._map_dis}

    def get_patient_phase(self):
        """
        :returns: The current phase of the patient and if the patient is in steady state or not.
        :rtype: tuple
        """
        return self._patient_phase, self._steady_state

    def get_patient_internal_states(self):
        """
        :returns: The current patient internal states
        :rtype: dict
        """
        return {'x_prop': self._x_prop, 'x_ce_prop': self._x_ce_prop, 'x_delay':self._x_delay,
                'x_remi': self._x_remi,'x_ce_remi': self._x_ce_remi, 'x_nore': self._x_nore,
                'x_nore_delayed': self._x_nore_delayed, 'x_rocu': self._x_rocu, 'x_ce_rocu': self._x_ce_rocu}

    def set_time(self, time: int):
        """
        Set the current time of the simulation
        """
        self._time = time

    def set_disturbance_model(self, disturbance_model: Disturbance):
        """
        Set the disturbance model
        """
        self._disturbance_model = disturbance_model

    def set_volume_status(self, volume_status: dict):
        """
        Set the volume status
        """
        self._volume_status = volume_status

    def clear_data(self):
        """
        Clear the patient data except for the last values. Needed to start the simulation in maintenance phase.
        """
        self._cp_prop = self._cp_prop[-1:]
        self._ce_prop = self._ce_prop[-1:]
        self._ce_del = self._ce_del[-1:]
        self._ce_wav = self._ce_wav[-1:]
        self._ce_bis = self._ce_bis[-1:]

        self._cp_remi = self._cp_remi[-1:]
        self._ce_remi = self._ce_remi[-1:]

        self._c_nore = self._c_nore[-1:]
        self._cp_rocu = self._cp_rocu[-1:]
        self._ce_rocu = self._ce_rocu[-1:]

        self._u_prop = self._u_prop[-1:]
        self._u_remi = self._u_remi[-1:]
        self._u_nore = self._u_nore[-1:]
        self._u_rocu = self._u_rocu[-1:]

        self._WAV = self._WAV[-1:]
        self._BIS = self._BIS[-1:]

        self._MAP = self._MAP[-1:]
        self._CO = self._CO[-1:]
        self._HR = self._HR[-1:]
        self._SV = self._SV[-1:]
        self._TPR = self._TPR[-1:]

        self._NMB_m0 = self._NMB_m0[-1:]
        self._NMB_m1 = self._NMB_m1[-1:]
        self._NMB_m2 = self._NMB_m2[-1:]
        self._NMB_m3 = self._NMB_m3[-1:]

    def _compute_plasma_concentration(self, u_prop: np.ndarray, u_remi: np.ndarray, u_nore: np.ndarray,
                                      u_rocu: np.ndarray, t: np.array):
        """
        Compute the plasma concentration of propofol, remifentanil, norepinephrine, and rocuronium.
        :param u_prop: The propofol input [mg * sec^-1]
        :type u_prop: np.ndarray
        :param u_remi: The remifentanil input [µg * sec^-1]
        :type u_remi: np.ndarray
        :param u_nore: The norepinephrine input [µg * sec^-1]
        :type u_nore: np.ndarray
        :param u_rocu: The rocuronium input [mg * sec^-1]
        :type u_rocu: np.ndarray
        :param t: The time vector. [second]
        :type t: np.ndarray

        :returns: The plasma concentration of propofol, remifentanil, norepinephrine, and rocuronium.
        :rtype: tuple of np.ndarray
        """
        _, cp_prop_sim, x_prop_sim = signal.lsim(self._pk_prop.ss, u_prop, t, X0=self._x_prop)
        _, cp_remi_sim, x_remi_sim = signal.lsim(self._pk_remi.ss, u_remi, t, X0=self._x_remi)
        _, c_nore_sim, x_nore_sim = signal.lsim(self._pk_nore.ss, u_nore, t, X0=self._x_nore)
        _, cp_rocu_sim, x_rocu_sim = signal.lsim(self._pk_rocu.ss, u_rocu, t, X0=self._x_rocu)

        cp_prop_sim[cp_prop_sim < 0] = 0
        cp_remi_sim[cp_remi_sim < 0] = 0
        c_nore_sim[c_nore_sim < 0] = 0
        cp_rocu_sim[cp_rocu_sim < 0] = 0

        # Update initial states
        self._x_prop = x_prop_sim[-1, :]
        self._x_remi = x_remi_sim[-1, :]
        self._x_nore = x_nore_sim[-1, :]
        self._x_rocu = x_rocu_sim[-1, :]

        # Append the plasma concentrations to the history
        self._cp_prop = np.append(self._cp_prop, cp_prop_sim[:])
        self._cp_remi = np.append(self._cp_remi, cp_remi_sim[:])
        self._c_nore = np.append(self._c_nore, c_nore_sim[:])
        self._cp_rocu = np.append(self._cp_rocu, cp_rocu_sim[:])

        return cp_prop_sim, cp_remi_sim, c_nore_sim, cp_rocu_sim

    def _compute_WAV(self, ce_delayed_prop, ce_remi_sim: np.ndarray, t: np.ndarray):
        """
        Compute the WAV using the pharmacodynamic model.
        :param ce_delayed_prop: The delayed effect-site concentration of propofol.
        :type ce_delayed_prop: np.ndarray
        :param ce_remi_sim: The effect-site concentration of remifentanil.
        :type ce_remi_sim: np.ndarray
        :param t: The time vector.
        :type t: np.ndarray

        :returns: The WAV values.
        :rtype: np.ndarray
        """
        # Apply the filter defined for wav monitor
        _, ce_wav_filtered, x_wav_filtered = signal.lsim(PharmacodynamicDoH.neurow_sensor_dynamics(), ce_delayed_prop,
                                                         t, self._x_wav_filtered)

        ce_wav_filtered[ce_wav_filtered < 0] = 0
        self._x_wav_filtered = x_wav_filtered[-1, :]
        self._ce_wav = np.append(self._ce_wav, ce_wav_filtered[:])
        # If there is no interaction between propofol and remifentanil, the hill function is applied; otherwise the surface interaction model is used.
        if self._interaction == Interaction.NO_INTERACTION:
            wav = self._pd_doh.hillfun(ce_wav_filtered)
        elif self._interaction == Interaction.SURFACE:
            wav = self._pd_doh.responseSurfaceModel(ce_wav_filtered, ce_remi_sim)
        else:
            raise NotImplementedError("This interaction model is not implemented yet.")

        return wav

    def _compute_BIS(self, ce_delayed_prop, ce_remi_sim: np.ndarray,
                     t: np.ndarray):
        """
        Compute the BIS using the pharmacodynamic model.
        :param ce_delayed_prop: The delayed effect-site concentration of propofol.
        :type ce_delayed_prop: np.ndarray
        :param ce_remi_sim: The effect-site concentration of remifentanil.
        :type ce_remi_sim: np.ndarray
        :param t: The time vector.
        :type t: np.ndarray

        :returns: The BIS values.
          :rtype: np.ndarray
        """
        # Apply the filter defined for bis monitor

        bis_delay, bis_lti = PharmacodynamicDoH.bis_sensor_dynamics()

        _, ce_bis_filtered, x_bis_lti = signal.lsim(bis_lti, ce_delayed_prop, t, self._x_bis_lti)
        _, ce_bis_filtered, x_bis_delay = signal.lsim(bis_delay, ce_bis_filtered, t, self._x_bis_delay)

        ce_bis_filtered[ce_bis_filtered < 0] = 0
        self._ce_bis = np.append(self._ce_bis, ce_bis_filtered[:])

        self._x_bis_lti = x_bis_lti[-1, :]
        self._x_bis_delay = x_bis_delay[-1]  # this is correct, x_bis is a 1D array

        # If there is no interaction between propofol and remifentanil, the hill function is applied; otherwise the surface interaction model is used.
        if self._interaction == Interaction.NO_INTERACTION:
            bis = self._pd_doh.hillfun(
                ce_bis_filtered)  # This is not an error, the ce_bis_filtered is a ndarray of shape (t,1)
        elif self._interaction == Interaction.SURFACE:
            bis = self._pd_doh.responseSurfaceModel(ce_bis_filtered, ce_remi_sim)
        else:
            raise NotImplementedError("This interaction model is not implemented yet.")

        return bis

    def _compute_hemodynamic_variables(self, cp_prop: np.ndarray, cp_remi: np.ndarray, c_nore: np.ndarray,
                                       t: np.ndarray):
        """
        Compute the hemodynamic variables using the pharmacodynamic model.
        :param cp_prop: The plasma concentration of propofol.
        :type cp_prop: np.ndarray
        :param cp_remi: The plasma concentration of remifentanil.
        :type cp_remi: np.ndarray
        :param c_nore: The blood concentration of norepinephrine.
        :type c_nore: np.ndarray
        :param t: The time vector.
        :type t: np.ndarray
        :returns: The cardiac output, mean arterial pressure, heart rate, and stroke volume.
        :rtype: tuple of np.ndarray
        """

        # Norepinephrine's effect on hemodynamic variables
        _, c_nore_delayed_sim, x_nore_delayed_sim = signal.lsim(self._pd_hemo.get_nore_delay_ss(), c_nore, t,
                                                                X0=self._x_nore_delayed)
        c_nore_delayed_sim[c_nore_delayed_sim < 0] = 0
        self._x_nore_delayed = x_nore_delayed_sim[-1]  # This is correct, x_nore_delayed is a 1D array

        map_nore = self._pd_hemo.hillfun(c_nore_delayed_sim, 'map')
        co_nore = self._pd_hemo.hillfun(c_nore_delayed_sim, 'co')

        # General pharmacodynamic interaction (GPDI) model between Propofol and Remifentanil

        # Initialize the hemodynamic parameters
        tpr_interval = [self._hemodynamic_variables[0]]
        sv_star_interval = [self._hemodynamic_variables[1]]
        hr_star_interval = [self._hemodynamic_variables[2]]
        tde_hr_interval = [self._hemodynamic_variables[3] + self._pd_hemo.base_hr * self._pd_hemo.ltde_hr]
        tde_sv_interval = [self._hemodynamic_variables[4] + self._pd_hemo.base_sv * self._pd_hemo.ltde_sv]

        # Solving the nonlinear state space model
        for i in range(len(t) - 1):
            t_span = [t[i], t[i + 1]]
            initial_conditions = self._hemodynamic_variables
            t_eval = np.linspace(t_span[0], t_span[1], 41)
            options = {'rtol': 1e-2, 'atol': 1e-4}
            sol = solve_ivp(self._pd_hemo.ode_prop_hemodynamic, t_span, initial_conditions,
                            args=(cp_prop[i], cp_remi[i]), t_eval=t_eval, **options)
            y_ode = sol.y[:, -1]  # This is correct, sol.y is a ndarray of shape (5, 41)

            tpr_interval.append(y_ode[0])
            sv_star_interval.append(y_ode[1])
            hr_star_interval.append(y_ode[2])
            tde_hr = y_ode[3] + self._pd_hemo.base_hr * self._pd_hemo.ltde_hr
            tde_sv = y_ode[4] + self._pd_hemo.base_sv * self._pd_hemo.ltde_sv
            tde_hr_interval.append(tde_hr)
            tde_sv_interval.append(tde_sv)
            self._hemodynamic_variables = y_ode[:5]

        hr_interval = [x + y for x, y in zip(hr_star_interval, tde_hr_interval)]
        hr_interval = np.array(hr_interval)
        sv_star_interval = np.array(sv_star_interval)
        tde_sv_interval = np.array(tde_sv_interval)

        sv_interval = (sv_star_interval + tde_sv_interval) * (
                1 - self._pd_hemo.hr_sv * np.log(hr_interval / (self._pd_hemo.base_hr * (1 + self._pd_hemo.ltde_hr))))

        # Adding the effect of disturbances
        if self._disturbance_model is not None:
            interval = min(len(t), len(self._hr_dis) - self._time)
            tpr_interval = (tpr_interval * hr_interval + self._map_dis[
                                                         self._time: self._time + interval] / sv_interval) / (
                                   hr_interval + self._hr_dis[self._time: self._time + interval])
            hr_interval = hr_interval + self._hr_dis[self._time: self._time + interval]

        co_interval = hr_interval * sv_interval / 1000 + co_nore
        map_interval = tpr_interval * co_interval * 1000 + map_nore

        # Add the effect of blood volume status on the stroke volume
        if self._volume_status is not None and self._time in self._volume_status.keys():
            self._volume_status_coeff = self._volume_status[self._time].value

        hr_interval *= self._volume_status_coeff['hr']
        sv_interval *= self._volume_status_coeff['sv']

        # Linear interaction with norepinephrine
        map_interval *= self._volume_status_coeff['map']
        co_interval *= self._volume_status_coeff['co']

        return co_interval, map_interval, hr_interval, sv_interval, tpr_interval

    def _compute_nmb(self, cp_rocu: np.ndarray, t: np.ndarray):
        """
        Compute the neuromuscular blockade using the pharmacodynamic model.
        :param cp_rocu: The plasma concentration of rocuronium.
        :type cp_rocu: np.ndarray
        :param t: The time vector.
        :type t: np.ndarray
        """
        # Compute the effect-site concentration of Rocuronium
        _, ce_rocu_sim, x_ce_rocu_sim = signal.lsim(self._pd_nmb.pd_ce, cp_rocu, t, X0=self._x_ce_rocu)

        ce_rocu_sim[ce_rocu_sim < 0] = 0
        self._x_ce_rocu = x_ce_rocu_sim[-1]  # This is correct, x_ce_rocu is a 1D array

        # Compute the probability for each category of NMB
        nmb = self._pd_nmb.hillfun(ce_rocu_sim)

        return nmb, ce_rocu_sim

    def step(self, u_prop: float, u_remi: float, u_nore: float, u_rocu: float, t_s: int = 5):
        """
        Simulate one step of the patient.
        :param u_prop: The propofol input. [mg * sec^-1]
        :param u_remi: The remifentanil input. [µg * sec^-1]
        :param u_nore: The norepinephrine input. [nmol * sec^-1]
        :param u_rocu: The rocuronium input. [mg * sec^-1]
        :param t_s: The sampling time in seconds. Default is 5.
        """

        t = np.arange(0, t_s, 1)  # time vector that goes from 0 to t_s seconds with a step of 1 second

        # Generate constant input arrays for the simulation that lasts t_s seconds
        # One step of the simulation is t_s seconds
        u_sim_prop = np.full(t.shape, u_prop)
        u_sim_remi = np.full(t.shape, u_remi)
        u_sim_nore = np.full(t.shape, u_nore)
        u_sim_rocu = np.full(t.shape, u_rocu)

        self._u_remi.extend(u_sim_remi)
        self._u_prop.extend(u_sim_prop)
        self._u_nore.extend(u_sim_nore)
        self._u_rocu.extend(u_sim_rocu)

        # Compute plasma concentrations
        cp_prop, cp_remi, c_nore, cp_rocu = self._compute_plasma_concentration(u_sim_prop, u_sim_remi,
                                                                               u_sim_nore, u_sim_rocu, t)

        # Compute the disturbances caused by stimulations
        if self._disturbance_model is not None:
            self._doh_dis, self._hr_dis, self._map_dis = self._disturbance_model.get_disturbances(self._time,
                                                                                                  cp_prop[-t_s],
                                                                                                  cp_remi[-t_s])

        # Compute effect-site concentrations of propofol and remifentanil
        _, ce_prop_sim, x_ce_prop_sim = signal.lsim(self._pd_doh.pd_prop_ce, cp_prop, t, X0=self._x_ce_prop)
        _, ce_remi_sim, x_ce_remi_sim = signal.lsim(self._pd_doh.pd_remi_ce, cp_remi, t, X0=self._x_ce_remi)

        ce_prop_sim[ce_prop_sim < 0] = 0
        ce_remi_sim[ce_remi_sim < 0] = 0

        # Compute depth of hypnosis
        if self._pd_doh._delay == 0:
            ce_delayed_prop = ce_prop_sim
            x_ce_delayed_prop = x_ce_prop_sim
        else:
            _, ce_delayed_prop, x_ce_delayed_prop = signal.lsim(self._pd_doh.get_delay_ss(), ce_prop_sim, t,
                                                                X0=self._x_delay)

        ce_delayed_prop[ce_delayed_prop < 0] = 0

        self._x_ce_prop = x_ce_prop_sim[-1]
        self._x_ce_remi = x_ce_remi_sim[-1]
        self._x_delay = x_ce_delayed_prop[-1]

        self._ce_prop = np.append(self._ce_prop, ce_prop_sim[:])
        self._ce_remi = np.append(self._ce_remi, ce_remi_sim[:])
        self._ce_del = np.append(self._ce_del, ce_delayed_prop[:])

        # Based on the monitor used for measuring depth of hypnois, wav or bis index is updated.
        if self._doh_measure == DoHMeasure.WAV or self._doh_measure == DoHMeasure.BOTH:
            wav_interval = np.array(self._compute_WAV(ce_delayed_prop, ce_remi_sim, t), dtype=np.float64)
            if self._disturbance_model is not None:
                interval = min(t_s, len(self._doh_dis) - self._time)
                wav_interval += self._doh_dis[self._time: self._time + interval]
                wav_interval = np.clip(wav_interval, 0, 100)
        else:
            wav_interval = np.zeros(len(t))

        self._WAV = np.append(self._WAV, wav_interval[:])

        if self._doh_measure == DoHMeasure.BIS or self._doh_measure == DoHMeasure.BOTH:
            bis_interval = np.array(self._compute_BIS(ce_delayed_prop, ce_remi_sim, t), dtype=np.float64)
            if self._disturbance_model is not None:
                dbis, _ = PharmacodynamicDoH.bis_sensor_dynamics()
                self._doh_dis = np.asarray(self._doh_dis)
                time = np.linspace(0, len(self._doh_dis) - 1, len(self._doh_dis))
                _, delayed_disturb, _ = signal.lsim(dbis, self._doh_dis, time)
                interval = min(t_s, len(delayed_disturb) - self._time)
                bis_interval += delayed_disturb[self._time: self._time + interval]
                bis_interval = np.clip(bis_interval, 0, 100)
        else:
            bis_interval = np.zeros(len(t))

        self._BIS = np.append(self._BIS, bis_interval[:])

        elems = wav_interval if self._doh_measure == DoHMeasure.WAV or self._doh_measure == DoHMeasure.BOTH else bis_interval

        # Check if the patient is in the maintenance phase: if the patient's wav of bis has been under 60  for 3 minutes
        if self._patient_phase == PatientPhase.INDUCTION:
            for elem in elems:
                self._maintenance_time = self._maintenance_time + 1 if elem < 60 else 0
                if self._maintenance_time > 180:
                    self._patient_phase = PatientPhase.MAINTENANCE
                    break

        # Check if the patient is in steady state: if the output is stable for 3 minutes around the target value (WAV or BIS = 50)
        # Start checking after 1000 seconds (visually estimated time for the patient to reach the target value) to improve the performance
        if self._patient_phase == PatientPhase.MAINTENANCE and self._time > 1000 and not self._steady_state:
            lower_bound = 47.5  # 50 - 5%
            upper_bound = 52.5  # 50 + 5%

            for elem in elems:
                # Check if the element is within the target range and close to the previous value
                if lower_bound < elem < upper_bound and (
                        not self._steady_state_values or abs(elem - self._steady_state_values[-1]) < 0.2):
                    self._steady_state_values.append(elem)
                else:
                    self._steady_state_values.clear()  # Reset if an outlier is detected

                # Check if the steady-state condition is met
                if len(self._steady_state_values) == 180 and abs(
                        self._steady_state_values[0] - self._steady_state_values[-1]) < 0.5:
                    self._steady_state = True
                    break

        # Compute the hemodynamic variables
        co_interval, map_interval, hr_interval, sv_interval, tpr_interval = self._compute_hemodynamic_variables(cp_prop,
                                                                                                                cp_remi,
                                                                                                                c_nore,
                                                                                                                t)
        # Compute the effect-site concentration of rocuronium and the probabilities of the NMB categories
        nmb, ce_rocu_sim = self._compute_nmb(cp_rocu, t)

        # Save the stats in the recording arrays
        self._CO = np.append(self._CO, co_interval[:])
        self._MAP = np.append(self._MAP, map_interval[:])
        self._HR = np.append(self._HR, hr_interval[:])
        self._SV = np.append(self._SV, sv_interval[:])
        self._TPR = np.append(self._TPR, tpr_interval[:])

        self._ce_rocu = np.append(self._ce_rocu, ce_rocu_sim[:])
        self._NMB_m0 = np.append(self._NMB_m0, nmb[0, :])
        self._NMB_m1 = np.append(self._NMB_m1, nmb[1, :])
        self._NMB_m2 = np.append(self._NMB_m2, nmb[2, :])
        self._NMB_m3 = np.append(self._NMB_m3, nmb[3, :])

        self._time += t_s
