import logging
import pandas as pd
import os

from .disturbance import Disturbance
from .patient import Patient
from .TCI import TCI
from .utils.pid import PID
from .utils.plot_simulation import plot_simulation
from .utils.read_file import read_data
from .utils.enums import Model, DoHMeasure, Interaction, TciMode, Drug, SimulatorMode, BloodSampling

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


class Simulator:
    def __init__(self):
        # List of the patients simulated up to now
        self._patients = []

        # List of lists to store patients' propofol and remifentanil infusion rates
        self._u_prop_all = []
        self._u_remi_all = []
        self._u_nore_all = []
        self._u_rocu_all = []

        # List of lists to store patients' propofol, remifentanil, norepinephrine, and rocuronium plasma concentrations
        # and effect site concentrations of propofol and remifentanil
        self._cp_prop_all = []
        self._ce_prop_all = []
        self._cp_remi_all = []
        self._ce_remi_all = []
        self._c_nore_all = []
        self._cp_rocu_all = []

        # Delayed version of ce and ce incorporating sensor dynamics (WAV and BIS)
        self._ce_del_all = []
        self._ce_wav_all = []
        self._ce_bis_all = []

        # List of lists to store the WAV and BIS values of patients
        self._WAV_all = []
        self._BIS_all = []

        # List of lists to store the hemodynamics variables of patients
        self._MAP_all = []
        self._CO_all = []
        self._HR_all = []
        self._SV_all = []
        self._TPR_all = []

        # List of lists to store the neuromuscular blockade values of patients
        self._NMB_m0_all = []
        self._NMB_m1_all = []
        self._NMB_m2_all = []
        self._NMB_m3_all = []

        # List for the current patient
        self._current_patient = None  # The current patient being simulated
        self._t_sim = None  # The simulation time in seconds for the current patient, should be divisible by t_s
        self._t_s = None  # The sampling time in seconds for the simulation of the current patient
        self._current_time = 0

        self.nore_molweight = 169.18  # Molecular weight of norepinephrine in [g/mol]

    @classmethod
    # If the simulator mode is not specified, the default mode is set to INFUSION
    def create(cls, mode: SimulatorMode = SimulatorMode.INFUSION):
        if mode == SimulatorMode.INFUSION:
            return cls()
        elif mode == SimulatorMode.CONCENTRATION:
            return SimulatorConcentration()
        else:
            raise ValueError("Invalid simulator mode")

    def reset(self):
        """
        Reset all the values of the simulator.
        """
        default_values = {
            '_patients': [],
            '_u_prop_all': [],
            '_u_remi_all': [],
            '_u_nore_all': [],
            '_u_rocu_all': [],
            '_cp_prop_all': [],
            '_cp_remi_all': [],
            '_c_nore_all': [],
            '_cp_rocu_all': [],
            '_ce_prop_all': [],
            '_ce_remi_all': [],
            '_ce_del_all': [],
            '_ce_wav_all': [],
            '_ce_bis_all': [],
            '_WAV_all': [],
            '_BIS_all': [],
            '_MAP_all': [],
            '_CO_all': [],
            '_HR_all': [],
            '_SV_all': [],
            '_TPR_all': [],
            '_NMB_m0_all': [],
            '_NMB_m1_all': [],
            '_NMB_m2_all': [],
            '_NMB_m3_all': [],
            '_current_patient': None,
            '_t_sim': None,
            '_t_s': None,
            '_current_time': 0,
        }

        for attr, value in default_values.items():
            setattr(self, attr, value)

    def get_patient_demographics(self):
        """
        :returns: The demographic parameters for the current patient: [age, height, weight, gender, bmi, lbm]
        :rtype: dict
        """
        return self._current_patient.get_patient_demographics()

    def get_patient_state(self):
        """
        :returns: The state of the current patient at the current time:
        [WAV, BIS, MAP, CO, HR, SV,NMB_m0, NMB_m1, NMB_m2, NMB_m3, cp_prop, cp_remi, c_nore, cp_rocu, ce_prop, ce_remi, ce_del, ce_wav, ce_bis]
        :rtype: dict
        """
        state = self._current_patient.get_patient_state()
        c_nore = [x * self.nore_molweight / 1000 for x in state['c_nore']] if isinstance(state['c_nore'], list) else \
            state['c_nore'] * self.nore_molweight / 1000
        state['c_nore'] = c_nore
        return state

    def get_patient_state_history(self):
        """
        :returns: The state history from the start of the simulation to the current time for the current patient:
        [WAV, BIS, MAP, CO, HR, SV,NMB_m0, NMB_m1, NMB_m2, NMB_m3, cp_prop, cp_remi, c_nore, cp_rocu, ce_prop, ce_remi, ce_del, ce_wav, ce_bis]
        :rtype: dict
        """
        state = self._current_patient.get_patient_state_history()
        c_nore = [x * self.nore_molweight / 1000 for x in state['c_nore']] if isinstance(state['c_nore'], list) else \
            state['c_nore'] * self.nore_molweight / 1000
        state['c_nore'] = c_nore
        return state

    def get_patient_input_history(self):
        """
        :returns: The patient input history from the start of the simulation to the current time for the current patient: [u_prop, u_remi, u_nore, u_rocu] in [mg * sec^-1], [µg * sec^-1], [µg * sec^-1], [µg * sec^-1]
        :rtype: dict
        """
        input_history = self._current_patient.get_patient_input_history()
        # Convert the norepinephrine infusion rates from [nmol * sec^-1] to [ng * sec^-1]
        u_nore = [x * self.nore_molweight / 1000 for x in input_history['u_nore']] if isinstance(
            input_history['u_nore'], list) else input_history['u_nore'] * self.nore_molweight / 1000
        input_history['u_nore'] = u_nore
        return input_history

    def get_patient_phase(self):
        """
        :returns: The phase of the current patient at the current time: [induction, maintenance, emergence]
        :rtype: PatientPhase
        """
        return self._current_patient.get_patient_phase()[0]

    def get_patient_results(self, num_simulation: int = None):
        """
        Get the patient records for the given index number of simulation.
        If no number is given, return all the patient records.
        :param num_simulation: The index of the simulation to get the data from. If None, return all the simulations.
        :type num_simulation: int, optional
        :returns: A dictionary of the requested variables or all patient data.
        :rtype: dict
        """
        if num_simulation is not None:
            if not 0 <= num_simulation < len(self._patients):
                raise IndexError("Simulation index out of range")
            # Get all the attributes ending with '_all' and that are lists and return the data for the given simulation
            data = {var: getattr(self, var)[num_simulation] for var in dir(self)
                    if var.endswith('_all') and isinstance(getattr(self, var), list)}
        else:
            data = {var: getattr(self, var) for var in dir(self)
                    if var.endswith('_all') and isinstance(getattr(self, var), list)}
        return data

    def _initialize_simulation(self, patient_data, t_sim, t_s=5,
                               opiates=None, blood_sampling=None,
                               internal_states=None,
                               output_init=None,
                               pk_models=None,
                               pd_models=None,
                               interaction=None,
                               doh_measure=None,
                               stimuli=None,
                               volume_status=None,
                               bring_to_maintenance=False,
                               seed=None):
        """
        Initialize the simulation with given patient data and simulation parameters.
        :param patient_data: The patient data [height, weight, age, gender, E0, k_d , delay, C50P, gamma]
        :type patient_data: list
        :param t_sim: The simulation time in seconds
        :type t_sim: int
        :param t_s: The sampling time in seconds. It should be a divisor of t_sim and greater than 1 second.
        :type t_s: int
        :param opiates: True if opiates are used, False otherwise
        :type opiates: bool
        :param blood_sampling: The type of blood sampling to use for the simulation (arterial or venous)
        :type blood_sampling: BloodSampling
        :param internal_states : The internal states of the patient
        :type internal_states: dict
        :param output_init: The initial values for map, hr, sv
        :type output_init: dict
        :param pk_models: The pharmacokinetic models for each drug
        :type pk_models: dict
        :param pd_models: The pharmacodynamic models for propofol and remifentanil
        :type pd_models: dict
        :param interaction: The interaction model to use for the simulation (surface or no_interaction)
        :type interaction: Interaction
        :param doh_measure: The depth of hypnosis measure to compute for the simulation (BIS, WAV or both)
        :type doh_measure: DoHMeasure
        :param stimuli: The list of stimuli for the simulation
        :type stimuli: dict
        :param volume_status: The list of volume status for the simulation
        :type volume_status: dict
        :param bring_to_maintenance: True if the patient should be brought to a maintenance state with PID control
        :type bring_to_maintenance: bool
        :param seed: seed for the variability of the patient hemodynamics
        :type seed: int
        """

        # Setup patient
        patient_data = self._refactor_patient_data(patient_data)

        # Initialize patient
        self._current_patient = Patient(
            patient_data,
            internal_states=internal_states,
            output_init=output_init,
            pk_models=pk_models,
            pd_models=pd_models,
            opiates=opiates,
            blood_sampling=blood_sampling,
            interaction=interaction,
            doh_measure=doh_measure,
            disturbance_model=None,
            volume_status=None,
            seed=seed
        )

        # Set simulation parameters
        self._patients.append(self._current_patient)
        self._t_sim = t_sim
        self._t_s = t_s
        self._current_time = 0

        if bring_to_maintenance:
            self.stimuli_starts = None
            self.volume_status_starts = None
            self._bring_patient_to_maintenance(t_s, doh_measure)

        disturbance_model = Disturbance(t_sim, stimuli, bring_to_maintenance) if stimuli is not None else None
        self._current_patient.set_disturbance_model(disturbance_model)
        self._current_patient.set_volume_status(volume_status)

        self.stimuli_starts = stimuli.keys() if stimuli is not None else None
        self.volume_status_starts = volume_status.keys() if volume_status is not None else None

    def _bring_patient_to_maintenance(self, t_s, doh_measure):
        """
        Helper function to bring the patient to a maintenance state with PID control.
        Note: The PID has been optimized for the 44 patients in the dataset using the PK-PD model of Eleveld et al. for both propofol and remifentanil.
        """
        # Initialize PID parameters
        kp, ki, kd = 1.187610442449273, 0.001215206797108504, 0.20413144686025936
        pid = PID(kp=kp, ki=ki, kd=kd, setpoint=0.5, sample_time=t_s)

        measure = 'WAV' if doh_measure in [DoHMeasure.WAV, DoHMeasure.Both] else 'BIS'

        # While the patient does not reach  a steady state
        while not self._current_patient.get_patient_phase()[1]:
            doh = self.get_patient_state()[measure]
            u_prop = pid.compute(doh / 100, t_s)
            self.one_step_simulation(u_prop, u_prop / 2, 0, 0)
            if self._current_time > 5000:
                raise ValueError("The patient could not reach a steady state. Please check the PID parameters.")

        # Reset for the maintenance phase
        self._current_patient.clear_data()  # Clear the data except the last state
        self._current_patient.set_time(0)  # Reset the time of the patient to 0
        self._current_time = 0

    def init_simulation_from_file(self, id_patient, t_sim, t_s=5,
                                  opiates=None, blood_sampling=None,
                                  internal_states=None,
                                  output_init=None,
                                  pk_models=None,
                                  pd_models=None,
                                  interaction=None,
                                  doh_measure=None,
                                  stimuli=None,
                                  volume_status=None,
                                  bring_to_maintenance=False,
                                  seed=None):
        """
        Initialize the simulation using patient data from the file in the 'parameters' folder.
        :param id_patient: The index of the patient to simulate
        :type id_patient: int
        :param t_sim: The simulation time in seconds
        :type t_sim: int
        :param t_s: The sampling time in seconds
        :type t_s: int
        :param opiates: True if opiates are used, False otherwise
        :type opiates: bool
        :param blood_sampling: The type of blood sampling to use for the simulation (arterial or venous)
        :type blood_sampling: BloodSampling
        :param internal_states : The internal states of the patient
        :type internal_states: dict
        :param output_init: The initial values for map, hr, sv
        :type output_init: dict
        :param pk_models: The pharmacokinetic models for each drug
        :type pk_models: dict
        :param pd_models: The pharmacodynamic models for propofol and remifentanil
        :type pd_models: dict
        :param interaction: The interaction model to use for the simulation (surface or no_interaction)
        :type interaction: Interaction
        :param doh_measure: The depth of hypnosis measure to use for the simulation (BIS, WAV or both)
        :type doh_measure: DoHMeasure
        :param stimuli: The list of stimuli for the simulation
        :type stimuli: dict
        :param volume_status: The list of volume status for the simulation
        :type volume_status: dict
        :param bring_to_maintenance: True if the patient should be brought to a maintenance state with PID control
        :type bring_to_maintenance: bool
        :param seed: seed for the variability of the patient hemodynamics
        :type seed: int
        """
        if id_patient < 0 or id_patient > 43:
            raise ValueError("Patient index out of range")

        pk_model_prop = pk_models.get('prop', Model.ELEVELD)
        file_name_patient = 'pd_data_schnider.csv' if pk_model_prop == Model.SCHNIDER else 'pd_data_eleveld.csv'

        current_path = os.path.dirname(os.path.abspath(__file__))
        path_patient = os.path.join(current_path, 'parameters')
        patient_data = read_data(file_name_patient, path_patient, row_index=id_patient)

        self._initialize_simulation(patient_data, t_sim, t_s,
                                    opiates, blood_sampling,
                                    internal_states,
                                    output_init,
                                    pk_models,
                                    pd_models,
                                    interaction,
                                    doh_measure,
                                    stimuli,
                                    volume_status,
                                    bring_to_maintenance,
                                    seed)

    def init_simulation_from_data(self, data, t_sim, t_s=5,
                                  opiates=None, blood_sampling=None,
                                  internal_states=None,
                                  output_init=None,
                                  pk_models=None,
                                  pd_models=None,
                                  interaction=None,
                                  doh_measure=None,
                                  stimuli=None,
                                  volume_status=None,
                                  bring_to_maintenance=False,
                                  seed=None):
        """
        Initialize the simulation with patient parameters directly.
        :param data: The patient data [height, weight, age, gender, bmi, lbm, E0, k_d , delay, C50P, gammaP, rms_nonlin]
        Note: The parameters [E0, k_d , delay, C50P, gammaP, rms_nonlin] need to be provided only the for PATIENT_SPECIFIC PD model for propofol.
        :type data: list
        :param t_sim: The simulation time in seconds
        :type t_sim: int
        :param t_s: The sampling time in seconds. It should be a divisor of t_sim and greater than 1 second
        :type t_s: int
        :param opiates: True if opiates are used, False otherwise
        :type opiates: bool
        :param blood_sampling: The type of blood sampling to use for the simulation (arterial or venous)
        :type blood_sampling: BloodSampling
        :param internal_states : The internal states of the patient
        :type internal_states: dict
        :param output_init: The initial values for map, hr, sv
        :type output_init: dict
        :param pk_models: The pharmacokinetic models for propofol and remifentanil
        :type pk_models: dict
        :param pd_models: The pharmacodynamic models for propofol and remifentanil
        :type pd_models: dict
        :param interaction: The interaction model to use for the simulation (surface or no_interaction)
        :type interaction: Interaction
        :param doh_measure: The depth of hypnosis measure to use for the simulation (BIS, WAV or both)
        :type doh_measure: DoHMeasure
        :param stimuli: The list of stimuli for the simulation
        :type stimuli: dict
        :param volume_status: The list of volume status for the simulation
        :type volume_status: dict
        :param bring_to_maintenance: True if the patient should be brought to a maintenance state with PID control
        :type bring_to_maintenance: bool
        :param seed: seed for the variability of the patient hemodynamics
        :type seed: int
        """

        if data is None:
            raise ValueError("No patient data found. Please add the data for a simulation first.")

        self._initialize_simulation(data, t_sim, t_s, opiates, blood_sampling,
                                    internal_states,
                                    output_init,
                                    pk_models,
                                    pd_models,
                                    interaction,
                                    doh_measure,
                                    stimuli,
                                    volume_status,
                                    bring_to_maintenance,
                                    seed)

    def run_complete_simulation(self, u_prop: list, u_remi: list, u_nore: list, u_rocu: list):
        """
        Simulates the patient for the given inputs from the start to the end of the simulation.
        :param u_prop: A list of propofol inputs [mg * sec^-1]
        :type u_prop: list
        :param u_remi: A list of remifentanil inputs [µg * sec^-1]
        :type u_remi: list
        :param u_nore: A list of norepinephrine inputs [µg * sec^-1]
        :type u_nore: list
        :param u_rocu: A list of rocuronium inputs [µg * sec^-1]
        :type u_rocu: list
        """
        self._check_run_inputs(u_prop, u_remi, u_nore, u_rocu)

        for time in range(self._t_sim // self._t_s):
            self.one_step_simulation(u_prop[time * self._t_s], u_remi[time * self._t_s], u_nore[time * self._t_s],
                                     u_rocu[time * self._t_s])

    def one_step_simulation(self, u_prop: float, u_remi: float, u_nore: float, u_rocu: float):
        """
        Simulates the patient for a single step with the given inputs.
        :param u_prop: The propofol input [mg * sec^-1]
        :type u_prop: float
        :param u_remi: The remifentanil input [µg * sec^-1]
        :type u_remi: float
        :param u_nore: The norepinephrine input [µg * sec^-1]
        :type u_nore: float
        :param u_rocu: The rocuronium input [µg * sec^-1]
        :type u_rocu: float
        """

        if self._current_patient is None:
            raise ValueError("No patient data found. Please add a simulation first.")

        u_nore = u_nore * 1000 / self.nore_molweight  # Norepinephrine in [nmol/s]

        self._current_patient.step(u_prop, u_remi, u_nore, u_rocu, self._t_s)

        self._current_time += self._t_s

    def is_simulation_finished(self):
        """
        Check if the simulation is finished.

        :returns: True if the simulation is finished, False otherwise
        :rtype: bool
        """
        return self._current_time >= self._t_sim

    def save_simulation(self):
        """
        Save the results of the current simulation to the list of all patients' data.
        """
        if self._current_patient is None:
            raise ValueError("No patient data found. Please add a simulation first.")

        inputs = self.get_patient_input_history()
        state = self.get_patient_state_history()
        self._add_result_to_all(inputs, state)

    def _add_result_to_all(self, inputs, state):
        """
        Add the current patient data to the list of all patients' data.
        """

        # Add the input data
        self._u_prop_all.append(inputs['u_prop'])
        self._u_remi_all.append(inputs['u_remi'])
        self._u_nore_all.append(inputs['u_nore'])
        self._u_rocu_all.append(inputs['u_rocu'])

        # Add the doh
        self._WAV_all.append(state['WAV'])
        self._BIS_all.append(state['BIS'])

        # Add the Hemodynamic variables
        self._MAP_all.append(state['MAP'])
        self._CO_all.append(state['CO'])
        self._HR_all.append(state['HR'])
        self._SV_all.append(state['SV'])
        self._TPR_all.append(state['TPR'])

        # Add the neuromuscular blockade
        self._NMB_m0_all.append(state['NMB_m0'])
        self._NMB_m1_all.append(state['NMB_m1'])
        self._NMB_m2_all.append(state['NMB_m2'])
        self._NMB_m3_all.append(state['NMB_m3'])

        # Add the plasma and effect site concentrations
        self._cp_prop_all.append(state['cp_prop'])
        self._cp_remi_all.append(state['cp_remi'])
        self._c_nore_all.append(state['c_nore'])
        self._cp_rocu_all.append(state['cp_rocu'])
        self._ce_prop_all.append(state['ce_prop'])
        self._ce_remi_all.append(state['ce_remi'])
        self._ce_del_all.append(state['ce_del'])
        self._ce_wav_all.append(state['ce_wav'])
        self._ce_bis_all.append(state['ce_bis'])

    def _refactor_patient_data(self, data):
        """
        Refactor the patient data to include the BMI and the lean body mass.
        :param data: The patient data
        :type data: list
        :return: The refactored patient data with the BMI and the lean body mass
        :rtype: list
        """
        # Calculate the BMI and the lean body mass
        height = data[0]
        weight = data[1]
        gender = data[3]  # 0 for male and 1 for female
        bmi = weight / (0.01 * height) ** 2

        # Fat Free Mass (FFM) calculation using the Janmahasatian formula
        # Janmahasatian S, Duffull SB, Ash S, Ward LC, Byrne NM, Green B.
        # Quantification of lean bodyweight. Clin Pharmacokinet. 2005;44(10):1051-65.
        # doi: 10.2165/00003088-200544100-00004. PMID: 16176118.
        num = 9.27 * 1000 * weight
        den = 6.68 * 1000 + 216 * bmi if gender == 0 else 8.78 * 1000 + 244 * bmi
        lbm = num / den
        patient_data = data[:4] + [bmi, lbm] + data[4:]
        return patient_data

    def _check_run_inputs(self, input_prop, input_remi, input_nore, input_rocu):
        """
        Check the inputs for the simulation and adjust the simulation time if necessary.
        """
        if self._current_patient is None:
            raise ValueError("No patient data found. Please add the data for a simulation first.")

        # Check if the input lists have the same length as the simulation time
        # If not, adjust the simulation time
        if len(input_prop) != self._t_sim or len(input_remi) != self._t_sim or len(input_nore) != self._t_sim or len(
                input_rocu) != self._t_sim:
            self._t_sim = min(len(input_prop), len(input_remi), len(input_nore), len(input_rocu), self._t_sim)

    def plot_simulation(self, num_simulation: int = None, titles: list = None, labels=None, use_cmap=False,
                        save_path=None):
        """
        Create the plots for the simulation results for the given simulation index. If no index is given, plot all the
        simulations.
        :params num_simulation: The index of the simulation to plot. If None, plot all the simulations.
        :type num_simulation: int
        :param titles: List of titles for the plots
        :param labels: List of labels for the plots
        :param use_cmap: Boolean to use colormap
        :param save_path: Path to save the plots
        :return: Figures of the simulation results [DoH, Hemodynamics, NMB]
        :rtype: tuple(matlplotlib.figure.Figure)
        """
        if num_simulation is None:
            # Pass the entire list of lists to plot_simulation_temp
            fig1, fig2, fig3 = plot_simulation(self._u_prop_all, self._u_remi_all,
                                               self._u_nore_all, self._u_rocu_all,
                                               self._cp_prop_all, self._cp_remi_all,
                                               self._c_nore_all, self._cp_rocu_all,
                                               self._ce_prop_all, self._ce_remi_all,
                                               self._WAV_all, self._BIS_all,
                                               self._MAP_all, self._CO_all,
                                               self._HR_all, self._SV_all, self._TPR_all,
                                               self._NMB_m0_all, self._NMB_m1_all,
                                               self._NMB_m2_all, self._NMB_m3_all,
                                               titles, labels,
                                               self.stimuli_starts, self.volume_status_starts,
                                               use_cmap, save_path)
        else:
            # Plot for a specific patient
            fig1, fig2, fig3 = plot_simulation(self._u_prop_all[num_simulation],
                                               self._u_remi_all[num_simulation],
                                               self._u_nore_all[num_simulation],
                                               self._u_rocu_all[num_simulation],
                                               self._cp_prop_all[num_simulation],
                                               self._cp_remi_all[num_simulation],
                                               self._c_nore_all[num_simulation],
                                               self._cp_rocu_all[num_simulation],
                                               self._ce_prop_all[num_simulation],
                                               self._ce_remi_all[num_simulation],
                                               self._WAV_all[num_simulation],
                                               self._BIS_all[num_simulation],
                                               self._MAP_all[num_simulation],
                                               self._CO_all[num_simulation],
                                               self._HR_all[num_simulation],
                                               self._SV_all[num_simulation],
                                               self._TPR_all[num_simulation],
                                               self._NMB_m0_all[num_simulation],
                                               self._NMB_m1_all[num_simulation],
                                               self._NMB_m2_all[num_simulation],
                                               self._NMB_m3_all[num_simulation],
                                               titles, labels,
                                               self.stimuli_starts, self.volume_status_starts,
                                               use_cmap, save_path)

        return fig1, fig2, fig3

    def save_to_csv(self, path, filename):
        """
        Save the results of the simulation to a CSV file.

        :param path: The directory path where the CSV file will be saved.
        :type path: str
        :param filename: The name of the CSV file.
        :type filename: str
        """

        # Generate the full path for the CSV file
        full_path = os.path.join(path, filename)
        # Ensure the directory exists
        os.makedirs(os.path.dirname(full_path), exist_ok=True)

        # Retrieve simulation data
        data = self.get_patient_results()

        # Define new column names
        columns_mapping = {
            '_u_prop_all': 'u_prop', '_u_remi_all': 'u_remi', '_u_nore_all': 'u_nore', '_u_rocu_all': 'u_rocu',
            '_cp_prop_all': 'cp_prop', '_cp_remi_all': 'cp_remi', '_c_nore_all': 'c_nore', '_cp_rocu_all': 'cp_rocu',
            '_ce_prop_all': 'ce_prop', '_ce_remi_all': 'ce_remi', '_ce_del_all': 'ce_del', '_ce_wav_all': 'ce_wav',
            '_ce_bis_all': 'ce_bis', '_WAV_all': 'wav', '_BIS_all': 'bis', '_MAP_all': 'map', '_CO_all': 'co',
            '_HR_all': 'hr', '_SV_all': 'sv', '_NMB_m0_all': 'nmb_m0', '_NMB_m1_all': 'nmb_m1', '_NMB_m2_all': 'nmb_m2',
            '_NMB_m3_all': 'nmb_m3'
        }

        # Prepare the data for the DataFrame
        df_data = {}

        for old_col_name, new_col_name in columns_mapping.items():
            if old_col_name in data:
                # Add each array as a column in the DataFrame, with indices for each time step
                if len(data[old_col_name]) > 1:
                    for i in range(len(data[old_col_name])):
                        df_data[f'{new_col_name}_{i + 1}'] = data[old_col_name][i]
                else:
                    df_data[new_col_name] = data[old_col_name][0]

        # Create the DataFrame
        df = pd.DataFrame(df_data.values(), index=df_data.keys()).T

        # Save the DataFrame to CSV
        df.to_csv(full_path, index=False)


class SimulatorConcentration(Simulator):
    def __init__(self):
        super().__init__()
        self.tci_prop = None
        self.tci_remi = None
        self.tci_nore = None
        self.tci_rocu = None

    def reset(self):
        """
        Reset all the values of the AReS.
        """
        super().reset()
        self.tci_prop = None
        self.tci_remi = None
        self.tci_nore = None
        self.tci_rocu = None

    def _initialize_simulation(self, patient_data, t_sim, t_s=5,
                               opiates=None, blood_sampling=None,
                               internal_states=None,
                               output_init=None,
                               pk_models=None,
                               pd_models=None,
                               interaction=None,
                               doh_measure=None,
                               stimuli=None,
                               volume_status=None,
                               bring_to_maintenance=False,
                               seed=None,
                               **kwargs):

        """
        Initialize the simulation with given patient data and simulation parameters.
        :param patient_data: The patient data [height, weight, age, gender, E0, k_d , delay, C50P, gammaP, rms_nonlin]ù
        :type patient_data: list
        :param t_sim: The simulation time in seconds
        :type t_sim: int
        :param t_s: The sampling time in seconds
        :type t_s: int
        :param opiates: True if opiates are used, False otherwise
        :type opiates: bool
        :param blood_sampling: The type of blood sampling to use for the simulation (arterial or venous)
        :type blood_sampling: BloodSampling
        :param internal_states : The internal states of the patient
        :type internal_states: dict
        :param output_init: The initial values for map, hr, sv
        :type output_init: dict
        :param pk_models: The pharmacokinetic models for propofol and remifentanil
        :type pk_models: dict
        :param pd_models: The pharmacodynamic models for propofol and remifentanil
        :type pd_models: dict
        :param interaction: The interaction model to use for the simulation (surface or no_interaction)
        :type interaction: Interaction
        :param doh_measure: The depth of hypnosis measure to use for the simulation (BIS, WAV or both)
        :type doh_measure: DoHMeasure
        :param stimuli: The list of stimuli for the simulation
        :type stimuli: dict
         :param volume_status: The list of volume status for the simulation
        :type volume_status: dict
        :param bring_to_maintenance: True if the patient should be brought to a maintenance state with PID control
        :type bring_to_maintenance: bool
        :param seed: seed for the variability of the patient hemodynamics
        :type seed: int
        :param kwargs: Additional keyword arguments to pass to the TCI initialization
        :type kwargs: dict
        """

        if not 'data_TCI' in kwargs:
            raise ValueError("No TCI data found. Please add the TCI data for a simulation first.")
        data_TCI = kwargs.get('data_TCI', [])
        data_TCI = self._refactor_patient_data(data_TCI)

        pk_models_TCI = kwargs.get('pk_models_TCI', {})
        pd_models_TCI = kwargs.get('pd_models_TCI', {})

        opiates = kwargs.get('opiates', True)
        blood_sampling = kwargs.get('blood_sampling', BloodSampling.ARTERIAL)

        limits_TCI = kwargs.get('limits_TCI', None)
        modes_TCI = kwargs.get('modes_TCI', {})

        kwargs = {'pk_models': pk_models_TCI, 'pd_models': pd_models_TCI,
                  'data': data_TCI, 'opiates': opiates,
                  'blood_sampling': blood_sampling}

        cp_limit_prop = 10 if limits_TCI is None else limits_TCI.get('cp_limit_prop', 10)
        infusion_limit_prop = 2 if limits_TCI is None else limits_TCI.get('infusion_limit_prop', 2)
        self.tci_prop = TCI.create(Drug.PROPOFOL, modes_TCI.get('prop', TciMode.EFFECT_SITE), cp_limit_prop,
                                   infusion_limit_prop, t_s, **kwargs)

        cp_limit_remi = 10 if limits_TCI is None else limits_TCI.get('cp_limit_remi', 10)
        infusion_limit_remi = 0.5 if limits_TCI is None else limits_TCI.get('infusion_limit_remi', 0.5)
        self.tci_remi = TCI.create(Drug.REMIFENTANIL, modes_TCI.get('remi', TciMode.EFFECT_SITE), cp_limit_remi,
                                   infusion_limit_remi, t_s, **kwargs)

        cp_limit_nore = 15 if limits_TCI is None else limits_TCI.get('cp_limit_nore', 15)
        infusion_limit_nore = 0.5 if limits_TCI is None else limits_TCI.get('infusion_limit_nore', 0.5)
        self.tci_nore = TCI.create(Drug.NOREPINEPHRINE, modes_TCI.get('nore', TciMode.PLASMA), cp_limit_nore,
                                   infusion_limit_nore, t_s, **kwargs)

        cp_limit_rocu = 10 if limits_TCI is None else limits_TCI.get('cp_limit_rocu', 10)
        infusion_limit_rocu = 0.1 if limits_TCI is None else limits_TCI.get('infusion_limit_rocu', 0.1)
        self.tci_rocu = TCI.create(Drug.ROCURONIUM, modes_TCI.get('rocu', TciMode.EFFECT_SITE), cp_limit_rocu,
                                   infusion_limit_rocu, t_s, **kwargs)

        super()._initialize_simulation(patient_data, t_sim, t_s,
                                       opiates, blood_sampling,
                                       internal_states,
                                       output_init,
                                       pk_models,
                                       pd_models,
                                       interaction,
                                       doh_measure,
                                       stimuli,
                                       volume_status,
                                       bring_to_maintenance,
                                       seed)

    def _bring_patient_to_maintenance(self, t_s, doh_measure):
        """
        Helper function to bring the patient to a maintenance state with PID control.
        Note: The PID has been optimized for the 44 patients in the dataset using the PK-PD model of Eleveld et al. for both propofol and remifentanil.
        """
        # Initialize PID parameters
        kp = 1.02
        ki = 0.03
        kd = 0.54

        pid = PID(kp=kp, ki=ki, kd=kd, setpoint=0.5, sample_time=t_s, output_limits=(0, 10))

        measure = 'WAV' if doh_measure in [DoHMeasure.WAV, DoHMeasure.Both] else 'BIS'
        while not self._current_patient.get_patient_phase()[1]:
            doh = self.get_patient_state()[measure]
            u_prop = pid.compute(doh / 100, t_s)
            self.one_step_simulation(u_prop, u_prop / 2, 0, 0)
            if self._current_time > 5000:
                raise ValueError("The patient could not reach a steady state. Please check the PID parameters.")

        # Reset for the maintenance phase
        self._current_patient.clear_data()
        self._current_patient.set_time(0)
        self._current_time = 0

    def init_simulation_from_file(self, id_patient, t_sim, t_s=5,
                                  opiates=None, blood_sampling=None,
                                  internal_states=None,
                                  output_init=None,
                                  pk_models=None,
                                  pd_models=None,
                                  interaction=None,
                                  doh_measure=None,
                                  stimuli=None,
                                  volume_status=None,
                                  bring_to_maintenance=False,
                                  seed=None,
                                  **kwargs):
        """
        Initialize the simulation using patient data from the file in the 'parameters' folder.
        :param id_patient: The index of the patient to simulate.
        :type id_patient: int
        :param t_sim: The simulation time in seconds.
        :type t_sim: int
        :param t_s: The sampling time in seconds. It should be a divisor of t_sim and greater than 1 second. Default is 5.
        :type t_s: int
        :param opiates: True if opiates are used, False otherwise.
        :type opiates: bool
        :param blood_sampling: The type of blood sampling. Default is 'arterial'.
        :type blood_sampling: BloodSampling
        :param internal_states: The internal states of the patient. Default is None.
        :type internal_states: dict
        :param output_init: The initial values for map, hr, sv. Default is None.
        :type output_init: dict
        :param pk_models: The pharmacokinetic models for propofol and remifentanil. Default is None.
        :type pk_models: dict
        :param pd_models: The pharmacodynamic models for propofol and remifentanil. Default is None.
        :type pd_models: dict
        :param interaction: The interaction model. Default is Surface Model.
        :type interaction: Interaction
        :param doh_measure: The depth of hypnosis measure. Default is Both.
        :type doh_measure: DoHMeasure
        :param stimuli: The list of stimuli for the simulation. Default is None.
        :type stimuli: dict
        :param volume_status:The list of volume status for the simulation. Default is None.
        :type volume_status: dict
        :param bring_to_maintenance: True if the patient should be brought to a maintenance state with PID control
        :type bring_to_maintenance: bool
        :param seed: seed for the variability of the patient hemodynamics
        :type seed: int
        :param kwargs: Additional keyword arguments to pass to the TCI initialization
        :type kwargs: dict
        """

        pk_models_TCI = kwargs.get('pk_models_TCI', {})
        pd_models_TCI = kwargs.get('pd_models_TCI', {})

        file_name_patient = 'pd_data_schnider.csv' if pk_models.get('prop',
                                                                    Model.ELEVELD) == Model.SCHNIDER else 'pd_data_eleveld.csv'
        file_name_TCI = 'pd_data_schnider.csv' if pk_models_TCI.get('prop',
                                                                    Model.ELEVELD) == Model.SCHNIDER else 'pd_data_eleveld.csv'

        current_path = os.path.dirname(os.path.abspath(__file__))
        path_patient = os.path.join(current_path, 'parameters')
        patient_data = read_data(file_name_patient, path_patient, row_index=id_patient)
        data_TCI = read_data(file_name_TCI, path_patient, row_index=id_patient)

        limits_TCI = kwargs.get('limits_TCI', None)
        modes_TCI = kwargs.get('modes_TCI', {})

        kwargs = {'pk_models_TCI': pk_models_TCI, 'pd_models_TCI': pd_models_TCI,
                  'data_TCI': data_TCI, 'limits_TCI': limits_TCI, 'modes_TCI': modes_TCI}

        self._initialize_simulation(patient_data, t_sim, t_s,
                                    opiates, blood_sampling,
                                    internal_states,
                                    output_init,
                                    pk_models,
                                    pd_models,
                                    interaction,
                                    doh_measure,
                                    stimuli,
                                    volume_status,
                                    bring_to_maintenance,
                                    seed,
                                    **kwargs)

    def init_simulation_from_data(self, data, t_sim, t_s=5,
                                  opiates=None, blood_sampling=None,
                                  internal_states=None,
                                  output_init=None,
                                  pk_models=None,
                                  pd_models=None,
                                  interaction=None,
                                  doh_measure=None,
                                  stimuli=None,
                                  volume_status=None,
                                  bring_to_maintenance=False,
                                  seed=None,
                                  **kwargs):
        """
        Initialize the simulation with patient parameters directly.
        :param data: The patient data [height, weight, age, gender, E0, k_d , delay, C50P, gammaP, rms_nonlin]
        Note: The parameters [E0, k_d , delay, C50P, gammaP, rms_nonlin] need to be provided only the for PATIENT_SPECIFIC PD model for propofol.
        :type data: list
        :param t_sim: The simulation time in seconds
        :type t_sim: int
        :param t_s: The sampling time in seconds. It should be a divisor of t_sim and greater than 1 second
        :type t_s: int
        :param opiates: True if opiates are used, False otherwise
        :type opiates: bool
        :param blood_sampling: The type of blood sampling to use for the simulation (arterial or venous)
        :type blood_sampling: BloodSampling
        :param internal_states : The internal states of the patient
        :type internal_states: dict
        :param output_init: The initial values for map, hr, sv
        :type output_init: dict
        :param pk_models: The pharmacokinetic models for propofol and remifentanil
        :type pk_models: dict
        :param pd_models: The pharmacodynamic models for propofol and remifentanil
        :type pd_models: dict
        :param interaction: The interaction model to use for the simulation (surface or no_interaction)
        :type interaction: Interaction
        :param doh_measure: The depth of hypnosis measure to use for the simulation (BIS, WAV or both)
        :type doh_measure: DoHMeasure
        :param stimuli: The list of stimuli for the simulation
        :type stimuli: dict
        :param volume_status: The list of volume status for the simulation
        :type volume_status: dict
        :param bring_to_maintenance: True if the patient should be brought to a maintenance state with PID control
        :type bring_to_maintenance: bool
        :param seed: seed for the variability of the patient hemodynamics
        :type seed: int
        :param kwargs: Additional keyword arguments to pass to the TCI initialization
        :type kwargs: dict
        """

        if data is None:
            raise ValueError("No patient data found. Please add the data for a simulation first.")

        data_TCI = kwargs.get('data_TCI', None)
        if data_TCI is None:
            raise ValueError("No TCI data found. Please add the data for a simulation first.")

        self._initialize_simulation(data, t_sim, t_s,
                                    opiates, blood_sampling,
                                    internal_states,
                                    output_init,
                                    pk_models,
                                    pd_models,
                                    interaction,
                                    doh_measure,
                                    stimuli,
                                    volume_status,
                                    bring_to_maintenance,
                                    seed,
                                    **kwargs)

    def run_complete_simulation(self, t_prop: list, t_remi: list, t_nore: list, t_rocu: list):
        """
        Simulates the patient for the given inputs from the start to the end of the simulation.
        :param t_prop: A list of propofol concentrations targets for the TCI in [µg/ml].
        :type t_prop: list
        :param t_remi: A list of remifentanil concentrations targets for the TCI in [ng/ml].
        :type t_remi: list
        :param t_nore: A list of norepinephrine concentrations targets for the TCI in [ng/ml].
        :type t_nore: list
        :param t_rocu: A list of rocuronium concentrations targets for the TCI in [ng/ml].
        :type t_rocu: list
        """

        self._check_run_inputs(t_prop, t_remi, t_nore, t_rocu)

        for time in range(self._t_sim // self._t_s):
            self.one_step_simulation(t_prop[time * self._t_s], t_remi[time * self._t_s], t_nore[time * self._t_s],
                                     t_rocu[time * self._t_s])

    def one_step_simulation(self, t_prop: float, t_remi: float, t_nore: float, t_rocu: float):
        """
        Simulates the patient for a single step with the given inputs.
        :param t_prop: The propofol target concentration for the TCI in [µg/ml].
        :type t_prop: float
        :param t_remi: The remifentanil target concentration for the TCI in [ng/ml].
        :type t_remi: float
        :param t_nore: The norepinephrine target concentration for the TCI in [ng/ml].
        :type t_nore: float
        :param t_rocu: The rocuronium target concentration for the TCI in [ng/ml].
        :type t_rocu: float
        """

        if self._current_patient is None:
            raise ValueError("No patient data found. Please add a simulation first.")

        t_nore = t_nore * 1000 / self.nore_molweight  # Norepinephrine in [nmol/L]

        u_prop = self.tci_prop.compute_infusion(1, t_prop)[0]
        u_remi = self.tci_remi.compute_infusion(1, t_remi)[0]
        u_nore = self.tci_nore.compute_infusion(1, t_nore)[0]
        u_rocu = self.tci_rocu.compute_infusion(1, t_rocu)[0]

        self._current_patient.step(u_prop, u_remi, u_nore, u_rocu, self._t_s)

        self._current_time += self._t_s
