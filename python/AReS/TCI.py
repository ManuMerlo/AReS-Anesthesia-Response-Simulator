from .utils.enums import Drug, Model, TciMode, BloodSampling
from .pharmacokinetics import Pharmacokinetic
from .pharmacodynamics import PharmacodynamicDoH, PharmacodynamicNMB

import numpy as np
import scipy.signal as signal

# Custom class for discrete state-space systems
class StateSpaceDiscrete:
    """
    Class to define discrete state-space systems
    """
    def __init__(self, A, B, C, D, dt):
        self.A = A
        self.B = B
        self.C = C
        self.D = D
        self.dt = dt

class TCI:
    """
        TCI: Target Controlled Infusion

        Based on the targeted compartments, the TCI mode can be either
        effect-site [Shafer et al. 1992] or plasma targeted [Bailey et al. 1991].

        The limits of plasma concentration and infusion rate are determined
        for safety considerations [Van Poucke et al. 2004]
    """
    def __init__(self, drug, mode, cp_limit, infusion_limit, time_step):
        self.drug = drug                       # Type of the drug. The infusion rate each four drug can be computed using the TCI method.
        self.mode = mode                       # Either effect-site or plasma targeted TCI
        self.time_step = time_step             # Discretization time step
        self.cp_limit = cp_limit               # Maximum allowed plasma concentration
        self.infusion_limit = infusion_limit   # Maximum implementable drug's infusion rate

        self.pkpd_model = None                 # A discrete state space model (four states)
        self.pk_model = None                   #  Continuous pk model object
        self.pd_model = None                   # Continuous pd model object

        self.u = []                            # Computed infusion rate
        self.target = None                     # Target value of effect-site or plasma concentration. The unit is either [µg/mL] or [ng/mL] depending on the type of the drug
        self.x_internal = None                 # Internal state of the PK-PD model
        self._zero_input_plasma = None         # Plasma concentration response to zero infusion rate where cp(t=0) = 1
        self._step_plasma = None               # Plasma concentration response to u(t) = 1
        self._zero_input_effect = None         # Effect-site concentration response to zero infusion rate where ce(t=0) = 1
        self._impulse_effect = None            # Effect-site concentration response to u(t=0) = 1

    @classmethod
    def create(self, drug: Drug, mode: TciMode, cp_limit: float, infusion_limit: float, time_step: int, **kwargs):
        """
        Create a TCI object based on the drug and mode of administration.
        :param drug: Drug to be administered (Propofol, Remifentanil, Norepinephrine, Rocuronium)
        :type drug: Drug
        :param mode: Mode of administration (Plasma, EffectSite)
        :type mode: TciMode
        :param cp_limit: target concentration limit for the drug
        :type cp_limit: float
        :param infusion_limit: infusion limit for the drug
        :type infusion_limit: float
        :param time_step: time step for the simulation
        :type time_step: float
        :param kwargs: additional parameters for pharmacokinetic and pharmacodynamic to be passed to the TCI object
        :returns: TCI object
        :rtype: TCI
        """
        tci = TCI(drug, mode, cp_limit, infusion_limit, time_step)

        data = kwargs['data']

        # Patient demographic parameters
        height = data[0]
        weight = data[1]
        age = data[2]
        gender = data[3]  # 0 for male and 1 for female
        bmi = data[4]
        lbm = data[5]

        opiates = kwargs.get('opiates', True)
        blood_sampling = kwargs.get('blood_sampling', BloodSampling.ARTERIAL)

        pk_models = kwargs.get('pk_models', {})
        pd_models = kwargs.get('pd_models', {})

        kwargs = {'age': age, 'weight': weight,
                  'height': height, 'gender': gender, 'bmi': bmi, 'lbm': lbm,
                  'opiates': opiates, 'blood_sampling': blood_sampling}

        if pd_models.get('prop', Model.PATIENT_SPECIFIC) == Model.PATIENT_SPECIFIC:
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

        # Based on the drug, different setting is used to initialize the TCI system.
        if drug == Drug.PROPOFOL:
            tci.tci_setting_propofol(pk_models=pk_models, pd_models=pd_models, **kwargs)
        elif drug == Drug.REMIFENTANIL:
            tci.tci_setting_remifentanil(pk_models=pk_models, pd_models=pd_models, **kwargs)
        elif drug == Drug.NOREPINEPHRINE:
            tci.tci_setting_norepinephrine(pk_models=pk_models, **kwargs)
        elif drug == Drug.ROCURONIUM:
            tci.tci_setting_rocuronium(pk_models=pk_models, **kwargs)
        # The prediction horizons to compute zero input, step, and impulse responses.
        pred_horizon_plasma = 7
        pred_horizon_effect = 50

        # For either of TCI mode, the internal model is defined, and the necessary concentration responses are computed.
        if tci.mode == TciMode.PLASMA:
            discrete_system = signal.cont2discrete((tci.pk_model.A, tci.pk_model.B, tci.pk_model.C, tci.pk_model.D),
                                                   tci.time_step)
            # Unpack the components and create a discrete model
            Ad, Bd, Cd, Dd, dt = discrete_system
            tci.pkpd_model = StateSpaceDiscrete(Ad, Bd, Cd, Dd, dt)
            # Response when no input is applied and the initial state is one
            tci._zero_input_plasma = tci.zeroInputResponse(tci.pkpd_model.C[0,:], pred_horizon_plasma)
            # Response when a step input of 1 is applied during the prediction horizon
            tci._step_plasma = tci.stepResponse(tci.pkpd_model.C[0,:], pred_horizon_plasma)
        elif tci.mode == TciMode.EFFECT_SITE:
            A = np.block([
                [tci.pk_model.A, np.zeros((3, 1))],
                [np.hstack([tci.pd_model.B, np.zeros((tci.pd_model.B.shape[0], 2))]), tci.pd_model.A]
            ])
            B = np.vstack([tci.pk_model.B, np.array([[0]])])
            C = np.eye(4)
            D = 0
            discrete_system = signal.cont2discrete((A, B, C, D), tci.time_step)
            Ad, Bd, Cd, Dd, dt = discrete_system
            tci.pkpd_model = StateSpaceDiscrete(Ad, Bd, Cd, Dd, dt)
            # Use the discrete system matrices for responses
            # Response when no input is applied and the initial state is one for the plasma compartment
            tci._zero_input_plasma = tci.zeroInputResponse(tci.pkpd_model.C[0,:], pred_horizon_plasma)
            # Response when a step input of 1 is applied during the prediction horizon for the plasma compartment
            tci._step_plasma = tci.stepResponse(tci.pkpd_model.C[0,:], pred_horizon_plasma)
            # Response when no input is applied and the initial state is one for the effect compartment
            tci._zero_input_effect = tci.zeroInputResponse(tci.pkpd_model.C[3,:], pred_horizon_effect)
            # Response when an impulse input is applied during the prediction horizon for the effect compartment
            tci._impulse_effect = tci.impulseResponse_ce(tci.pkpd_model.C[3,:], pred_horizon_effect)

        # Initialize the states of the internal model
        tci.x_internal = np.zeros((tci.pkpd_model.A.shape[1], 1))
        return tci

    def tci_setting_propofol(self, pk_models, pd_models, **kwargs):
        """
        Set the pharmacokinetic and pharmacodynamic models for propofol.
        :param pk_models: pharmacokinetic model for propofol
        :type pk_models: dict
        :param pd_models: pharmacodynamic model for propofol
        :type pd_models: dict
        :param kwargs: additional parameters for pharmacokinetic and pharmacodynamic models
        """
        pk_model = pk_models.get('prop', Model.ELEVELD)
        pd_model = pd_models.get('prop', Model.ELEVELD)

        # Pharmacokinetic models
        self.pk_model = Pharmacokinetic.create(Drug.PROPOFOL, pk_model, **kwargs).ss
        self.pd_model = PharmacodynamicDoH(model_prop=pd_model, **kwargs).pd_prop_ce

    def tci_setting_remifentanil(self, pk_models, pd_models, **kwargs):
        """
        Set the pharmacokinetic and pharmacodynamic models for remifentanil.
        :param pk_models: pharmacokinetic model for remifentanil
        :type pk_models: dict
        :param pd_models: pharmacodynamic model for remifentanil
        :type pd_models: dict
        :param kwargs: additional parameters for pharmacokinetic and pharmacodynamic models
        :type kwargs: dict
        """
        pk_model = pk_models.get('remi', Model.ELEVELD)
        pd_model = pd_models.get('remi', Model.ELEVELD)

        self.pk_model = Pharmacokinetic.create(Drug.REMIFENTANIL, pk_model, **kwargs).ss
        self.pd_model = PharmacodynamicDoH(model_remi=pd_model, **kwargs).pd_remi_ce

    def tci_setting_rocuronium(self, pk_models, **kwargs):
        """
        Set the pharmacokinetic and pharmacodynamic models for rocuronium.
        :param pk_models: only one type of pk-pd model exists for rocuronium [Model]
        :type pk_models: dict
        :param kwargs: additional parameters for pharmacokinetic and pharmacodynamic models
        :type kwargs: dict
        """
        pk_model = pk_models.get('rocu', Model.DAHE)
        self.pk_model = Pharmacokinetic.create(Drug.ROCURONIUM, pk_model, **kwargs).ss
        self.pd_model = PharmacodynamicNMB().pd_ce

    def tci_setting_norepinephrine(self, pk_models, **kwargs):
        """
        Set the pharmacokinetic model for norepinephrine.
        :param pk_models: only one type of pk model exists for norepinephrine [Model]
        :type pk_models: dict
        :param kwargs: additional parameters for pharmacokinetic model
        :type kwargs: dict
        """
        pk_model = pk_models.get('nore', Model.JOACHIM)
        self.pk_model = Pharmacokinetic.create(Drug.NOREPINEPHRINE, pk_model, **kwargs).ss
        # Since no effect-site concentration is defined for norepinephrine, the only available mode for this drug is plasma targeted.
        self.mode = TciMode.PLASMA

    def reset_state(self):
        """
        Reset the internal state of the TCI object.
        """
        self.x_internal = np.zeros(self.pkpd_model.A.shape[1])

    def compute_infusion(self, n_step, target):
        """
        Compute the infusion rate for the drug based on the target concentration.
        :param n_step: number of steps for the simulation
        :type n_step: int
        :param target: target concentration for the drug
        :type target: float
        :return: list of infusion rates
        :rtype: list
        """
        self.target = target
        self.u = []
        for i in range(n_step):
            u_interv = self.tci_interv(self._zero_input_plasma, self._step_plasma, self._zero_input_effect,
                                       self._impulse_effect, self.x_internal)
            # Apply infusion rate limit
            if u_interv > self.infusion_limit:
                u_interv = self.infusion_limit

            # Update the internal states
            self.x_internal = self.pkpd_model.A @ self.x_internal + self.pkpd_model.B * u_interv

            self.u.append(u_interv)

        return self.u

    def tci_interv(self, _zero_input_plasma, step_input_plasma, _zero_input_effect, _impulse_effect, x0):
        """
        Compute the infusion rate for the drug based on the target concentration.
        :param _zero_input_plasma: zero input response for the plasma compartment
        :type: numpy.ndarray
        :param step_input_plasma: step response for the plasma compartment
        :type: numpy.ndarray
        :param _zero_input_effect: zero input response for the effect compartment
        :type: numpy.ndarray
        :param _impulse_effect: impulse response for the effect compartment
        :type: numpy.ndarray
        :param x0: initial state of the system
        :type: numpy.ndarray
        :return: infusion rate
        :type: float
        """

        # if the concentration of the effect site is within 5% of the target, use the plasma control algorithm
        if self.mode == TciMode.PLASMA or (
                self.mode == TciMode.EFFECT_SITE and abs(x0[3] - self.target) < self.target * 0.05):
            yfree = _zero_input_plasma @ x0
            target_arr = self.target * np.ones((len(yfree), 1))
            # The infusion rate is computed based on the difference between the target and the concentration of the effect site
            # We want: target = yfree + step_response * u_interv
            # Since the equation may not have an exact solution we use the least-squares method to find the best u that minimizes the error.
            u_interv = 1 / (step_input_plasma.T @ step_input_plasma) * step_input_plasma.T @ (target_arr - yfree)
            while isinstance(u_interv, list):
                u_interv = u_interv[0]
        elif self.mode == TciMode.EFFECT_SITE:
            yfree = _zero_input_effect @ x0
            peakfree = np.argmax(yfree)
            peak = yfree[peakfree]
            # if the peak concentration is greater than the target, then the infusion rate is 0
            if peak > self.target:
                u_interv = 0
            else:
                # Simultaneously solve for tpeak and u as proposed in [Shafer et al. 1992]
                tpeak = peakfree   # tpeak is the time at which the peak concentration is reached
                tpeakold = -1      # To take into account the mismatch in the indexing of matlab and python
                infusion = 0
                while tpeak != tpeakold:
                    tpeakold = tpeak
                    infusion = (self.target - yfree[tpeak]) / _impulse_effect[tpeak]
                    if isinstance(infusion, list):
                        infusion = infusion[0]
                    tpeak = np.argmax(yfree + _impulse_effect * infusion)
                u_interv = infusion
                x_new = self.pkpd_model.A @ x0 + self.pkpd_model.B * u_interv
                c_p0 = x0[0]
                c_p_max = x_new[0]
                if c_p_max > self.cp_limit:
                    u_interv = u_interv * ((self.cp_limit - c_p0) / (c_p_max - c_p0))
        else:
            raise ValueError(f'The target concentration {self.mode} is not supported')
        u_interv = u_interv if u_interv > 0 else 0

        return u_interv

    def impulseResponse_ce(self, c, predictionhorizon):
        """
        Compute the impulse response of the system.
        :param c: effect-site concentration part of C matrix from state-space model
        :type c: numpy.ndarray
        :param predictionhorizon: prediction horizon for the simulation
        :type predictionhorizon: int
        :return: impulse response
        :rtype: numpy.ndarray
        """
        # The input is an impulse function u[k] = 1 for k = 0 and u[k] = 0 for k > 0
        # x = A @ x + B @ 1 = A @ x + B
        # y = C @ x + D @ 1 = C @ x + D but D = 0

        # x[0] = A @ x[0] + B @ 1 = B
        # y[0] = C @ x[0] = C @ B
        # y[k] = C @ A^k @ B

        a = self.pkpd_model.A
        b = self.pkpd_model.B
        vDtemp = b
        impulse_resp = np.zeros((predictionhorizon, 1))
        impulse_resp[0] = c @ b
        for i in range(1, predictionhorizon):
            vDtemp = a @ vDtemp
            impulse_resp[i,0] = c @ vDtemp
        return impulse_resp[:-1]  # starting at the response at the next step

    def stepResponse(self, c, predictionhorizon):
        """
        Compute the step response of the system.
        :param c: plasma concentration part of C matrix from state-space model
        :type c: numpy.ndarray
        :param predictionhorizon: prediction horizon for the simulation
        :type predictionhorizon: int
        :return: step response
        :rtype: numpy.ndarray
        """
        # The input is a step function u[k] = 1 for all k >= 0
        # x = A @ x + B @ 1 = A @ x + B
        # y = C @ x + D @ 1 = C @ x + D but D = 0

        a = self.pkpd_model.A
        b = self.pkpd_model.B
        stepResp = np.zeros((predictionhorizon, 1))
        x = np.zeros((a.shape[1], 1))
        for i in range(predictionhorizon):
            x = a @ x + b @ np.array([[1]])
            stepResp[i,0] = c @ x
        return stepResp[:-1]  # starting at the response at the next step

    def zeroInputResponse(self, c, predictionhorizon):
        """
        Compute the zero input response of the system.
        :param c: effect-site or plasma concentration part of C matrix from state-space model
        :type c: numpy.ndarray
        :param predictionhorizon: prediction horizon for the simulation
        :type predictionhorizon: int
        :return: zero input response
        :rtype: numpy.ndarray
        """
        # The input is zero
        # x = A @ x + B @ 0 = A @ x
        # y = C @ x + D @ 0 = C @ x

        # suppose x[0] = 1, then
        # x[1] = A @ x[0] = A
        # x[2] = A @ x[1] = A @ A
        # x[k] = A^k
        # y = C @ A^k

        a = self.pkpd_model.A
        zero_resp = np.zeros((predictionhorizon, a.shape[1]))
        zero_resp[0] = c
        for i in range(1, predictionhorizon):
            zero_resp[i] = zero_resp[i - 1,:] @ a
        return zero_resp[1:,:]  # starting with the response at the next step.

