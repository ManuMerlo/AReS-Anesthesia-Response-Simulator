# Third party library imports
import scipy.signal as signal
import numpy as np
import control as ctrl

# Local imports
from .utils.enums import Drug, Model, BloodSampling

from math import sqrt, exp

class PharmacodynamicDoH:
    def __init__(self, model_prop: Model = Model.ELEVELD, model_remi: Model = Model.ELEVELD, **kwargs):
        self.model_prop = model_prop
        self.model_remi = model_remi
        self._e0 = 100                    # The base value of the depth of hypnosis
        self._delay = 0                   # Response delay introduced in the PD model of the propofol  [seconds]

        self.pd_prop_ce = None            # The linear part of propofol PD model used to calculate the effect-site concentration
        self._gam_high_ce = None          # Steepness of the effect-concentration curve for high values of effect-site concentration
        self._gam = None                  # Steepness of the effect-concentration curve for low values of effect-site concentration
        self._ke0_prop = None             # Flow rates from effect-site compartment to the central one for porpofol [1/min]
        self._ec50_prop = None            # Propofol effect-site concentration at 50 % effect [µg/mL]

        self.pd_remi_ce = None            # The linear part of remifentanil PD model used to calculate the effect-site concentration
        self._ec50_remi = None            # Remifentanil effect-site concentration at 50 % effect [ng/mL]
        self._ke0_remi = None             # Flow rates from effect-site compartment to the central one for remifentanil [1/min]

        self.create_propofol_model(**kwargs)
        self.create_remifentanil_model(**kwargs)

    def set_e0(self, e0):
        self._e0 = e0

    def get_e0(self):
        return self._e0

    @staticmethod
    def neurow_sensor_dynamics():
        """
        Create the NeuroSENSE monitor dynamics.
        :return: NeuroSENSE monitor dynamics state space model.
        :rtype: signal.StateSpace
        """
        af = np.array([[0.875000000000001, -0.882496902584599],
                       [1.133148453066824, -1.124999999999994]])
        bf = np.array([[0.067872406960552],
                       [-0.073771149672801]])
        cf = np.array([[0.154464292547276, -0.044802755452488]])
        df = np.array([[0.007190984592330]])
        gns_sim = signal.StateSpace(af, bf, cf, df)
        return gns_sim

    @staticmethod
    def bis_sensor_dynamics():
        """
        Create the Bispectral sensor dynamics.
        :return: delay model and linear time invariant model for the Bispectral sensor dynamics.
        :rtype: tuple of tuple
        """
        bis_delay = 5  # 5 seconds delay
        num_pade, den_pade = ctrl.pade(bis_delay, 1)
        dbis = signal.tf2ss(num_pade, den_pade)
        num = [-0.049, 0.066, 0.976, 0.056]
        den = [1, 2.27, 6.353, 1.487, 0.056]
        bis_lti = signal.tf2ss(num, den)
        return dbis, bis_lti

    @staticmethod
    def get_state_space(ke0):
        """
        Get the state space model for the pharmacodynamic model.
        :param ke0: elimination rate constant.
        :type ke0: float
        :return: state space model for the pharmacodynamic model.
        :rtype: signal.StateSpace
        """
        A = -ke0 / 60
        B = ke0 / 60
        C = 1
        D = 0
        pd_ce = signal.StateSpace(A, B, C, D)
        return pd_ce

    def create_propofol_model(self, **kwargs):
        """
        Create the propofol pharmacodynamic model.
        :param kwargs: additional parameters.
        :type kwargs: dict
        """
        age = kwargs.get('age', -1)
        if age < 0:
            raise ValueError("Age is missing for the PK model for propofol.")

        if self.model_prop == Model.PATIENT_SPECIFIC:
            # Propofol patient-specific PD models (Hosseinirad et al. 2023)
            if 'pd_patient_specific_parameter' not in kwargs:
                raise ValueError("PD patient-specific parameters are missing for the PK model for propofol.")
            data = kwargs.get('pd_patient_specific_parameter', None)
            self._e0 = data['e0']
            self._ke0_prop = data['ke0']
            self._delay = data['delay']
            self._ec50_prop = data['ec50']
            self._gam = data['gamma']
            self._gam_high_ce = self._gam

        elif self.model_prop == Model.SCHNIDER:
            self._ke0_prop = 0.456               # [min ^ (-1)]
            self._ec50_prop = 2.9 - 0.022 * age  # [µg ml ^ (-1)]
            self._delay = 0
            self._gam = 1.43
            self._gam_high_ce = self._gam

        elif self.model_prop == Model.ELEVELD:
            weight = kwargs.get('weight', -1)
            blood_sampling = kwargs.get('blood_sampling', None)

            if weight < 0 or blood_sampling is None:
                raise ValueError("Weight or blood sampling site is missing for the PK model for propofol.")

            thetaPD_1 = 3.08                # ec50 [µg mL^-1]
            thetaPD_2 = 0.146               # ke0 for arterial samples [min^-1]
            thetaPD_3 = 93                  # Baseline BIS value
            thetaPD_4 = 1.47                # PD sigmoid slope (Ce > Ce50)
            thetaPD_7 = -0.00635            # Decrease in Ce50 with age
            thetaPD_8 = 1.24                # ke0 for venous samples [min^-1]
            thetaPD_9 = 1.89                # PD sigmoid slope (Ce < Ce50)

            # The effect of the blood sampling site on the volume and clearance rates
            if blood_sampling == BloodSampling.ARTERIAL:
                theta_ke0 = thetaPD_2
            elif blood_sampling == BloodSampling.VENOUS:
                theta_ke0 = thetaPD_8
            else:
                raise ValueError(f"The {blood_sampling} is not valid.")

            f_aging_theta7 = np.exp(thetaPD_7 * (age - 35))
            self._ec50_prop = thetaPD_1 * f_aging_theta7           # [µg ml ^ (-1)]
            self._ke0_prop = theta_ke0 * (weight / 70) ** (-0.25)  # [min ^ (-1)]
            self._e0 = thetaPD_3
            self._delay = 0
            self._gam = thetaPD_9
            self._gam_high_ce = thetaPD_4
        else:
            raise ValueError(f"Model not supported: {self.model_prop}")

        self.pd_prop_ce = PharmacodynamicDoH.get_state_space(self._ke0_prop)

    def create_remifentanil_model(self, **kwargs):
        """
        Create the remifentanil pharmacodynamic model.
        :param kwargs: additional parameters.
        :type kwargs: dict
        """
        age = kwargs.get('age', -1)
        if age < 0:
            raise ValueError("Age is missing for the PK model for remifentanil.")

        if self.model_remi is None or self.model_remi == Model.MINTO:
            # Remifentanil parameters (Minto et al. 1997)
            self._ke0_remi = 0.595 - 0.007 * (age - 40)
            self._ec50_remi = 13.1 - 0.148 * (age - 40)

        elif self.model_remi == Model.ELEVELD:
            # Remifentanil parameters (Eleveld et al. 2017)
            theta1 = -0.0289
            self._ec50_remi = 12.7  # [ng/ml]
            self._ke0_remi = 1.09 * np.exp(theta1 * (age - 35))  # [1/min]
        else:
            raise ValueError(f"Model not supported: {self.model_remi}")

        self.pd_remi_ce = PharmacodynamicDoH.get_state_space(self._ke0_remi)

    def hillfun(self, x: np.ndarray) -> np.ndarray:
        """
        Hill function for the pharmacodynamic model.
        :param x: input concentration.
        :type x: np.ndarray
        :return: effect of the drug.
        :rtype: np.ndarray
        """
        x = np.array([i if i > 0 else 0 for i in x])
        e = np.zeros_like(x)
        for i in range(len(x)):
            gamma = self._gam if x[i] < self._ec50_prop else self._gam_high_ce
            e[i] = x[i] ** gamma / (self._ec50_prop ** gamma + x[i] ** gamma)
        e = self._e0 - self._e0 * e
        return e

    def responseSurfaceModel(self, x_prop, x_remi):
        """
        Response surface model for the interaction between propofol and remifentanil introduced by Bouillon et al. (2004).
        :param x_prop: propofol concentration.
        :param x_remi: remifentanil concentration.
        :return: effect of the drug.
        :rtype: np.ndarray
        """
        # Beta value is not the one reported by Bouillon et al.(2004), and
        # this value is empirically found by trying different combinations
        # of pk-pd model for propofol and remifentanil such that the
        # difference between the computed doH with interaction and
        # without the interaction is not more than 20% of the value of
        # DoH if we disregard the interactions
        beta = 1.0
        x_prop = np.array(x_prop)
        x_remi = np.array(x_remi)
        inter_prop = x_prop / self._ec50_prop
        inter_remi = x_remi / self._ec50_remi
        theta = inter_prop / (inter_prop + inter_remi + np.finfo(float).eps)  # eps to avoid division by zero
        inter = (inter_prop + inter_remi) / (1 - beta * theta + beta * theta ** 2)
        inter[inter < 0] = 0
        inter = np.array(inter)
        e = self._e0 - self._e0 * (inter ** self._gam / (1 + inter ** self._gam))
        return e

    def get_delay_ss(self):
        """
        Get the state space model for the delay.
        :return: state space model for the delay.
        :rtype: tuple
        """
        num, den = ctrl.pade(self._delay, 1)
        pd_ce_delayed = signal.tf2ss(num, den)
        return pd_ce_delayed


class PharmacodynamicHemo:
    """
    TPR: total peripheral resistance [mmHg mL^-1 min]
    SV: stroke volume [mL]
    HR: heart rate [beats min^-1]
    MAP: mean arterial pressure [mmHg]
    CO: cardiac output [L/min]
    TDE: time-dependent effect

    Please refer to the following two papers for this PD model:
    1) Pharmacodynamic mechanism-based interaction model for the haemodynamic effects of remifentanil and propofol
       in healthy volunteers. Su et al., Br J Anaesthesia, 2023.
    2) Design of a pharmacokinetic/pharmacodynamic model for administration of low dose peripheral norepinephrine during
       general anaesthesia. Joachim et a., Br J Clin Pharmacol, 2024

    """

    def __init__(self, age):
        super().__init__()
        # Private properties
        self._emax_tpr_prop = -0.778                                       # maximum effect of propofol on TPR
        self._ec50_tpr_prop = 3.21                                         # c_e that produce half of the maximal propofol effect of TPR [µg ml^-1]
        self._emax_sv_age = 0.0333
        self._emax_sv_prop = -0.154 * exp(self._emax_sv_age * (age - 35))  # maximum effect of propofol on SV
        self._ec50_sv_prop = 0.44                                          # c_e that produce half of the maximal propofol effect of SV [µg ml^-1]
        self._emax_tpr_remi = -1                                           # maximum effect of remifentanil on TPR
        self._ec50_tpr_remi = 4.59                                         # c_e that produce half of the maximal remifentanil effect of TPR [ng ml^-1]
        self._sl_sv_remi = 0.0581                                          # slope of remifentanil effect on SV [ng ml^-1]
        self._sl_hr_remi = 0.0327                                          # slope of remifentanil effect on HR [ng ml^-1]

        # maximum magnitude changes of slope of remifentanil on tpr SV and HR caused by propofol
        self._int_tpr = 1
        self._int_hr = -0.119           # [ng ml^-1]
        self._int_sv = -0.212           # [ng ml^-1]
        self._ec50_int_hr = 0.196       # Interaction potency of propofol on the effect of remifentanil on HR [microg ml^-1]

        # Hill coefficient for TPR sigmoid
        self._gamma_tpr_prop = 1.83
        self._gamma_tpr_remi = 1

        self._fb = 0.661  # magnitude of the feedback
        self._kout = 0.072 / 60         # first-order dissipation rate constant for SV, HR, and TPR [second^-1]
        self._k_tde = 0.067 / 60        # first-order dissipation rate constant for TDE [second^-1]

        # Inter - patient variability
        self._omega_base_sv = sqrt(0.0328)
        self._omega_base_tpr = sqrt(0.0528)
        self._omega_base_hr = sqrt(0.0242)
        self._omega_ec50_tpr_prop = sqrt(0.44)
        self._omega_emax_tpr_remi = sqrt(0.449)
        self._omega_sl_hr_remi = sqrt(0.00382)
        self._omega_sl_sv_remi = sqrt(0.00868)

        # Hill function parameters for norepinephrine effect on hemodynamics variables
        self._ec50_nore_map = 15.27         # [nmol l^-1]
        self._ec50_nore_co = 36             # [nmol l^-1]
        self._gamma_nore_map = 1.46
        self._gamma_nore_co = 2.3
        self._t_lag = 25                    # [second]

        # Public properties
        self.hr_sv = 0.312                  # magnitude of the inverse effect of HR on SV
        self.ltde_hr = 0.121                # percentage of increased baseline HR caused by the time-dependent effect
        self.ltde_sv = 0.0899               # percentage of increased baseline HR caused by the time-dependent effect
        self.base_sv = 82.2                 # base value of SV [ml]
        self.base_hr = 56.1                 # base value of HR [beats min^-1]
        self.base_tpr = 0.0163              # base value of TPR [mmHg ml^-1 min]

    def interpatient_variability(self, seed):
        """
        Create model 's parameters with interpatient_variability
        :param seed: seed for random number generator
        :return:
        """
        if seed is not None:
            np.random.seed(seed)  # seed for random number generator to reproduce the same results
            # if the seed is not provided, the random number generator will use the system time
        var_rand = np.random.rand(1, 1) * 0.5
        var_rand = var_rand[0][0]
        self.base_sv = self.base_sv * exp(self._omega_base_sv * var_rand)
        self.base_hr = self.base_hr * exp(self._omega_base_hr * var_rand)
        self.base_tpr = self.base_tpr * exp(self._omega_base_tpr * var_rand)
        self._ec50_tpr_prop = self._ec50_tpr_prop * exp(self._omega_ec50_tpr_prop * var_rand)
        self._emax_tpr_remi = self._emax_tpr_remi + self._omega_emax_tpr_remi * var_rand
        self._sl_hr_remi = self._sl_hr_remi + self._omega_sl_hr_remi * var_rand
        self._sl_sv_remi = self._sl_sv_remi + self._omega_sl_sv_remi * var_rand

    def hillfun(self, x: np.ndarray, output: str):
        """
        Hill function for the pharmacodynamic model.
        :param x: input concentration.
        :type x: np.ndarray
        :param output: output of the pharmacodynamic model, either MAP or CO.
        :type output: str
        :return: effect of the drug on the output.
        :rtype: np.ndarray
        """
        # Norepinephrine effect on MAP and CO.
        if output == 'map':
            gamma = self._gamma_nore_map
            ec50 = self._ec50_nore_map                                                                      # [mmHg]
            emax = self.base_hr * (1 + self.ltde_hr) * self.base_sv * (1 + self.ltde_sv) * self.base_tpr    # [mmHg]
        elif output == 'co':
            gamma = self._gamma_nore_co
            ec50 = self._ec50_nore_co                                                                       # [L/min]
            emax = self.base_hr * (1 + self.ltde_hr) * self.base_sv * (1 + self.ltde_sv) / 1000             # [L/min]
        else:
            raise ValueError(f"Norepinephrine is not considered to affect: {output}")

        x = np.array([i if i > 0 else 0 for i in x])
        # Hill function
        e = x ** gamma / (ec50 ** gamma + x ** gamma)
        e = emax * e
        return e

    def prop_remi_inter(self, cp_prop: float, cp_remi: float):
        """
        Propofol and remifentanil interaction.
        :param cp_prop: plasma concentration of propofol.
        :type cp_prop: float
        :param cp_remi: plasma concentration of remifentanil.
        :type cp_remi: float
        :return: propofol and remifentanil interaction.
        :rtype: tuple of float
        """
        tpr_prop = ((self._emax_tpr_prop + self._int_tpr * cp_remi /
                     (self._ec50_tpr_remi + cp_remi)) * cp_prop ** self._gamma_tpr_prop /
                    (self._ec50_tpr_prop ** self._gamma_tpr_prop + cp_prop ** self._gamma_tpr_prop))

        sv_prop = (self._emax_sv_prop * cp_prop) / (self._ec50_sv_prop + cp_prop)

        tpr_remi = (self._emax_tpr_remi * cp_remi ** self._gamma_tpr_remi /
                    (self._ec50_tpr_remi ** self._gamma_tpr_remi + cp_remi ** self._gamma_tpr_remi))

        sv_remi = (self._sl_sv_remi + self._int_sv * cp_prop / (self._ec50_sv_prop + cp_prop)) * cp_remi

        hr_remi = (self._sl_hr_remi + self._int_hr * cp_prop / (self._ec50_int_hr + cp_prop)) * cp_remi

        tpr_prop = np.max([tpr_prop, -0.999])
        sv_prop = np.max([sv_prop, -0.999])
        tpr_remi = np.min([tpr_remi, 0.999])
        sv_remi = np.min([sv_remi, 0.999])
        hr_remi = np.min([hr_remi, 0.999])

        return tpr_prop, tpr_remi, sv_prop, sv_remi, hr_remi

    def ode_prop_hemodynamic(self, t, y: np.ndarray, cp_prop: float = 0, cp_remi: float = 0) -> np.ndarray:
        """
        ODE for nonlinear model of the evolution of hemodynamic parameters.
        :param t: time ( not used in the function, but required to match the signature of the function).
        :type t: float
        :param y: state variables.
        :type y: np.ndarray
        :param cp_prop: plasma concentration of propofol.
        :type cp_prop: float
        :param cp_remi: plasma concentration of remifentanil.
        :type cp_remi: float
        :return: derivative of the state variables.
        :rtype: np.ndarray
        """
        dydt = np.zeros(5)

        tpr_prop, tpr_remi, sv_prop, sv_remi, hr_remi = self.prop_remi_inter(cp_prop, cp_remi)

        k_in_hr = self._kout * self.base_hr
        k_in_sv = self._kout * self.base_sv
        k_in_tpr = self._kout * self.base_tpr

        tde_hr = y[3] + self.base_hr * self.ltde_hr
        tde_sv = y[4] + self.base_sv * self.ltde_sv

        hr = y[2] + tde_hr
        sv = (y[1] + tde_sv) * (1 - self.hr_sv * np.log(hr / (self.base_hr * (1 + self.ltde_hr))))

        rmap = y[0] * hr * sv / (
                self.base_hr * (1 + self.ltde_hr) * self.base_sv * (1 + self.ltde_sv) * self.base_tpr)

        dydt[0] = k_in_tpr * rmap ** (-self._fb) * (1 + tpr_prop) - self._kout * y[0] * (1 - tpr_remi)
        dydt[1] = k_in_sv * rmap ** (-self._fb) * (1 + sv_prop) - self._kout * y[1] * (1 - sv_remi)
        dydt[2] = k_in_hr * rmap ** (-self._fb) - self._kout * y[2] * (1 - hr_remi)
        dydt[3] = -self._k_tde * y[3]
        dydt[4] = -self._k_tde * y[4]
        # y[0] = TPR, y[1] = SV, y[2] = HR, y[3] = TDE_HR, y[4] = TDE_SV
        return dydt

    def get_nore_delay_ss(self):
        """
        Get the state space model for the norepinephrine delay corresponding to the depot compartment.
        :return: state space model for the norepinephrine delay.
        :rtype: tuple
        """
        num, den = ctrl.pade(self._t_lag, 1)
        pd_delayed = signal.tf2ss(num, den)
        return pd_delayed


class PharmacodynamicNMB:
    """
    NMB: Nueromuscular Blockade
    m = 0 (NMB is moderate)
    m = 1 (NMB is deep)
    m = 2 (NMB is profound)
    m = 3 (NMB is very profound)

    Please refer to the following paper for this PD model:

    Comparison of two pharmacokinetic–pharmacodynamic models of rocuronium bromide
    during profound neuromuscular block: analysis of estimated and measured
    post-tetanic count effect, M. Couto et al., BJ Anaesthesia (2022)
    """
    def __init__(self, drug: Drug = Drug.ROCURONIUM):
        self._drug = drug

        if self._drug == Drug.ROCURONIUM:
            self._ke0 = 0.134                  # Flow rates from effect-site compartment to the central one for rocuronium [1/min]
            self._ec50_m1 = 1.10               # Rocuronium effect-site concentration associated with 50% of probability for the m=1 category [µg/mL]
            self._ec50_m2 = 1.5                # Rocuronium effect-site concentration associated with 50% of probability for the m=2 category [µg/mL]
            self._ec50_m3 = 2.2                # Rocuronium effect-site concentration associated with 50% of probability for the m=3 category [µg/mL]
            self._gam = 4.5                    # Steepness of the effect-concentration curve
        else:
            raise ValueError(f"Drug not supported: {self._drug}")

        self.pd_ce = self.get_state_space()    # Linear part of PD model to calculate the rocuronium effect-site concentration

    def hillfun(self, x: np.ndarray) -> np.ndarray:
        """
        Hill function for the pharmacodynamic model.
        :param x: input concentration.
        :type x: np.ndarray
        :return: effect of the drug.
        :rtype: np.ndarray
        """

        x = [i if i > 0 else 0 for i in x]
        x = np.array(x)

        p1 = x ** self._gam / (self._ec50_m1 ** self._gam + x ** self._gam)
        p2 = x ** self._gam / (self._ec50_m2 ** self._gam + x ** self._gam)
        p3 = x ** self._gam / (self._ec50_m3 ** self._gam + x ** self._gam)

        p_m0 = 1 - p1    # Probability for the m=0 category
        p_m1 = p1 - p2   # Probability for the m=1 category
        p_m2 = p2 - p3   # Probability for the m=2 category
        p_m3 = p3        # Probability for the m=3 category
        return np.array([p_m0, p_m1, p_m2, p_m3])

    def get_state_space(self):
        """
        Get the state space model for the pharmacodynamic model.
        :return: state space model for the pharmacodynamic model.
        :rtype: signal.StateSpace
        """
        A = -self._ke0 / 60
        B = self._ke0 / 60
        C = 1
        D = 0
        pd_ce = signal.StateSpace(A, B, C, D)
        return pd_ce
