import math
import numpy as np
import scipy.signal as signal

from .utils.enums import Drug, Model, BloodSampling


class Pharmacokinetic:
    def __init__(self, drug: Drug, model: Model):
        self.drug = drug
        self.model = model
        self.ss = None

        # Three compartment model parameters
        self._v1comp3 = None    # Volume of compartment 1 [L]
        self._v2comp3 = None    # Volume of compartment 2 [L]
        self._v3comp3 = None    # Volume of compartment 3 [L]
        self._cl1comp3 = None   # Clearance rate of compartment 1 [L/min]
        self._cl2comp3 = None   # Clearance rate of compartment 2 [L/min]
        self._cl3comp3 = None   # Clearance rate of compartment 3 [L/min]

        self._k10 = None        # Drug flow rate from compartment 1 to outside of body [1/min]
        self._k12 = None        # Drug flow rate from compartment 1 to 2 [1/min]
        self._k13 = None        # Drug flow rate from compartment 1 to 3 [1/min]
        self._k21 = None        # Drug flow rate from compartment 2 to 1 [1/min]
        self._k31 = None        # Drug flow rate from compartment 2 to 3 [1/min]
        self._ka = None         # Drug flow rate from depot compartment to central one [1/min]

    # Factory method to create the pharmacokinetic model for each drug
    @classmethod
    def create(cls, drug: Drug, model: Model, **kwargs):
        """
        Create a pharmacokinetic model for a given drug and model.
        :param drug: Drug enum (Propofol, Remifentanil, Norepinephrine, Rocuronium)
        :type drug: Drug
        :param model: Model enum
        :type model: Model
        :param kwargs: additional parameters for the model
        :type kwargs: dict
        :return: Pharmacokinetic model
        """
        pk_model = cls(drug=drug, model=model)
        if drug == Drug.PROPOFOL:
            pk_model.create_propofol(model, **kwargs)
        elif drug == Drug.REMIFENTANIL:
            pk_model.create_remifentanil(model, **kwargs)
        elif drug == Drug.NOREPINEPHRINE:
            pk_model.create_norepinephrine(model, **kwargs)
        elif drug == Drug.ROCURONIUM:
            pk_model.create_rocuronium(model, **kwargs)
        else:
            raise ValueError(f"Drug {drug} not supported")

        return pk_model

    # The helper functions defined by Eleveld et al.(2018) for propfol
    # pk model to account for demographic parameters
    @staticmethod
    def f_aging(x, age):
        return np.exp(x * (age - 35))

    @staticmethod
    def f_sigmoid(x, e50, gamma):
        return x ** gamma / (x ** gamma + e50 ** gamma)

    def create_propofol(self, model: Model, **kwargs):
        """
        Create a pharmacokinetic model for propofol.
        :param model: Model enum (Schnider, Eleveld)
        :type model: Model
        :param kwargs: additional parameters for the model
        :type kwargs: dict
        """
        age = kwargs.get('age', -1)
        weight = kwargs.get('weight', -1)

        if any(x == -1 for x in [age, weight]):
            raise ValueError("Age or weight is missing for Propofol PK model")

        if model == Model.SCHNIDER:
            # Propofol parameters (Schnider et al. 1998)
            height = kwargs.get('height', -1)
            lbm = kwargs.get('lbm', -1)

            if any(x is None for x in [height, lbm]):
                raise ValueError("Height or lean body mass is missing for Schnider PK model for Propofol")

            self._v1comp3 = 4.27  # [l]
            self._v2comp3 = 18.9 - 0.391 * (age - 53)  # [l]
            self._v3comp3 = 238  # [l]
            self._cl1comp3 = 1.89 + 0.0456 * (weight - 77) - 0.0681 * (lbm - 59) + 0.0264 * (height - 177)
            self._cl2comp3 = 1.29 - 0.024 * (age - 53)
            self._cl3comp3 = 0.836  # [l / min]

        elif model == Model.ELEVELD:
            # Propofol parameters (Eleveld et al. 2018)
            gender = kwargs.get('gender', -1)
            bmi = kwargs.get('bmi', -1)
            opiates = kwargs.get('opiates', None)
            blood_sampling = kwargs.get('blood_sampling', None)

            if any(x == -1 for x in [gender, bmi]) or opiates is None or blood_sampling is None:
                raise ValueError("Gender, BMI, opiates or blood sampling is missing for Eleveld PK model for Propofol")

            theta_1 = 6.28  # V1_ref[litre]
            theta_2 = 25.5  # V2_ref[litre]
            theta_3 = 273  # V3_ref[litre]
            theta_4 = 1.79  # CL_ref(male)[Litre / min]
            theta_5 = 1.75  # Q2_ref[Litre / min]
            theta_6 = 1.11  # Q3_ref[Litre / min]
            # theta_7 = 0.191  # Typical residual error
            theta_8 = 42.3  # CL maturation E50[weeks]
            theta_9 = 9.06  # CL maturation E50 slope
            theta_10 = -0.0156  # Smaller V2 with age
            theta_11 = -0.00286  # Lower CL with age
            theta_12 = 33.6  # Weight for 50% of maximal V1[kg]
            theta_13 = -0.0138  # Smaller V3 with age
            theta_14 = 68.3  # Maturation of Q3[weeks]
            theta_15 = 2.10  # CL_ref(female)[Litre / min]
            theta_16 = 1.3  # Higher Q2 for maturation of Q3
            theta_17 = 1.42  # V1 venous samples(children)
            theta_18 = 0.68  # Higher Q2 venous samples

            # Intermediate Coefficients
            pma = age * 365.25 / 7 + 40  # Post Menstural Age in weeks
            f_agingTheta10 = Pharmacokinetic.f_aging(theta_10, age)
            f_centralWGT = weight / (weight + theta_12)
            f_centralWGTRef = 70 / (70 + theta_12)
            f_CLmaturation = pma ** theta_9 / (pma ** theta_9 + theta_8 ** theta_9)
            f_CLmaturation_ref = 1866 ** theta_9 / (1866 ** theta_9 + theta_8 ** theta_9)
            f_Q3maturation = (age * 365.25 / 7 + 40) / (age * 365.25 / 7 + 40 + theta_14)
            f_Q3maturation_ref = (35 * 365.25 / 7 + 40) / (35 * 365.25 / 7 + 40 + theta_14)

            # The impact of gender on the intermediate coefficients
            if gender:  # female = 1
                f_AlSallami = (1.11 + (1 - 1.11) / (1 + (age / 7.1) ** (-1.1))) * (
                        9270 * weight / (8780 + 244 * bmi))
                theta_CL = theta_15
            else:  # male = 0
                f_AlSallami = (0.88 + (1 - 0.88) / (1 + (age / 13.4) ** (-12.7))) * (
                        9270 * weight / (6680 + 216 * bmi))
                theta_CL = theta_4

            f_AlSallami_ref = (0.88 + (1 - 0.88) / (1 + (35 / 13.4) ** (-12.7))) * (9270 * 70 / (6680 + 216 * 24.2))

            # The impact of opiates' presence on the intermediate coefficients
            if opiates:
                f_opiates_theta13 = math.exp(theta_13 * age)
                f_opiates_theta11 = math.exp(theta_11 * age)
            else:
                f_opiates_theta13 = 1
                f_opiates_theta11 = 1

            # Volume of compartments [Litre]
            V1_arterial = theta_1 * f_centralWGT / f_centralWGTRef
            V1_venous = V1_arterial * (1 + theta_17 * (1 - f_centralWGT))
            self._v2comp3 = theta_2 * (weight / 70) * f_agingTheta10
            self._v3comp3 = theta_3 * (f_AlSallami / f_AlSallami_ref) * f_opiates_theta13

            # Elimination and inter - compartment clearance rates [Litre / min]
            self._cl1comp3 = theta_CL * (
                    weight / 70) ** 0.75 * f_CLmaturation * f_opiates_theta11 / f_CLmaturation_ref
            Q2_arterial = theta_5 * (self._v2comp3 / theta_2) ** 0.75 * (1 + theta_16 * (1 - f_Q3maturation))
            Q2_venous = Q2_arterial * theta_18
            self._cl3comp3 = theta_6 * (self._v3comp3 / theta_3) ** 0.75 * f_Q3maturation / f_Q3maturation_ref

            # The effect of the blood sampling site on the volume and clearance rates
            if blood_sampling == BloodSampling.ARTERIAL:
                self._v1comp3 = V1_arterial
                self._cl2comp3 = Q2_arterial
            elif blood_sampling == BloodSampling.VENOUS:
                self._v1comp3 = V1_venous
                self._cl2comp3 = Q2_venous
            else:
                raise ValueError(f"Blood sampling site {blood_sampling} not supported for Propofol")

        else:
            raise ValueError(f"Model {model} not supported for drug {Drug.PROPOFOL}")
        self.ss = self.getSS_3comp()

    def create_remifentanil(self, model: Model = Model.MINTO, **kwargs):
        """
        Create a pharmacokinetic model for remifentanil.
        :param model: Model enum
        :type model: Model
        :param kwargs: additional parameters for the model
        :type kwargs: dict
        """
        age = kwargs.get('age', -1)

        if age == -1:
            raise ValueError("Age is missing for Remifentanil PK model")

        if model is None or model == Model.MINTO:
            # Remifentanil parameters (Minto et al. 1997)
            lbm = kwargs.get('lbm', -1)
            if lbm == -1:
                raise ValueError("Lean body mass is missing for Minto PK model for Remifentanil")

            self._v1comp3 = 5.1 - 0.0201 * (age - 40) + 0.072 * (lbm - 55)
            self._v2comp3 = 9.82 - 0.0811 * (age - 40) + 0.108 * (lbm - 55)
            self._v3comp3 = 5.42
            self._cl1comp3 = 2.6 - 0.0162 * (age - 40) + 0.0191 * (lbm - 55)
            self._cl2comp3 = 2.05 - 0.0301 * (age - 40)
            self._cl3comp3 = 0.076 - 0.00113 * (age - 40)

        elif model == Model.ELEVELD:
            # Remifentanil parameters (Eleveld et al. 2017)
            weight = kwargs.get('weight', -1)
            gender = kwargs.get('gender', -1)
            bmi = kwargs.get('bmi', -1)

            if any(x == -1 for x in [weight, gender, bmi]):
                raise ValueError("Weight or gender or bmi is missing for Eleveld PK model for Remifentanil")

            v1_ref = 5.81  # [L]
            v2_ref = 8.82  # [L]
            v3_ref = 5.03  # [L]
            cl_ref = 2.58  # [L/min]
            q2_ref = 1.72  # [L/min]
            q3_ref = 0.124  # [L/min]
            theta1 = 2.88
            theta2 = -0.00554
            theta3 = -0.00327
            theta4 = -0.0315
            theta5 = 0.47
            theta6 = -0.0260

            if gender:  # female = 1
                ffm = (1.11 + (1 - 1.11) / (1 + (age / 7.1) ** (-1.1))) * (9270 * weight / (8780 + 244 * bmi))
            else:
                ffm = (0.88 + (1 - 0.88) / (1 + (age / 13.4) ** (-12.7))) * (9270 * weight / (6680 + 216 * bmi))

            ffm_ref = (1.11 + (1 - 1.11) / (1 + (35 / 7.1) ** (-1.1))) * (9270 * 70 / (8780 + 244 * 24.2))
            size_wg = ffm / ffm_ref

            kmat = Pharmacokinetic.f_sigmoid(weight, theta1, 2)
            kmat_ref = Pharmacokinetic.f_sigmoid(70, theta1, 2)
            if gender:
                ksex = 1 + theta5 * Pharmacokinetic.f_sigmoid(age, 12, 6) * (1 - Pharmacokinetic.f_sigmoid(age, 45, 6))
            else:
                ksex = 1

            self._v1comp3 = v1_ref * size_wg * Pharmacokinetic.f_aging(theta2, age)
            self._v2comp3 = v2_ref * size_wg * Pharmacokinetic.f_aging(theta3, age) * ksex
            self._v3comp3 = v3_ref * size_wg * Pharmacokinetic.f_aging(theta4, age) * np.exp(theta6 * weight - 70)
            self._cl1comp3 = cl_ref * size_wg ** 0.75 * (kmat / kmat_ref)
            self._cl2comp3 = q2_ref * (self._v2comp3 / v2_ref) ** 0.75 * Pharmacokinetic.f_aging(theta2, age) * ksex
            self._cl3comp3 = q3_ref * (self._v3comp3 / v3_ref) ** 0.75 * Pharmacokinetic.f_aging(theta2, age)

        else:
            raise ValueError(f"Model {model} not supported for drug {Drug.Remifentanil}")
        self.ss = self.getSS_3comp()

    def create_norepinephrine(self, model: Model = Model.JOACHIM, **kwargs):
        """
        Create a pharmacokinetic model for norepinephrine.
        :param model: Model enum
        :type model: Model
        :param kwargs: additional parameters for the model
        :type kwargs: dict
        """
        if model is None or model == Model.JOACHIM:
            # Norepinephrine parameters (Joachim et al. 2024)
            lbm = kwargs.get('lbm', -1)
            if lbm == -1:
                raise ValueError("Lean body mass is missing for Joachim PK model for Norepinephrine")
            self._k10 = np.exp(0.64 * np.log(lbm) - 5.52)  # [17s]
            self._ka = 0.02  # [1/s]
            self._k12 = 0.06  # [1/s]
            self._k21 = 0.04  # [1/s]
            self._v1comp3 = 0.49  # [L]
        else:
            raise ValueError(f"Model {model} not supported for drug {Drug.Norepinephrine}")
        self.ss = self.getSS_3comp()

    def create_rocuronium(self, model: Model = Model.DAHE, **kwargs):
        """
        Create a pharmacokinetic model for rocuronium.
        :param model: Model enum (Dahe)
        :type model: Model
        :param kwargs: additional parameters for the model
        :type kwargs: dict
        """
        if model is None or model == Model.DAHE:
            # Rocuronium parameters (Da Haes et al. 2002)
            weight = kwargs.get('weight', -1)
            if weight == -1:
                raise ValueError("Weight is missing for Dahe PK model for Rocuronium")

            self._v1comp3 = 42 * weight / 1000  # [L]
            self._v2comp3 = 40 * weight / 1000
            self._v3comp3 = 69 * weight / 1000
            self._cl1comp3 = 3.2 * weight / 1000  # [L/min]
            self._cl2comp3 = 5.2 * weight / 1000
            self._cl3comp3 = 0.9 * weight / 1000
        else:
            raise ValueError(f"Model {model} not supported for drug {Drug.Rocuronium}")
        self.ss = self.getSS_3comp()

    def getSS_3comp(self):
        """
        Get the state space model for a 3-compartment pharmacokinetic model.
        :return: The state space model for the 3-compartment model
        :rtype: signal.StateSpace
        """
        ss = None

        if self.drug in [Drug.PROPOFOL, Drug.REMIFENTANIL, Drug.ROCURONIUM]: # Three compartment model
            k10 = self._cl1comp3 / self._v1comp3
            k12 = self._cl2comp3 / self._v1comp3
            k13 = self._cl3comp3 / self._v1comp3
            k21 = self._cl2comp3 / self._v2comp3
            k31 = self._cl3comp3 / self._v3comp3

            A = np.array([[-(k10 + k12 + k13), k21, k31],
                          [k12, -k21, 0],
                          [k13, 0, -k31]]) / 60 # simulation is in seconds

            B = np.array([[1 / self._v1comp3],
                          [0],
                          [0]])

            C = np.array([[1, 0, 0]])
            D = np.array([[0]])
            ss = signal.StateSpace(A, B, C, D)

        elif self.drug == Drug.NOREPINEPHRINE:  # two compartment + depot PK model
            A = np.array([[-(self._k10 + self._k12), self._k21, self._ka],
                          [self._k12, -self._k21, 0],
                          [0, 0, -self._ka]])
            B = np.array([[0], [0], [1 / self._v1comp3]])
            C = np.array([[1, 0, 0]])
            D = np.array([[0]])

            ss = signal.StateSpace(A, B, C, D)

        return ss
