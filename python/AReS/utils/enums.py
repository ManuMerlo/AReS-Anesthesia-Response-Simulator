from enum import Enum, auto, unique

@unique
class DisturbanceType(Enum):
    INTUBATION = auto()
    INCISION = auto()
    SKIN_MANIPULATION = auto()
    SUTURE = auto()


@unique
class Interaction(Enum):
    NO_INTERACTION = "No_interaction"
    SURFACE = "Surface"


@unique
class Drug(Enum):
    PROPOFOL = "Propofol"
    REMIFENTANIL = "Remifentanil"
    NOREPINEPHRINE = "Norepinephrine"
    ROCURONIUM = "Rocuronium"


@unique
class Model(Enum):
    ELEVELD = "Eleveld"
    MINTO = "Minto"
    SCHNIDER = "Schnider"
    JOACHIM = "Joachim"
    DAHE = "Dahe"
    PATIENT_SPECIFIC = "PatientSpecific"


@unique
class Outputs(Enum):
    DOH = "DoH"
    HEMODYNAMICS = "Hemodynamics"
    NMB = "NMB"


@unique
class DoHMeasure(Enum):
    BIS = 0
    WAV = 1
    BOTH = 2


@unique
class TciMode(Enum):
    EFFECT_SITE = "EffectSite"
    PLASMA = "Plasma"


@unique
class SimulatorMode(Enum):
    INFUSION = "Infusion"
    CONCENTRATION = "Concentration"


@unique
class VolumeStatus(Enum):
    HYPOVOLEMIA = {"sv": 0.74, "co": 0.868, "map": 1.04, "hr": 1.16}
    NORMOVOLEMIA = {"sv": 1.0, "co": 1.0, "map": 1.0, "hr": 1.0}


@unique
class PatientPhase(Enum):
    INDUCTION = "Induction"
    MAINTENANCE = "Maintenance"
    EMERGENCE = "Emergence"

@unique
class BloodSampling(Enum):
    ARTERIAL = "Arterial"
    VENOUS = "Venous"
