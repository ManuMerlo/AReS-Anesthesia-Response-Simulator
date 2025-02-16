from enum import Enum, auto, unique

# Types of disturbances considered in ARes
@unique
class DisturbanceType(Enum):
    INTUBATION = auto()
    INCISION = auto()
    SKIN_MANIPULATION = auto()
    SUTURE = auto()

# Types of interactions considered between the effects of remifentanil and propofol
@unique
class Interaction(Enum):
    NO_INTERACTION = "No_interaction"
    SURFACE = "Surface"

# Four drugs considered in ARes
@unique
class Drug(Enum):
    PROPOFOL = "Propofol"
    REMIFENTANIL = "Remifentanil"
    NOREPINEPHRINE = "Norepinephrine"
    ROCURONIUM = "Rocuronium"

# The names of the PK or PD model considered in ARes
@unique
class Model(Enum):
    ELEVELD = "Eleveld"
    MINTO = "Minto"
    SCHNIDER = "Schnider"
    JOACHIM = "Joachim"
    DAHE = "Dahe"
    PATIENT_SPECIFIC = "PatientSpecific"

# Types of outputs in AReS simulator
@unique
class Outputs(Enum):
    DOH = "DoH"
    HEMODYNAMICS = "Hemodynamics"
    NMB = "NMB"

# Types of depth of hypnosis measurement the user decide to simulate
@unique
class DoHMeasure(Enum):
    BIS = 0
    WAV = 1
    BOTH = 2

# The compartments that can be targeted when TCI is on
@unique
class TciMode(Enum):
    EFFECT_SITE = "EffectSite"
    PLASMA = "Plasma"

# Modes of ARes based on the inputs of the simulator
@unique
class SimulatorMode(Enum):
    INFUSION = "Infusion"
    CONCENTRATION = "Concentration"

# Hypovolemia indicates the status of fluid loss, and Normovolemia represents no fluid loss
@unique
class VolumeStatus(Enum):
    HYPOVOLEMIA = {"sv": 0.74, "co": 0.868, "map": 1.04, "hr": 1.16}
    NORMOVOLEMIA = {"sv": 1.0, "co": 1.0, "map": 1.0, "hr": 1.0}

# The anesthesia phase that the patient is under
@unique
class PatientPhase(Enum):
    INDUCTION = "Induction"
    MAINTENANCE = "Maintenance"
    EMERGENCE = "Emergence"

# Types of blood sampling used as a variable for Eleveld PK model
@unique
class BloodSampling(Enum):
    ARTERIAL = "Arterial"
    VENOUS = "Venous"
