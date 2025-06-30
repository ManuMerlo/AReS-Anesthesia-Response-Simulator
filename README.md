# AReS - Anesthesia Response Simulator

This repository contains the source code for AReS (Anesthesia Response Simulator), a comprehensive anesthesia patient simulator designed to mimic the response of a patient undergoing total intravenous anesthesia (TIVA) with four commonly used drugs: propofol, remifentanil, norepinephrine, and rocuronium. The simulator is available in both Python and MATLAB versions.
The methods for developing AReS and several analysis to examine its performance are provided in the following reference: 

Hosseinirad, Sara and Merlo, Manuela and Trovò, Francesco and Tognoli, Emiliano and Metelli, Alberto Maria and Dumont, Guy A., Ares: A Patient Simulator to Facilitate Testing of Automated Anesthesia. Available at Journal of Computer Methods and Programs in Biomedicine: https://doi.org/10.1016/j.cmpb.2025.108901
## Features

* **Pharmacokinetic/Pharmacodynamic (PK/PD) Modeling**: AReS utilizes population-based and subject-specific PK/PD models to predict drug effects on various physiological parameters.
* **Target-Controlled Infusion (TCI)**: The simulator incorporates a TCI module to control drug concentrations at the target site (plasma or effect-site).
* **Surgical Stimuli**: AReS can simulate various surgical stimuli with varying dynamics and durations, allowing for the investigation of the patient's response to surgical stress.
* **Fluid Loss**: The simulator can also model the effect of fluid loss on the patient's physiological response.
* **Output Variables**: AReS provides predictions for a range of vital physiological parameters, including:
    * Depth of hypnosis (NeuroWave index and Bispectral Index)
    * Cardiac output
    * Mean arterial pressure
    * Heart rate
    * Stroke volume
    * Neuromuscular blockade

## Repository Structure

* **Python**: Source code for the Python version of AReS, implemented as a Python package.
* **Matlab**: Source code for the MATLAB version of AReS.


## Contributing

Contributions to AReS are welcome! Please submit issues or pull requests on GitHub.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

[//]: # (## Citation)

[//]: # ()
[//]: # (If you use AReS in your research, please cite the associated manuscript.)