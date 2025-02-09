# AReS - Anesthesia Response Simulator

This repository contains the source code for AReS (Anesthesia Response Simulator), a comprehensive anesthesia patient simulator designed to mimic the response of a patient undergoing total intravenous anesthesia (TIVA) with four commonly used drugs: propofol, remifentanil, norepinephrine, and rocuronium. The simulator is available in both Python and MATLAB versions.

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

## Getting Started

### Python Version

1. **Prerequisites**: Make sure you have Python 3.9 or higher installed.
2. **Installation**: Install the AReS package using `pip`:
   ```bash
   pip install AReS
   ```
3. **Usage**: Import the AReS module in your Python script and use its functions to create and run simulations. Refer to the documentation and examples in the Python folder for details.

### MATLAB Version

1. **Prerequisites**: Make sure you have MATLAB installed.
2. **Installation**: Add the MATLAB folder to your MATLAB path.
3. **Usage**: Use the MATLAB functions provided in the MATLAB folder to create and run simulations. Refer to the documentation and examples in the MATLAB folder for details.


## Examples

Example scripts demonstrating how to use AReS for different simulation scenarios are provided in notebooks folder.


## Contributing

Contributions to AReS are welcome! Please submit issues or pull requests on GitHub.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

[//]: # (## Citation)

[//]: # ()
[//]: # (If you use AReS in your research, please cite the associated manuscript.)