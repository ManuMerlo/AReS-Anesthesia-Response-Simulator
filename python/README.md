# AReS – Anesthesia Response Simulator

**AReS** (Anesthesia Response Simulator) is an open-source tool designed to simulate patient responses under anesthesia. 

You can use AReS as a standalone Python package for simulations or explore and modify it locally.

---

## 📦 Option 1: Install AReS as a Python Package

Use this option if you want to **integrate AReS in your own Python projects** without modifying the source code.

### Prerequisites

* Python 3.9 or higher
* Git

### Installation

```bash
pip install git+https://github.com/ManuMerlo/AReS-Anesthesia-Response-Simulator.git@main#subdirectory=python
```

### Usage

Import the package in your Python scripts and call its functions as needed.
Example notebooks demonstrating usage are available in the `python/notebooks/` folder of this repository.

---
## 🔧 Option 2: Local Development Setup

Use this option to **run notebooks and explore or modify the source code** in a local environment.

### Prerequisites

* [Conda (Miniconda or Anaconda)](https://docs.conda.io/en/latest/miniconda.html)
* Git

---

## 🖥️ Terminal-Based Workflow

Use this workflow if you prefer working in a terminal or command line interface.

### 1. Clone the Repository

```bash
git clone https://github.com/ManuMerlo/AReS-Anesthesia-Response-Simulator.git
cd AReS-Anesthesia-Response-Simulator/python
```

### 2. Create the Conda Environment

```bash
conda env create -f environment.yml
```

### 3. Activate the Environment

```bash
conda activate AReS
```
Open and run notebooks located in the `notebooks/` folder (note: some notebooks are outside the `python` folder).

---

## 🖥️ IDE-Based Workflow

If you prefer to work inside an Integrated Development Environment (IDE), you can run everything there — notebooks and code — by configuring your IDE to use the Conda environment.

### 1. Clone the Repository and Open the Project

#### <img src="https://code.visualstudio.com/assets/favicon.ico" alt="VS Code" width="20"/>  **VS Code**

1. Open the **Command Palette** (`Ctrl+Shift+P` or `Cmd+Shift+P`)
2. Type and select `Git: Clone`
3. Paste the repo URL: `https://github.com/ManuMerlo/AReS-Anesthesia-Response-Simulator.git`
4. Choose a destination folder
5. Once cloned, open the **root folder**:
   `AReS-Anesthesia-Response-Simulator`

---

#### <img src="https://resources.jetbrains.com/storage/products/pycharm/img/meta/pycharm_logo_300x300.png" alt="PyCharm" width="20"/>  **PyCharm**

1. On the **Welcome screen**, select **Clone Repository**
2. Paste the repo URL: `https://github.com/ManuMerlo/AReS-Anesthesia-Response-Simulator.git`
3. Choose your destination folder and clone
4. Open the **root folder**:
   `AReS-Anesthesia-Response-Simulator`

---

### 2. Create the Conda Environment

* Open the terminal inside your IDE or a system terminal.

* Navigate into the `python` folder (where `environment.yml` is located):

  ```bash
  cd python
  ```

* Create the Conda environment:

  ```bash
  conda env create -f environment.yml
  ```

* Activate the environment:

  ```bash
  conda activate AReS
  ```
---

### 3. Configure the IDE to Use the Conda Environment

#### <img src="https://code.visualstudio.com/assets/favicon.ico" alt="VS Code" width="20"/> **VS Code**

* Open the **Command Palette** (`Ctrl+Shift+P` or `Cmd+Shift+P`)
* Select `Python: Select Interpreter`
* Choose the Conda environment named `AReS`
* You can now run notebooks and scripts directly from the root folder or any subfolder

---

#### <img src="https://resources.jetbrains.com/storage/products/pycharm/img/meta/pycharm_logo_300x300.png" alt="PyCharm" width="20"/> **PyCharm**

* Navigate to **Settings > Project > Python Interpreter**
* Click **Add Interpreter** → **Add Local Interpreter** → **Conda Environment**
* Select **Use existing environment** and choose `AReS`
* ✅ PyCharm **Professional** supports running Jupyter notebooks natively

---

