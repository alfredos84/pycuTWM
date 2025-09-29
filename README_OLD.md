# pycuTWM

*GPU-based Three-Wave Mixing processes in a nonlinear crystal using focused Gaussian beams in space-time domain.*

## Package Description

**pycuTWM** is a GPU-accelerated toolkit that simulates the coupled wave equations (CWEs) describing three-wave mixing (TWM) processes in second-order nonlinear media. The package is written to be as general as possible, including diffraction and dispersion effects within a single simulation. This means that the model is based on a (3+1)D physical problem (three spatial dimensions and one temporal dimension). The model incorporates terms for diffraction, dispersion, and linear absorption.

**pycuTWM** is built on a highly efficient C++/CUDA backend, with a Python wrapper that makes simulations user-friendly and easy to control.

The simulation outcomes are the pump, signal, and idler complex electric fields in the temporal, spatial, and spectral domains, as required.

C++/CUDA programming allows you to implement parallel computing in order to speed up calculations that typically require considerable computational demand. The Python wrapper makes using the package simple and accessible.

With this package, users can:
- Calculate the electric fields involved in second harmonic generation (SHG) and optical parametric generation (OPG) processes.
- Simulate both continuous-wave (cw) and pulsed pumping, using either focused Gaussian or plane-wave beams.
- For pulsed cases (femtosecond and picosecond regimes), the package efficiently simulates the simultaneous effects of dispersion and diffraction ((3+1)D problem).
- With the Python functions provided, users can easily run parameter sweeps over physical quantities of interest and perform simple postprocessing.

---

## Installation

To run this code, it is necessary to have an NVIDIA GPU installed on your computer, along with the CUDA drivers and the CUDA Toolkit.  
To install the CUDA driver and toolkit on your system, please visit the official documentation:  
https://docs.nvidia.com/cuda/

- Since the core engine uses a JSON configuration file, it is essential to install the library provided at [https://github.com/nlohmann/json](https://github.com/nlohmann/json).  
  **Before continuing with the installation, it is crucial to verify that the library can be included and used on your system.**  
  Please check that the developer's example works in your environment:
  ```cpp
  #include <nlohmann/json.hpp>
  ```
  If this test compiles and runs successfully, you are ready to proceed.

- Next, install the Python package:
  ```bash
  pip install pycutwm
  ```

- Finally, clone this repository to download all C++/CUDA files to your system:
  ```bash
  git clone https://github.com/alfredos84/pycuTWM.git
  ```

After completing these steps, you will have a directory with one main file and two subdirectories:

- `config.json`: This is the configuration file containing all parameters for the simulation.  
  The user can edit this file before running their scripts.  
  In the provided examples, you will notice that this file is copied with a new name and stored in a new folder. This copy will have the specific parameters of each simulation run.  
  This file is also used to run batch simulations, where a particular parameter (such as beam waist, pump power, etc.) is swept over multiple values.

- `src/`: This folder contains the main file `cuTWM.cu`, which runs the C++/CUDA core engine.  
  Inside this directory, you will also find the `headers/` folder with all relevant header files.

- `examples/`: In this folder, you will find several example and utility files to help you learn how to use this package.