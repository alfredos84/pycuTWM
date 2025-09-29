# pycuTWM: High-Performance Three-Wave Mixing Simulation Package for simulations in (3+1)D

## Overview

pycuTWM is a GPU-accelerated simulation suite for modeling nonlinear three-wave mixing (TWM) processes in optical crystals. The core engine is written in C++/CUDA, using cuFFT, THRUST, nlohmann/json, and HDF5 for efficient computation and flexible configuration. A Python package (`pycutwm`) is provided for seamless workflow: compilation, execution, parameter sweeps, and a basic post-processing.

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

### System Requirements

- An **NVIDIA GPU** with supported drivers already installed. Check with `nvidia-smi`: it should list your GPU.
- At least 8 GB RAM recommended
- Linux-Ubuntu or Windows 10/11 (x64)
- GCC (Linux-Ubuntu) or MSVC (Windows)
- Libraries `JSON` and H5 for `C++`

## 1. Using the nlohmann/json C++ Library with Any C++/CUDA Project

This guide explains how to install the [nlohmann/json](https://github.com/nlohmann/json) C++ library (modern, header-only JSON for C++), and configure your system so you can compile and run C++ (or CUDA) code using JSON from any directory, **without needing to manually specify include paths each time**.

---

### a. Install the nlohmann/json Library (Ubuntu/Debian)

The package is header-only and available as `nlohmann-json-dev` in Ubuntu repos:
```bash
sudo apt-get update
sudo apt-get install nlohmann-json-dev
```
---

### b. Identify JSON Include Path

Typical location is (**Linux-Ubuntu**):

- /usr/include/nlohmann/json.hpp

You can check with:
```bash
dpkg -L nlohmann-json-dev | grep json.hpp
```
---

### c. Add JSON Include Path to Environment Variable

Add this to your `~/.bashrc` (**Linux-Ubuntu**):
```bash
export CPATH=/usr/include/nlohmann:$CPATH
```
Then reload your shell:
```bash
source ~/.bashrc
```
---

### d. Minimal Example: C++ Code Using nlohmann/json (works with nvcc)

Create a file called `test_json_nvcc.cpp`:
```cpp
#include <iostream>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

int main() {
    // Create JSON object
    json j;
    j["pi"] = 3.14159;
    j["answer"]["everything"] = 42;
    j["name"] = "Alfredo";
    j["is_cuda"] = true;

    // Serialize to string and print
    std::cout << j.dump(4) << std::endl;

    // Parse from string
    std::string s = R"({"greeting":"Hello, world!","year":2025})";
    json j2 = json::parse(s);
    std::cout << "Greeting: " << j2["greeting"] << ", Year: " << j2["year"] << std::endl;

    return 0;
}
```
---

### e. Compile With g++ or nvcc (no need to specify -I paths)

With the environment variable set, just run:
```bash
g++ test_json_nvcc.cpp -o test_json_nvcc
```
Or, if you want to check nvcc works (for host code):
```bash
nvcc test_json_nvcc.cpp -o test_json_nvcc
```
---

### f. Run the Example
**Linux-Ubuntu**
```bash
./test_json_nvcc 
```
or 
**Windows**
```bash
./test_json_nvcc.exe
```

Expected output:
```bash
{
    "answer": {
        "everything": 42
    },
    "is_cuda": true,
    "name": "Alfredo",
    "pi": 3.14159
}
Greeting: Hello, world!, Year: 2025
```
---

Now you can compile C++/CUDA code with nlohmann/json from any folder without worrying about include paths!



## 2. Using HDF5 with CUDA/nvcc on Ubuntu

This guide shows how to install the HDF5 C++ library and configure your system so you can compile and run C++/CUDA code with HDF5 from any directory using `nvcc`, without needing to manually specify include/library paths each time.

---

### a. Install the HDF5 Library

Open a terminal and run:
```bash
sudo apt-get update
sudo apt-get install libhdf5-dev
```
---

### b. Identify HDF5 Include and Library Paths

Typical Ubuntu/Debian paths are:

- Include: /usr/include/hdf5/serial
- Libs: /usr/lib/x86_64-linux-gnu/hdf5/serial

You can verify with:
```bash
dpkg -L libhdf5-dev | grep include
dpkg -L libhdf5-dev | grep lib
```
---

### c. Add HDF5 Paths to Environment Variables

Add the following lines to your ~/.bashrc:

```bash
export CPATH=/usr/include/hdf5/serial:$CPATH
export LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial:$LIBRARY_PATH
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial:$LD_LIBRARY_PATH
```
Then reload your shell:
```bash
source ~/.bashrc
```
---

### d. Minimal Example: C++ Code Using HDF5 (works with nvcc)

Create a file called `test_h5_nvcc.cpp`:

```cpp
#include <iostream>
#include <H5Cpp.h>
using namespace H5;

const H5std_string FILE_NAME("test_file.h5");
const H5std_string DATASET_NAME("scalar");

int main() {
    double value_out = 3.14159;
    try {
        H5File file(FILE_NAME, H5F_ACC_TRUNC);
        hsize_t dim[] = {1};
        DataSpace dataspace(1, dim);
        DataSet dataset = file.createDataSet(DATASET_NAME, PredType::NATIVE_DOUBLE, dataspace);
        dataset.write(&value_out, PredType::NATIVE_DOUBLE);
        file.close();
        std::cout << "File created: " << FILE_NAME << std::endl;
    } catch (FileIException& error) { error.printErrorStack(); return -1; }
    catch (DataSetIException& error) { error.printErrorStack(); return -1; }

    double value_in = 0.0;
    try {
        H5File file(FILE_NAME, H5F_ACC_RDONLY);
        DataSet dataset = file.openDataSet(DATASET_NAME);
        dataset.read(&value_in, PredType::NATIVE_DOUBLE);
        std::cout << "Value read: " << value_in << std::endl;
        file.close();
    } catch (FileIException& error) { error.printErrorStack(); return -1; }
    catch (DataSetIException& error) { error.printErrorStack(); return -1; }
    return 0;
}
```
---

### e. Compile With nvcc (no need to specify paths)

With the environment variables set, simply run:
```bash
nvcc test_h5_nvcc.cpp -o test_h5_nvcc -lhdf5_cpp -lhdf5
```
---

### f. Run the Example
**Linux-Ubuntu**
```bash
./test_h5_nvcc
```
or
**Windows**
```bash
./test_h5_nvcc.exe
```
Expected output:

File created: test_file.h5
Value read: 3.14159

---

Now you can compile HDF5 C++ code with CUDA/nvcc from any folder without worrying about include/library paths!



## 3. Clone the repository to download all C++/CUDA files to your system:
  ```bash
  git clone https://github.com/alfredos84/pycuTWM.git
  ```

## 4. Create an isolated Conda environment: `pytwm-env`

This guide explains how to set up a **fully reproducible, portable environment** for compiling and running CUDA/C++ scientific software with Python wrappers using [conda](https://docs.conda.io/).

#### Why use conda?

- **Works on Linux-Ubuntu and Windows**.
- No need for root/administrator permissions.
- All dependencies (CUDA Toolkit, nvcc, HDF5, cuFFT, Thrust, FFTW, curl, Python libraries, C++ toolchain) are managed inside the environment.
- Maximum reproducibility: users can clone your environment exactly.
- Prerequisites: [Miniconda or Anaconda](https://docs.conda.io/en/latest/miniconda.html) installed.

1. Create the environment

Open a terminal (**Linux-Ubuntu**) or Anaconda Prompt (**Windows**) and run:
  ```bash
    conda create -n pytwm-env python=3.12
    conda activate pytwm-env
  ```
2. Install all dependencies from conda-forge
  ```bash
    conda install -c conda-forge cudatoolkit=12.4 cudatoolkit-dev hdf5 libcurl fftw zlib gxx_linux-64 gcc_linux-64
    conda activate pytwm-env
    pip install pycutwm
  ```
- **Please check the suitable** `cudatoolkit=XX` for your system.
- Package `pycutwm` includes a set functions to execute the pyTWM package.
- On **Windows**, use the same command (conda will pick the appropriate toolchains for your OS).
- You can add additional Python packages as needed (e.g. numpy, scipy, h5py).

3. Compiling your CUDA/C++ code inside conda

**Always** compile using the tools and libraries provided by the conda environment.

Execute the test example included in the repository:
  ```bash
    conda activate pytwm-env
    python test_example.py
  ```
If you see on screen **Compilation succeeded.**, congratulations!.

4. Best practices

- **Never mix system and conda libraries:** Always compile and run inside your conda environment.
- **Document your CUDA architecture flags** (`compute_XX`, `sm_XX`) are authomatically set.
- **Ensure drivers are installed** (users must have an NVIDIA driver compatible with the CUDA version).
- **Windows:** Adapt include/lib paths and linker flags as needed.

5. Troubleshooting

- **Missing `.so` or `.dll` errors:**  
  Ensure `LD_LIBRARY_PATH` (Linux-Ubuntu) or the equivalent (Windows) is set as above.
- **`nvidia-smi` does not show your GPU:**  
  (Re)install the correct NVIDIA driver for your operating system.
- **Compiler errors (`gcc`):**  
  Conda provides its own toolchains (`gcc_linux-64`, `gxx_linux-64`).

Summary

This approach ensures:
- Maximum portability (Linux-Ubuntu and Windows)
- No root/admin required
- Seamless Python/CUDA/C++ integration
- Easy environment sharing and reproducibility

For advanced examples, further support, or questions, open an issue on the project repository.

## License

[Specify your license, e.g., MIT, GPL-3.0, etc.]
