# pycuTWM: High-Performance Three-Wave Mixing Simulation Package for simulations in (3+1)D

## Overview
**pycuTWM** is a GPU-accelerated simulation suite for modeling nonlinear three-wave mixing (TWM) processes in optical crystals. The core engine is written in `C++/CUDA`, using the CUDA libraries `cuFFT` (the standard `fftw` Fast-Fourier Transform library for CUDA), `THRUST` (the Standard Template Library (STL) for CUDA). `C++/CUDA` programming allows you to implement parallel computing in order to speed up calculations that typically require considerable computational demand. 

Input data are managed by a single JSON file (`nlohmann/json` is required), while output are stored in`HDF5` files for efficient computation and flexible configuration and portability. 

A Python wrapper makes using the package simple and accessible. This Python package (`pycutwm`) is provided for seamless workflow: compilation, execution, parameter sweeps, and a basic post-processing. 

The simulation outcomes are the pump, signal, and idler complex electric fields in the temporal and spatial domains ($x,y,z,t$ coordinates). The output are stored in HDF5 files for efficient and flexible portability.


## Package Description

**pycuTWM** is a GPU-accelerated toolkit that simulates the coupled wave equations (CWEs) describing three-wave mixing (TWM) processes in second-order nonlinear media. The physics solved in the package is written to be as general as possible, including diffraction and dispersion effects within a single simulation. This means that the model is based on a (3+1)D physical problem (three spatial dimensions and one temporal dimension). The model incorporates terms for diffraction, dispersion, and linear absorption.

With this package, users can:
- Calculate the (3+1)D-electric fields involved, $A_{\lambda} = A_{\lambda}(x,y,z,t)$, in Sum Frequency Generation (SFG) (keep in mind that second harmonic generation (SHG) is a particular case of SFG) and optical parametric generation (OPG) processes.

$$ \frac{\partial A_{p}}{\partial z} = i\kappa_p A_{s} A_{i}e^{\mp i\Delta k z} + \left(\hat{\mathcal{D}}^{(\tau)}_{p}+\hat{\mathcal{D}}^{(xy)}_{p} - \frac{\alpha_p}{2} \right)A_{p} $$
$$ \frac{\partial A_{s}}{\partial z} = i\kappa_s A_{p} A_{i}^*e^{\pm i\Delta k z} + \left(\hat{\mathcal{D}}^{(\tau)}_{s}+\hat{\mathcal{D}}^{(xy)}_{s} - \frac{\alpha_s}{2} \right)A_{s} $$
$$ \frac{\partial A_{i}}{\partial z} = i\kappa_i A_{p} A_{s}^*e^{\pm i\Delta k z} + \left(\hat{\mathcal{D}}^{(\tau)}_{i}+\hat{\mathcal{D}}^{(xy)}_{i} - \frac{\alpha_i}{2} \right)A_{i}$$
where
$$\hat{\mathcal{D}}^{(\tau)}_{\lambda} = -\left[ \frac{\alpha_{\lambda}}{2}+ \left(\frac{1}{\nu_s} - \frac{1}{\nu_{\lambda}}\right) \frac{\partial}{\partial \tau}+i\frac{k^{''}_{\lambda}}{2}\frac{\partial^2}{\partial \tau^2} + i\frac{k^{'''}_{\lambda}}{3}\frac{\partial^3}{\partial \tau^3} \right]$$
and
$$  \hat{\mathcal{D}}^{(xy)}_{\lambda} = - \frac{i}{2k_{\lambda}}\left(\frac{\partial^2}{\partial x^2}+\frac{\partial^2}{\partial y^2}\right) $$
- Simulate both continuous-wave (cw) and pulsed pumping, using either focused Gaussian or plane-wave beams.
- For pulsed cases (femtosecond and picosecond regimes), the package efficiently simulates the simultaneous effects of dispersion and diffraction ((3+1)D problem).
- With the Python functions provided, users can easily run parameter sweeps over physical quantities of interest and perform simple postprocessing.

---

## Installation

### System Requirements

- An **NVIDIA GPU** with supported drivers already installed. Check with `nvidia-smi`: it should list your GPU.
- At least 8 GB RAM recommended
- The package was tested on Linux-Ubuntu
- `gcc`,`g++`, and `nvcc` compilers
- Libraries `JSON` and H5 for `C++`
- The package requires at least `C++14`

## 1. Using the nlohmann/json C++ Library with Any C++/CUDA Project

This guide explains how to install the [nlohmann/json](https://github.com/nlohmann/json) C++ library (modern, header-only JSON for C++), and configure your system so you can compile and run C++ (or CUDA) code using JSON from any directory, **without needing to manually specify include paths each time**.

---

### a. Install the nlohmann/json Library (Ubuntu/Debian)

The package is header-only and available as `nlohmann-json-dev` in Ubuntu repos:
```bash
sudo apt-get update
sudo apt-get install nlohmann-json3-dev
```
---

### b. Identify JSON Include Path

Typical location is:

- /usr/include/nlohmann/json.hpp

You can check with:
```bash
dpkg -L nlohmann-json3-dev | grep json.hpp
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

### e. Compile With g++ or nvcc

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

```bash
./test_json_nvcc 
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

Add the following lines to your `~/.bashrc`:

```bash
export CPATH=/usr/include/hdf5/serial:$CPATH
export LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial:$LIBRARY_PATH
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH
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
g++ test_h5_nvcc.cpp -o test_h5_nvcc -lhdf5_cpp -lhdf5
```
---

### f. Run the Example

```bash
./test_h5_nvcc
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

- Prerequisites: [Miniconda or Anaconda](https://docs.conda.io/en/latest/miniconda.html) installed.

1. Create the environment

Open a terminal and run:
  ```bash
    conda create -n pytwm-env python=3.8 # or higher
    conda activate pytwm-env
    conda install -c conda-forge libstdcxx-ng
  ```
2. Install all dependencies from conda-forge
  ```bash
    conda activate pytwm-env
    pip install pycutwm
  ```
- **Please check the suitable** `cudatoolkit=XX` for your system.
- Package `pycutwm` includes a set functions to execute the pyTWM package.
- You can add additional Python packages as needed (e.g. numpy, scipy, h5py).

3. Compiling your CUDA/C++ code inside conda

**Always** compile using the tools and libraries provided by the conda environment.

Execute the test example included in the repository (folder `examples/test_example`):
  ```bash
    conda activate pytwm-env
    python test_example.py
  ```
If you see on screen **Compilation succeeded.**, congratulations!.

4. Included examples

In the `examples/` folder, users will find several ready-to-run illustrative examples. These examples are designed to help you quickly reproduce typical experiments and to gain a deeper understanding of the underlying physics beyond what measurements alone can reveal.

I encourage curious users to explore this folder and discover the full range of possibilities offered by `pycuTWM`.
Please do not hesitate to contact me if you have any questions or suggestions!

5. Best practices

- **Never mix system and conda libraries:** Always compile and run inside your conda environment.
- **Document your CUDA architecture flags** (`compute_XX`, `sm_XX`) are authomatically set.
- **Ensure drivers are installed** (users must have an NVIDIA driver compatible with the CUDA version).

6. Troubleshooting

- **`nvidia-smi` does not show your GPU:**  
  (Re)install the correct NVIDIA driver for your operating system.
- **Compiler errors (`gcc`):**  
  Conda provides its own toolchains (`gcc_linux-64`, `gxx_linux-64`).

Summary

This approach ensures:
- No root/admin required
- Seamless Python/CUDA/C++ integration
- Easy environment sharing and reproducibility

For advanced examples, further support, or questions, open an issue on the project repository.

## License

MIT Licence
