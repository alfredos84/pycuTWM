# pytwm-env: Universal CUDA/Python Scientific Environment

This guide explains how to set up a **fully reproducible, portable environment** for compiling and running CUDA/C++ scientific software with Python wrappers using [conda](https://docs.conda.io/).

## Why use conda?

- **Works on Linux and Windows**.
- No need for root/administrator permissions.
- All dependencies (CUDA Toolkit, nvcc, HDF5, cuFFT, Thrust, FFTW, curl, Python libraries, C++ toolchain) are managed inside the environment.
- Maximum reproducibility: users can clone your environment exactly.

## Prerequisites

- An **NVIDIA GPU** with supported drivers already installed.  
  (Check with `nvidia-smi` — it should list your GPU.)
- [Miniconda or Anaconda](https://docs.conda.io/en/latest/miniconda.html) installed.

## 1. Create the environment

Open a terminal (**Linux**) or Anaconda Prompt (**Windows**) and run:

    conda create -n pytwm-env python=3.12
    conda activate pytwm-env

## 2. Install all dependencies from conda-forge

    conda install -c conda-forge cudatoolkit=XX cudatoolkit-dev hdf5 libcurl fftw zlib gxx_linux-64 gcc_linux-64
    conda activate pytwm-env
    pip install pycutwm

- Please check the suitable `cudatoolkit=XX` for your system.
- Package `pycutwm` includes a set functions to execute the pyTWM package.
- On **Windows**, use the same command (conda will pick the appropriate toolchains for your OS).
- You can add additional Python packages as needed (e.g. numpy, scipy, h5py).

## 3. Compiling your CUDA/C++ code inside conda

**Always** compile using the tools and libraries provided by the conda environment.

Execute the test example included in the repository:

    conda activate pytwm-env
    python test_example.py

If you see on screen **Hello pyTWM**, congratulations!.

## 4. Best practices

- **Never mix system and conda libraries:** Always compile and run inside your conda environment.
- **Document your CUDA architecture flags** (`compute_XX`, `sm_XX`) are authomatically set.
- **Ensure drivers are installed** (users must have an NVIDIA driver compatible with the CUDA version).
- **Windows:** Adapt include/lib paths and linker flags as needed.

## 5. Troubleshooting

- **Missing `.so` or `.dll` errors:**  
  Ensure `LD_LIBRARY_PATH` (Linux) or the equivalent (Windows) is set as above.
- **`nvidia-smi` does not show your GPU:**  
  (Re)install the correct NVIDIA driver for your operating system.
- **Compiler errors (`gcc`):**  
  Conda provides its own toolchains (`gcc_linux-64`, `gxx_linux-64`).

## Summary

This approach ensures:
- Maximum portability (Linux and Windows)
- No root/admin required
- Seamless Python/CUDA/C++ integration
- Easy environment sharing and reproducibility

For advanced examples, further support, or questions, open an issue on the project repository.