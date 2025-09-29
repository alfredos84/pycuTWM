# Using HDF5 with CUDA/nvcc on Ubuntu (Universal Setup, Single Copy-Paste)

This guide shows how to install the HDF5 C++ library and configure your system so you can compile and run C++/CUDA code with HDF5 from any directory using `nvcc`, without needing to manually specify include/library paths each time.

---

## 1. Install the HDF5 Library

Open a terminal and run:
```bash
sudo apt-get update
sudo apt-get install libhdf5-dev
```
---

## 2. Identify HDF5 Include and Library Paths

Typical Ubuntu/Debian paths are:

- Include: /usr/include/hdf5/serial
- Libs: /usr/lib/x86_64-linux-gnu/hdf5/serial

You can verify with:
```bash
dpkg -L libhdf5-dev | grep include
dpkg -L libhdf5-dev | grep lib
```
---

## 3. Add HDF5 Paths to Environment Variables

Add the following lines to your ~/.bashrc (or ~/.zshrc if you use zsh):

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

## 4. Minimal Example: C++ Code Using HDF5 (works with nvcc)

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

## 5. Compile With nvcc (no need to specify paths)

With the environment variables set, simply run:
```bash
nvcc test_h5_nvcc.cpp -o test_h5_nvcc -lhdf5_cpp -lhdf5
```
---

## 6. Run the Example
```
./test_h5_nvcc
```
Expected output:

File created: test_file.h5
Value read: 3.14159

---

Now you can compile HDF5 C++ code with CUDA/nvcc from any folder without worrying about include/library paths!
