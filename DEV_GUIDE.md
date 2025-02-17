# Developer Guide

## Overall code structure

The software is written as a relatively lightweight python package. It is intended to be used in python scripts.
There is currently a strong dependence on RDKit which is used for reading and writing files, and various molecular operations,
e.g. assigning pharmacophore features. 

The Roshambo2 code treats molecules as a set of 3D points, it does not know about molecular graphs etc.

The code has a Python front end which has the following functionalities:
- interface with RDkit for SDF I/O.
- interface with RDkit for feature assignment.
- uses h5py for Roshambo2 HDF5 file format I/O.
- initial 3D preparation stage (align to center of mass and PCA axes).
- data manipulation from lists of 3D molecules to batched numpy arrays of many 3D configurations.
- passing numpy arrays to the C++/CUDA backend code.
- sorting the results and data post-processing (e.g. alignment of best fit structures for output).

The code has a C++ and CUDA backend which has the following functionalities:
- uses the gaussian overlap method to find the optimal rotation and translation to best align two 3D configurations. (many in parallel)
- compute the shape and color overlap volumes and the corresponding tanimoto scores of the best align configurations.
- takes the input molecule coordinates as numpy arrays, returns the scores and best fit transformation as a numpy array. 

## Python code 

The API documentation lists all of the standard user facing classes and functions.
There are other modules and code that are not documented there. 
These are the internal functions that are not designed to by directly used by a typical user.

Currently there are:

- roshambo2/utils.py: Contains various python utility functions
- roshambo2/backends/_pytorch_backend.py: The PyTorch implementation of Gaussian shape overlay. Used in the testing framework.
- roshambo2/backends/cpp_backend.py:  The python wrapper around the cpp code.
- roshambo2/backends/cuda_backend.py: The python wrapper around the cuda code.


## PyTorch code

There is a PyTorch backend which can be used (for shape only and start_mode=0 only):

```
roshambo2_calculator.compute(backend="_pytorch") 
```

It is easy to modiy the optimizer settings with this backend, e.g. in `roshambo2/backends/_pytorch_backend.py`
you can change this line to a different optimizer:
```
    optimizer = torch.optim.Adagrad([{'params': q, 'lr': 0.1}, {'params': t, 'lr': 0.1}], )
```

The pytorch backend was implemented first to be used as the testing frame work for ensuring the gradients and optimization are correct in the C++/CUDA code.

## C++

The C++ backend is all in the file `roshambo2/backends/cpp_src/cpp_functions.cpp`.

The main function is `optimize_overlap_color(..)`. This is the only function that is called by the main roshambo2 python code. It is wrapped with pybind. The input arguments are numpy arrays, ints, or floats.
It is important that the dtype of the numpy array matches the dtype of the C++ code. `dtype=np.float32` should be used in the numpy side and `float` should be used in the C++ side.

A few other functions are wrapped with pybind for the testing framework. The binding of functions for use in python is done at the end of the file:
```c++
////////////////////////////////////////////////////////////////////////////////
/// Bindings for Python
////////////////////////////////////////////////////////////////////////////////
PYBIND11_MODULE(_roshambo2_cpp, m) { 
    m.def("optimize_overlap_color", &optimize_overlap_color, "computes overlap of ref mol A with fit mols B with color");
    m.def("test_overlap", &test_overlap, "computes overlap of ref mol A with fit mols B");
    m.def("test_overlap_single", &test_overlap_single, "single overlap for testing");
    m.def("test_gradient", &test_gradient, "quaternion gradient of single overlap for testing");
}
```

The c++ program flow is:

- (loop 1) for each query config:
    - (loop 2) over each start mode
        - (loop 3) over each dataset config
            - (loop 4) run N steps of optimization
    
    - save the results for the best start mode

The loop number 3, over the dataset configs, is parallel with OpenMP threads.


## CUDA

The CUDA backend is in `roshambo2/backends/cuda_src`. There are currently three files:
- roshambo2/backends/cuda_src/cuda_interface.cpp: contains the c++ function that is wrapped with pybind to be used in the python code. The function in this file has the same arguments and name as the `optimize_overlap_color(..)` in the c++ backend. The function copies the data to the GPU and launches the cuda code.
- `roshambo2/backends/cuda_src/cuda_functions.cuh`. Headers and constants for the cuda code.
- `roshambo2/backends/cuda_src/cuda_functions.cu`. Cuda kernel code.

The functions are all very similar to the c++ ones.

The program flow is:

- (loop1) for N GPUs split the dataset configs into N chunks
    - (loop2) for each query config:
        - (kernel launch) run each query-dataset overlap on a seperate thread
            - (per thread):
                - loop over each start mode
                - run N steps of optimization
                - save the results for best start mode


## PyBind11

PyBind11 (https://pybind11.readthedocs.io/en/stable/) is used to create the python bindings of the c++ code.

Specifically the numpy direct access mode is used https://pybind11.readthedocs.io/en/stable/advanced/pycpp/numpy.html#direct-access



## Build process

The c++ and CUDA code is compiled by CMAKE. When you do `pip install .` the setup.py script contains the commands to
launch the cmake process with the correct options.

If you are developing the c++/cuda code it can be easier if you do the cmake part yourself.

### Manual cmake build
from inside the top level roshambo2 folder run these commands:
```bash
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../ # this will install to the same location that pip install -e . does, you can also do -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX
# you can check the options with the gui ccmake (if installed)
ccmake ..
# If everything looks good you can build and install with
make install
```

If the configure steps fail it usually means you need to set some more paths to the CUDA libraries.
If the make step fails it usually means there are errors in the code, or that the linker paths have not been set correctly





