# ROSHAMBO 2.0

Python program for molecular shape overlay calculation.

Uses the Gaussian description of molecular shape (Grant and Pickup 1995 https://doi.org/10.1021/j100011a016)

MIT License. Copyright (c) 2025

## Install

This software is distributed in source code form. It can be build with pip inside an appropriate conda environment. Please follow the steps below.

1. clone this repo or (extract from tar archive/zip folder).

```
git clone <repo>
```

2. create a conda environment with the provided yaml file

```
cd roshambo2
conda env create -n roshambo2 -f environment.yaml
conda activate roshambo2
```

3. install this package

```
pip install .
```

## Documentation install

To build and view the documentation please go into the `doc` folder.
You will then need to install the extra dependencies to your conda environment:

```
conda install  -c sphinx pydata-sphinx-theme myst-parser
```

You can then build the html docs with:

```
make html
```

And you can view them in your browser:

```
xdg-open build/html/index.html
```

## Requirements

A CMake build system is used to compile the C++ and CUDA code. You will need to have a working C++ and nvcc compiler. Cmake will need to be able to find your CUDA installation.

## Usage

See the examples folder: "examples/README.md"

And the user guide: "USER_GUIDE.md"

## Developers

See the developer guide: "DEV_GUIDE.md"

## Testing

Test can be run using pytest in the test folder. You will need to install `pytest` and `pytorch`.

```
cd test
pytest
```

Note that the PyTorch is used in the testing framework.
(Autograd can compute the overlap gradients which are compared with
the hand-coded ones in the ROSHAMBO2 C++ and CUDA code.)
