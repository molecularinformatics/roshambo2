# Installation

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
You will then need to install the extra dependencies (these are listed in doc/env.yaml) to your conda environment:

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

## Troubleshooting

If cmake complains about not finding cuda you will need to make sure `nvcc` can be found on your path.
You might need to do something like this:
```bash
export PATH=/usr/local/cuda/bin:${PATH}
export LD_LIBRARY_PATH=/usr/local/cuda/lib64:${LD_LIBRARY_PATH}
```