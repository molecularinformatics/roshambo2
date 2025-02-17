#!/bin/bash

# The slowest parts of the workflow are reading molecules from SDF, preparing them and assigning color.
# For large scale datasets is it best to do this in an initial preparation stage.
# We use the provided script "prepare_dataset_from_sdf.py" to create Roshambo2 a HDF5 file.
# The H5 dataset can be read in very quickly 

# 1. prepare. This will take a few minutes.
# Note: Assigning color is a command line flag
# with color:
prepare_dataset_from_sdf.py --color '../test/confs.sdf' processed_dataset.h5

# without color:
#prepare_dataset_from_sdf.py large_dataset.sdf processed_dataset.h5

# you can specify the number of CPUs to use (default is all of them), e.g. 4
#prepare_dataset_from_sdf.py --color '../test/confs.sdf' processed_dataset.h5 -n 4

# 2. compute scores
python large_run.py
