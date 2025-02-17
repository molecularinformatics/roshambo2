# Basic usage of roshambo2 program
#
# You need to have:
#  1. query molecule in a 3D sdf file: query.sdf
#  2. dataset molecules in a 3D sdf files: dataset.sdf
# we have provided example files in this folder
# 
# The steps are explained in the comments


# Import the Roshambo2 class from the roshambo2 library
from roshambo2 import Roshambo2
 
# Create an instance of the Roshambo2 class using your query and dataset files.
# In this example we will not compute color features.
roshambo2_calculator = Roshambo2("query.sdf", "dataset.sdf", color=False, verbosity=2)

# Calculate the overlap scores with the compute method.
# In this example we choose the 'cpp' backend because the dataset is small.
# We use optim_mode='shape' because we did not assign color features. 
# We set reduce_over_conformers to False so we get an overlap score for each conformer in the dataset
# We choose start mode 1 which balances accuracy and speed.
scores = roshambo2_calculator.compute(backend='cuda', optim_mode='shape', reduce_over_conformers=False, start_mode=1, keep_order=True) 

# the returned scores object is a dictionary of pandas dataframes. The key is the name of the query molecule
# In this case with one query it is a dictionary with one query name and the corresponding dataframe.
for key in scores:
    print(key)
    df = scores[key]
    print(df)

# we can then write the structures with the best fit overlap transformation to SDF files. 
# There will be one file per query
# Optionally you can append the query molecule to the start of the SDF file
roshambo2_calculator.write_best_fit_structures(append_query=True)
