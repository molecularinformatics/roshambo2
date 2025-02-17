# Basic usage of roshambo2 program
# The same as basic_run.py but with color used in the optimizer

from roshambo2 import Roshambo2
 
# The color flag is true here so that color features are assigned when the data is loaded and prepared.
roshambo2_calculator = Roshambo2("query.sdf", "dataset.sdf", color=True, verbosity=1, remove_Hs_before_color_assignment=False)

# to optimize with color we set the optim_mode to combination and the combination_param to 0.5
scores = roshambo2_calculator.compute(backend='cuda', optim_mode='combination', 
   combination_param=0.5, reduce_over_conformers=True, tanimoto_threshold=0.1, max_results=1000000, start_mode=1) 

for key in scores:
    print(key)
    df = scores[key]
    print(df)

roshambo2_calculator.write_best_fit_structures()
