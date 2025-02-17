from roshambo2 import Roshambo2

# read in the processed dataset. The query can still be a sdf, but could also be a processed H5 file.
roshambo2_calculator = Roshambo2("query.sdf", ["processed_dataset.h5"], color=True)

# compute scores, keep only the top 100
scores = roshambo2_calculator.compute(backend='cuda', optim_mode='combination', max_results=100)

print(scores)
roshambo2_calculator.write_best_fit_structures()