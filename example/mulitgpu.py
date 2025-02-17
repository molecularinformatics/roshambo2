
import glob
from roshambo2 import Roshambo2
 
dataset="dataset.sdf"
print(dataset)
roshambo2_calculator = Roshambo2("query.sdf", dataset, color=False, verbosity=2, data_mode="in_memory")

n_gpus=4
scores = roshambo2_calculator.compute(backend='cuda', optim_mode='shape', start_mode=1, n_gpus=n_gpus) 


for key in scores:
   print(key)
   df = scores[key]
   print(df)

