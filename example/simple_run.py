
from roshambo2 import Roshambo2
 
roshambo2_calculator = Roshambo2("query.sdf", "dataset.sdf")
scores = roshambo2_calculator.compute() 

print(scores)
