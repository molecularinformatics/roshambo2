from roshambo2 import Roshambo2
import os
import pytest

def test_query_name():
    with pytest.raises(ValueError):
        roshambo2_calculator = Roshambo2(["query.sdf", "query.sdf"], "dataset.sdf", color=False)
        scores = roshambo2_calculator.compute(backend='cpp')

if __name__ == "__main__":
    test_query_name()