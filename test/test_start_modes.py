from roshambo2 import Roshambo2
import os
import pytest


@pytest.mark.parametrize('color',[True, False])
@pytest.mark.parametrize('mode',[0, 1, 2])
def test_start_modes(color, mode):
    scores_list = []
    for backend in ['cpp', 'cuda']:
        roshambo2_calculator = Roshambo2("query.sdf", "dataset.sdf", color=color)

        scores_dict = roshambo2_calculator.compute(backend=backend, color_scores=color, start_mode=mode)

        name = [key for key in scores_dict.keys()][0]

        scores = scores_dict[name]

        scores_list.append(scores)

    # make sure top 3 are the same
    names = [ scores['name'][:3].to_list() for scores in scores_list  ]
    print(names)
    assert( names[0] == names[1])





if __name__ == "__main__":
    test_start_modes(False, 0)
    test_start_modes(False, 1)
    test_start_modes(False, 2)
    test_start_modes(True, 0)
    test_start_modes(True, 1)
    test_start_modes(True, 2)