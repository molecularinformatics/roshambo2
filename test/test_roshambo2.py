from roshambo2 import Roshambo2
import os
import pytest


# for shape only and start mode 0 we can run the pytorch version
# all three should be the same
def test_frontend_single_query_start_0_shape():
    scores_list = []
    for backend in ['_pytorch', 'cpp', 'cuda']:
        roshambo2_calculator = Roshambo2("query.sdf", "dataset.sdf", color=False)

        scores_dict = roshambo2_calculator.compute(backend=backend, start_mode=0)

        name = [key for key in scores_dict.keys()][0]

        scores = scores_dict[name]

        scores_list.append(scores)

    # make sure top N are the same
    names = [ scores['name'][:10].to_list() for scores in scores_list  ]
    print(names)
    assert( names[0] == names[1] == names[2])

def test_frontend_single_query_start_0_shape_h5():
    scores_list = []
    for backend in ['_pytorch', 'cpp', 'cuda']:
        roshambo2_calculator = Roshambo2("query.sdf", "dataset.h5", color=False)

        scores_dict = roshambo2_calculator.compute(backend=backend, start_mode=0)

        name = [key for key in scores_dict.keys()][0]

        scores = scores_dict[name]

        scores_list.append(scores)

    # make sure top N are the same
    names = [ scores['name'][:10].to_list() for scores in scores_list  ]
    print(names)
    assert( names[0] == names[1] == names[2])


# test color scores are the same for cpp and cuda
@pytest.mark.parametrize('start_mode',[0,1,2])
def test_frontend_single_color(start_mode):
    scores_list = []
    for backend in ['cpp', 'cuda']:
        roshambo2_calculator = Roshambo2("query.sdf", "dataset.sdf", color=True)

        scores_dict = roshambo2_calculator.compute(backend=backend, start_mode=start_mode)
        print(scores_dict)
        name = [key for key in scores_dict.keys()][0]

        scores = scores_dict[name]

        scores_list.append(scores)

    # make sure top N are the same
    names = [ scores['name'][:10].to_list() for scores in scores_list  ]
    print(names)
    assert( names[0] == names[1])


# test that when there are more than one query the cpp and cuda give the same answers
def test_frontend_multiple_querys():
    for backend in ['cpp', 'cuda']:
        roshambo2_calculator = Roshambo2("querys.sdf", "dataset.sdf", color=False)

        scores_dict = roshambo2_calculator.compute(backend=backend)

        for name in scores_dict:
            scores = scores_dict[name]
            print(name)
            print(scores)
            top_scores = scores['name'].to_list()[:5]

            print(name, top_scores)
            assert(name.split('_')[0] == top_scores[0].split('_')[0])




if __name__ == "__main__":
    test_frontend_single_query_start_0_shape()
    test_frontend_single_query_start_0_shape_h5()
    test_frontend_single_color(0)
    test_frontend_single_color(1)
    test_frontend_single_color(2)
    test_frontend_multiple_querys()












    