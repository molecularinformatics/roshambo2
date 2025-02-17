from roshambo2 import Roshambo2
import os
import pytest


def test_frontend_multi_gpu_list():
    roshambo2_calculator = Roshambo2("query.sdf", ["dataset_split_0.sdf","dataset_split_1.sdf","dataset_split_2.sdf"], color=False)

    ref_scores = roshambo2_calculator.compute(backend='cuda', n_gpus=1)

    print(ref_scores)

    roshambo2_calculator = Roshambo2("query.sdf", ["dataset_split_0.sdf","dataset_split_1.sdf","dataset_split_2.sdf"], color=False, verbosity=2)

    test_scores1 = roshambo2_calculator.compute(backend='cuda', n_gpus=4)

    print(test_scores1)

    name = [key for key in ref_scores][0]

    # make sure top N are the same
    N=10
    names1 = ref_scores[name]['name'][:10].to_list()
    test_names2 = test_scores1[name]['name'][:10].to_list()

    print(names1, test_names2)

    assert(names1 == test_names2)

def test_frontend_multi_gpu_single():
    roshambo2_calculator = Roshambo2("query.sdf", ["dataset.sdf"], color=False)

    ref_scores = roshambo2_calculator.compute(backend='cuda', n_gpus=1)

    roshambo2_calculator = Roshambo2("query.sdf", ["dataset.sdf"], color=False, verbosity=2)

    test_scores1 = roshambo2_calculator.compute(backend='cuda', n_gpus=4)


    name = [key for key in ref_scores][0]

    # make sure top N are the same
    N=10
    names1 = ref_scores[name]['name'][:10].to_list()
    test_names2 = test_scores1[name]['name'][:10].to_list()

    print(names1, test_names2)

    assert(names1 == test_names2)


if __name__ == "__main__":
    test_frontend_multi_gpu_list()
    test_frontend_multi_gpu_single()
 












    
