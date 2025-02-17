from roshambo2 import Roshambo2
import os
import pytest


@pytest.mark.parametrize('color',   [True,False])
@pytest.mark.parametrize('backend', ['cpp', 'cuda'])
@pytest.mark.parametrize('groupsize', [10,100,1000])
def test_prepare(color, backend, groupsize):

    print(color)

    if color:
        os.system(f'prepare_dataset_from_sdf.py --color dataset.sdf test.h5 --max_mols_per_group={groupsize}')
    else:
        os.system(f'prepare_dataset_from_sdf.py  dataset.sdf test.h5 --max_mols_per_group={groupsize}')


    roshambo2_calculator1 = Roshambo2("query.sdf", "dataset.sdf", color=color)

    ref_scores = roshambo2_calculator1.compute(backend=backend, color_scores=color)

    roshambo2_calculator2 = Roshambo2("query.sdf", "test.h5", color=color, verbosity=2)

    test_scores = roshambo2_calculator2.compute(backend=backend, color_scores=color)
    print(test_scores)
    
    name = [key for key in ref_scores][0]

    # make sure top N are the same
    names1 = ref_scores[name]['name'][:3].to_list()
    names2 = test_scores[name]['name'][:3].to_list()

    print(names1, names2)

    assert(names1 == names2)

if __name__ == "__main__":
    test_prepare(True, 'cpp', 10)
    test_prepare(True, 'cuda', 10)
    test_prepare(False, 'cpp', 10)
    test_prepare(False, 'cuda', 10)
    test_prepare(True, 'cpp', 100)
    test_prepare(True, 'cuda', 100)
    test_prepare(False, 'cpp', 100)
    test_prepare(False, 'cuda', 100)
    test_prepare(True, 'cpp', 1000)
    test_prepare(True, 'cuda', 1000)
    test_prepare(False, 'cpp', 1000)
    test_prepare(False, 'cuda', 1000)




