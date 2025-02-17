from roshambo2 import Roshambo2
from roshambo2.prepare import combine_datasets

import os
import pytest


@pytest.mark.parametrize('backend', ['cpp', 'cuda'])
def test_prepare(backend):

    os.system('prepare_dataset_from_sdf.py  dataset_split_0.sdf test_chunk_0.h5')
    os.system('prepare_dataset_from_sdf.py  dataset_split_1.sdf test_chunk_1.h5')
    os.system('prepare_dataset_from_sdf.py  dataset_split_2.sdf test_chunk_2.h5')

    
    combine_datasets(['test_chunk_0.h5', 'test_chunk_1.h5', 'test_chunk_2.h5'], 'test_combined.h5')



    roshambo2_calculator1 = Roshambo2("query.sdf", "dataset.sdf")

    ref_scores = roshambo2_calculator1.compute(backend=backend)

    roshambo2_calculator2 = Roshambo2("query.sdf", ["test_chunk_0.h5", "test_chunk_1.h5", "test_chunk_2.h5"])

    test_scores1 = roshambo2_calculator2.compute(backend=backend)


    roshambo2_calculator3 = Roshambo2("query.sdf", ["test_combined.h5"])

    test_scores2 = roshambo2_calculator3.compute(backend=backend)
    

    
    name = [key for key in ref_scores][0]

    # make sure top N are the same
    names1 = ref_scores[name]['name'][:3].to_list()
    test_names2 = test_scores1[name]['name'][:3].to_list()
    test_names3 = test_scores2[name]['name'][:3].to_list()

    print(names1, test_names2, test_names3)

    assert(names1 == test_names2 == test_names3)

if __name__ == "__main__":
    test_prepare('cpp')
    test_prepare('cuda')


