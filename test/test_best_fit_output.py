from roshambo2 import Roshambo2
import os
import pytest
from rdkit.Chem import SDMolSupplier
from _roshambo2_cpp import test_overlap_single
import numpy as np
import torch

@pytest.mark.parametrize('backend',['_pytorch','cpp', 'cuda'])
@pytest.mark.parametrize('start_mode',[0,1,2])
@pytest.mark.parametrize('data_file', ['dataset.sdf', 'dataset.h5'])
def test_best_fit_output(backend, start_mode, data_file):
    if backend=='_pytorch' and start_mode !=0:
        return
    
    roshambo2_calculator = Roshambo2("query.sdf", data_file, color=False)

    scores = roshambo2_calculator.compute(backend=backend, optim_mode='shape', reduce_over_conformers=False, start_mode=start_mode)

    print(scores)

    roshambo2_calculator.write_best_fit_structures()

    # read in the hits and check
    with SDMolSupplier('query.sdf') as supplier:
        assert (len(supplier) == 1)
        for mol in supplier:
            rq = np.array(mol.GetConformer().GetPositions(), dtype=np.float32)
            print(rq.shape)
    overlaps = []
    with SDMolSupplier('hits_for_query_CHEMBL221029.sdf') as supplier:
        for mol in supplier:
            r = np.array(mol.GetConformer().GetPositions(), dtype=np.float32)
            print(r.shape)

            mola = np.ones((len(rq), 4))
            molb = np.ones((len(r), 4))

            mola[:,:3] = rq
            molb[:,:3] = r

            overlap = test_overlap_single(mola,molb)
            overlaps.append(overlap)

    
    key = list(scores.keys())[0]
    ref_vs = scores[key]['overlap_volume'].to_list()
    print(ref_vs)
    print(overlaps)

    for i, (ref_v, test_v) in enumerate(zip(ref_vs, overlaps)):
        diff = test_v - ref_v
        print(i, test_v, ref_v, diff)
        assert(np.fabs(diff)<0.1)

if __name__ == "__main__":
    test_best_fit_output('_pytorch', 0, 'dataset.sdf')
    test_best_fit_output('cuda', 0, 'dataset.sdf')
    test_best_fit_output('cuda', 1, 'dataset.sdf')
    test_best_fit_output('cuda', 2, 'dataset.sdf')
    test_best_fit_output('cpp',  0, 'dataset.sdf')
    test_best_fit_output('cpp',  1, 'dataset.sdf')
    test_best_fit_output('cpp',  2, 'dataset.sdf')
    test_best_fit_output('cuda', 0, 'dataset.h5')
    test_best_fit_output('cuda', 1, 'dataset.h5')
    test_best_fit_output('cuda', 2, 'dataset.h5')
    test_best_fit_output('cpp',  0, 'dataset.h5')
    test_best_fit_output('cpp',  1, 'dataset.h5')
    test_best_fit_output('cpp',  2, 'dataset.h5')