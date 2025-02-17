import h5py
import numpy as np
from roshambo2.classes import Roshambo2Data
from roshambo2 import Roshambo2


def large_size_make():
    N=10000
    L=100

    # make synthetic dataset
    # use __new__ to skip __init__ which will throw an error when trying to make an empty class
    dataset = Roshambo2Data.__new__(Roshambo2Data)

    names = [f'xxx_{i}' for i in range(N)]
    x = np.random.rand(N,L,4)
    x[:,:,3] = 1.0

    x_types  = np.zeros((N,L))
    x_n      = np.array([L for _ in range(N)])
    x_n_real = np.array([L for _ in range(N)])

    
    dataset.f_names = names
    dataset.f_x = x
    dataset.f_n = x_n
    dataset.f_n_real = x_n_real
    dataset.f_types = x_types
    dataset.color = False
    dataset.f_smiles = ['C'*L for _ in range(N)]
    dataset.f_smiles_order = [np.array(list(range(L))) for _ in range(N)]

    dataset.save_to_h5('test.h5')


def large_size_run():

    calc = Roshambo2("query.sdf", "test.h5")

    scores = calc.compute(backend='cuda', reduce_over_conformers=False, write_scores=False)
    print(scores)

def test_large():
    large_size_make()
    large_size_run()

if __name__ == "__main__":
    test_large()


