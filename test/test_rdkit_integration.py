import pytest


from rdkit import Chem
from rdkit.Chem import rdDistGeom, AllChem
from tqdm import tqdm
from roshambo2 import Roshambo2
from rdkit.Chem import SDMolSupplier
from _roshambo2_cpp import test_overlap_single
import numpy as np
from roshambo2.prepare import prepare_from_rdkitmols

def generate_molecule(smiles, Nconfs=1):
    """ turn a smiles into 3D configuration with Nconfs conformers
    """
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        print("failed for ", smiles)
        return None

    else:
        mol = Chem.AddHs(mol)
        conformers = rdDistGeom.EmbedMultipleConfs(mol, numConfs=Nconfs)
        return mol

def test_rdkit():
    

    # 1. Generate 3D molecules from smiles for the dataset.
    #    We will generate 5 conformers for each molecule.abs

    dataset_file = "example.smi"  # Path to input file containing SMILES strings

    with open(dataset_file, 'r') as file:
        smiles_list = file.readlines()

    dataset_mols = []
    for smiles in tqdm(smiles_list[:10]):
        smiles = smiles.rstrip('\n')

        mol = generate_molecule(smiles, Nconfs=5)
        dataset_mols.append(mol)


    # 2. Generate 3D query molecule (single configuration)
    # query molecule (molecule 0 in the example list of SMILES)
    query_smiles='CC(C1=CC=C2C=C(CCN3N=CC4=C3N=C(N)N3N=C(C5=CC=CO5)N=C43)C=CC2=N1)N1CCOCC1 CHEMBL516753'
    query_mol = generate_molecule(query_smiles, Nconfs=1)

    # write the structure
    writer = Chem.SDWriter('rdkit_query.sdf')
    writer.write(query_mol)


    # also write to file to test this
    roshambo2_data = prepare_from_rdkitmols(dataset_mols, color=True)
    
    roshambo2_data.save_to_h5('from_rdkit.h5')

    # color can be True or False. If True the setup will be a bit slower.
    roshambo2_calculator1 = Roshambo2(query_mol, dataset_mols, color=True)

    # compute scores
    scores1 = roshambo2_calculator1.compute(backend='cpp',reduce_over_conformers=True, optim_mode='combination', write_scores=False)
    print(scores1)
   
   
    # color can be True or False. If True the setup will be a bit slower.
    roshambo2_calculator2 = Roshambo2(query_mol, 'from_rdkit.h5', color=True)

    # compute scores
    scores2 = roshambo2_calculator2.compute(backend='cpp',reduce_over_conformers=True, optim_mode='combination', write_scores=False)
    print(scores2)

    # get the best fit molecules aligned to the query molecule(s) sorted by score
    best_confs_aligned = roshambo2_calculator1.get_best_fit_structures()
    
    # this is a python dict where the keys are the name of the query molecule
    qname = list(best_confs_aligned.keys())[0]
    mols = best_confs_aligned[qname]


    # create a sdf file
    writer = Chem.SDWriter('rdkit_fitted_mols.sdf')
    for mol in mols:
        writer.write(mol)


    # load it back in and check volumes are what we expect
    # read in the hits and check
    with SDMolSupplier('rdkit_query.sdf') as supplier:
        assert (len(supplier) == 1)
        for mol in supplier:
            rq = np.array(mol.GetConformer().GetPositions(), dtype=np.float32)
            print(rq.shape)
    overlaps = []
    with SDMolSupplier('rdkit_fitted_mols.sdf') as supplier:
        for mol in supplier:
            r = np.array(mol.GetConformer().GetPositions(), dtype=np.float32)
            print(r.shape)

            mola = np.ones((len(rq), 4))
            molb = np.ones((len(r), 4))

            mola[:,:3] = rq
            molb[:,:3] = r

            overlap = test_overlap_single(mola,molb)
            overlaps.append(overlap)

    
    key = list(scores1.keys())[0]
    ref_vs1 = scores1[key]['overlap_volume'].to_list()
    ref_vs2 = scores2[key]['overlap_volume'].to_list()
    print(ref_vs1)
    print(ref_vs2)
    print(overlaps)

    for i, (ref_v1, ref_v2, test_v) in enumerate(zip(ref_vs1, ref_vs2, overlaps)):
        diff1 = test_v - ref_v1
        diff2 = test_v - ref_v2
        print(i, diff1, diff2)
        assert(np.fabs(diff1)<0.1)
        assert(np.fabs(diff2)<0.1)




if __name__ == "__main__":
    test_rdkit()