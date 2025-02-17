# example of direct integration with RDKit


from rdkit import Chem
from rdkit.Chem import rdDistGeom, AllChem
from tqdm import tqdm
from roshambo2 import Roshambo2
from rdkit.Chem import Draw

def generate_molecule(smiles, Nconfs=1):
    """ turn a smiles into 3D configuration with Nconfs conformers
    """
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        print("failed for ", smiles)
        return None

    else:
        mol = Chem.AddHs(mol)
        conformers = rdDistGeom.EmbedMultipleConfs(mol, numConfs=Nconfs, useRandomCoords=True)
        return mol

def main():
    

    # 1. Generate 3D molecules from smiles for the dataset.
    #    We will generate 5 conformers for each molecule.abs

    dataset_file = "example.smi"  # Path to input file containing SMILES strings

    with open(dataset_file, 'r') as file:
        smiles_list = file.readlines()

    dataset_mols = []
    for smiles in tqdm(smiles_list):
        smiles = smiles.rstrip('\n')

        mol = generate_molecule(smiles, Nconfs=5)
        dataset_mols.append(mol)


    # 2. Generate 3D query molecule (single configuration)
    # query molecule (molecule 0 in the example list of SMILES)
    query_smiles='CC(C1=CC=C2C=C(CCN3N=CC4=C3N=C(N)N3N=C(C5=CC=CO5)N=C43)C=CC2=N1)N1CCOCC1 CHEMBL516753'
    query_mol = generate_molecule(query_smiles, Nconfs=1)

    # write the structure
    writer = Chem.SDWriter('my_query.sdf')
    writer.write(query_mol)


    # 3. Normal Roshambo2 script but we can directly pass the RDkit molecules for both query and dataset
   
   
    # color can be True or False. If True the setup will be a bit slower.
    roshambo2_calculator = Roshambo2(query_mol, dataset_mols, color=True)

    # compute scores
    scores = roshambo2_calculator.compute(backend='cpp',reduce_over_conformers=False, optim_mode='combination', write_scores=False,)
    print(scores)

    # get the best fit molecules aligned to the query molecule(s) sorted by score
    best_confs_aligned = roshambo2_calculator.get_best_fit_structures()
    
    # this is a python dict where the keys are the name of the query molecule
    qname = list(best_confs_aligned.keys())[0]
    mols = best_confs_aligned[qname]
    
    # create a sdf file
    writer = Chem.SDWriter('fitted_mols.sdf')
    for mol in mols:
        writer.write(mol)

    # Now you can visualize to compare.
    # Note that due to the conformer generation process the top config will not have the exact same config as the query.
    # Even though they are the same molecules


if __name__ == "__main__":
    main()