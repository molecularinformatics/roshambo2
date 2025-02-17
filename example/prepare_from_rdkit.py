# example of creating a roshambo2 h5 file from a list of rdkit molecules


from rdkit import Chem
from rdkit.Chem import rdDistGeom, AllChem
from tqdm import tqdm
from roshambo2.prepare import prepare_from_rdkitmols
from roshambo2 import Roshambo2


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


    # 2. convert to roshambo2 dataset object
    roshambo2_data = prepare_from_rdkitmols(dataset_mols, color=True)
    
  
    # 3. save to h5
    roshambo2_data.save_to_h5('from_rdkit.h5')


    # 4. run a search using it

    # query
    query_smiles='CC(C1=CC=C2C=C(CCN3N=CC4=C3N=C(N)N3N=C(C5=CC=CO5)N=C43)C=CC2=N1)N1CCOCC1 CHEMBL516753'
    query_mol = generate_molecule(query_smiles, Nconfs=1)


    roshambo2_calculator = Roshambo2(query_mol, 'from_rdkit.h5', color=True)

    # compute scores
    scores = roshambo2_calculator.compute(backend='cpp',reduce_over_conformers=True, optim_mode='combination', write_scores=False)
    print(scores)



if __name__ == "__main__":
    main()