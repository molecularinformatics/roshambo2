import h5py
from rdkit import Chem
from rdkit.Chem import rdDepictor

with h5py.File('dataset_split_0.h5', 'r') as f:
    group = f['molecules_0']
    smiles_data = group['smiles'][:]
    names_data = group['names'][:]

smiles_list = [s.decode('utf-8') for s in smiles_data]
names_list = [n.decode('utf-8') for n in names_data]

writer = Chem.SDWriter('dataset_split_0.sdf')

for i, smiles in enumerate(smiles_list):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        rdDepictor.Compute2DCoords(mol)
        mol.SetProp("_Name", names_list[i])
        mol.SetProp("SMILES", smiles)
        writer.write(mol)

writer.close()