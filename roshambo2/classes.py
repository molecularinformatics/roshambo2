# MIT License
# 
# Copyright (c) 2025 molecularinformatics  
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.



import time
import os
import sys
import numpy as np
import h5py
import logging

from rdkit import Chem
from rdkit.Geometry import Point3D

from roshambo2.utils import _pad_with_types, _pad_BN, _pad_BN_int, _pad_BND

logger = logging.getLogger(__name__)

class Roshambo2Mol:
    """ Roshambo2 Molecule object


        On creation an Roshambo2 molecule is translated to its COM and PCA frame.
        If a color generator is given then dummy color atoms are created and stored.
        No H atoms are stored.
        If requested the SMILES string is kept, this enables the molecule to be fully rebuilt as an RDKit molecule with a
        specific 3D configuration.
    
        Attributes:
            r (np.array): Coordinates

    """
    def __init__(self, mol, color_generator=None, keep_smiles=True, remove_Hs_before_color_assignment=False):
        """

            Args:
                mol (rdkit.Chem.rdchem.Mol): RDKit molecule or binary version of an RDkit molecule. 
                color_generator (PharmacophoreGenerator): used to assign color dummy atoms
                keep_smiles (bool): flag to keep the smiles string
                remove_Hs_before_color_assignment (bool): If True H atoms will be removed before the color feature 
                    assignment is done. If False the H atoms will be used to assign color features and them removed.
                    Note that H atoms are always removed for the shape calculations.
        """

        if isinstance(mol, bytes): # if rdkit binary molecule
            mol = Chem.Mol(mol)

        # there should be only one conformer in the mol
        assert(mol.GetNumConformers() == 1)


        # before we do anything we store the original coordinates without Hs
        self.original_coords = np.array(Chem.RemoveHs(mol).GetConformer().GetPositions(), dtype=np.float32)
        
        if remove_Hs_before_color_assignment:
            mol = Chem.RemoveHs(mol)    


        # assign the features and get dummy color atoms
        if color_generator is not None:
            self.original_color_coords, self.color_type = color_generator.generate_color_atoms(mol)
            self.has_color=True
        else:
            self.color_coords, self.color_type = None, None
            self.has_color=False
        
        # get the coordinates
        r = np.array(mol.GetConformer().GetPositions(), dtype=np.float32)

        # centre and align
        # important: the COM and PCA are calculated ignoring Hs
        
        centroid = np.mean(self.original_coords, axis=0)
        r = r-centroid # center

        cvm = np.cov(self.original_coords, rowvar=False)
        _,_,v = np.linalg.svd(cvm)
        r = np.dot(r, v.T) # align to principal axis

        # center and align the color dummy atoms
        if self.has_color:
            self.color_coords = self.original_color_coords - centroid
            self.color_coords = np.dot(self.color_coords, v.T) 

        # update the coords of the rdkit mol
        conf = mol.GetConformer()
        for j in range(len(r)):
            x = [float(v) for v in r[j,:]]
            conf.SetAtomPosition(j,Point3D(*x))

        # now we remove the H atoms (unless we already did)
        if not remove_Hs_before_color_assignment:
            mol = Chem.RemoveHs(mol)

        # and now get the coordinates without the Hs
        r = np.array(mol.GetConformer().GetPositions(), dtype=np.float32)

        self.coords = r
        self.n_atoms = r.shape[0] # number of real atoms (not counting Hs)

        if self.has_color:
            # color coords are just dummy atoms so we can join them with the current atoms
            self.coords = np.concatenate((self.coords, self.color_coords))
            self.types = np.concatenate((np.zeros(self.n_atoms), self.color_type))
        else:
            self.types = np.zeros(self.n_atoms)

        self.name = mol.GetProp('_Name')
    
        if keep_smiles: # Save the SMILES and the order so we can rebuild the rdkit molecule later and assign the coordinates in the correct order.
            self.smiles = Chem.MolToSmiles(mol)
            self.smiles_order = list(mol.GetPropsAsDict(True,True)["_smilesAtomOutputOrder"])
        else:
            self.smiles = None
            self.smiles_order = None
        

    def __str__(self):
        attrs = '\n'.join(f"{key}:{value}" for key, value in vars(self).items())
        return f'{self.__class__.__name__}({attrs})'


# TODO: make base DataReader class and these should be subclasses
class Roshambo2DataReaderh5:
    """ Roshambo2 H5 data reader

        This class is a reader for files in Roshambo2 H5 format.

        It does not actually load the data into memory until the get_data() method is called.
    """
    def __init__(self, fnames):
        """

        Args:
            fnames (List): List of filenames

        """
        
        # make a dict where the key is the file name. 
        # the item is a list of the groups
        self.fs_keys = {}

        for fname in fnames:
            with h5py.File(fname,'r') as f:
                molkeys = list(f.keys())
                # print(molkeys)
                logger.debug(f"Scanning {fname} with contains {len(molkeys)} groups")
                for key in molkeys:
                    group = f[key]
                    data = group['r_padded']
                    l = data.shape[0]
                    m = data.shape[1]
                    logger.debug(f"Group <{key}> contains {l} configs with padded length {m}")

            
                self.fs_keys[fname] = molkeys

    # the data is actually read in this generator function which loops over all files and all groups
    def get_data(self, read_smiles=False):
        """Generator that loops over all files and groups, yielding an Roshambo2Data object for each group.

            Args:
                read_smiles (bool): If True, the SMILES of each molecule will be read. 
                                    This is necessary if full RDKit molecules are to be returned.

            Yields:
                Roshambo2Data: The next Roshambo2Data object corresponding to each group.
        """
        for fname, molkeys in self.fs_keys.items(): 
            for key in molkeys:
                yield Roshambo2Data(fname, key=key, keep_smiles=read_smiles)

    def __len__(self):
        # this returns the total number of molecule groups! Not total number of files
        return sum([len(molkeys) for fname,molkeys in self.fs_keys.items()])

class Roshambo2DataReaderSDF:
    """ Roshambo2 SDF data reader

        This class is a reader for SDF files

        This loads and prepares the SDFs on initialization.
    """
    def __init__(self, fnames, color_generator, quiet=False, keep_original_coords=False, **kwargs):
        """

        Args:
            fnames (List): List of SDF file names.
            color_generator: Color generator to be used by roshambo2.prepare.prepare_from_file.
            quiet (bool):
            keep_original_coords (bool):
            kwargs (dict, optional): extra keyword arguments are passed to roshambo2.prepare.prepare_from_file()

        """

        # this init function reads and prepares all the data from the sdf files
        # the datasets are stored in a list
        # import is here to avoid circular dependency
        from roshambo2.prepare import prepare_from_file

        self.datasets = []
        self.fname_dict = {}

        for i,fname in enumerate(fnames):
            dataset = prepare_from_file(fname, color_generator=color_generator, return_direct=True, quiet=quiet, keep_original_coords=keep_original_coords, **kwargs)
            dataset.source_file = fname
            dataset.group_name = 'None'
            self.datasets.append(dataset)
            self.fname_dict[fname] = i
  
       
    def get_data(self, **kwargs):
        """Generator that loops over all files, yielding an Roshambo2Data object for each file

            Yields:
                Roshambo2Data: The next Roshambo2Data object corresponding to each file
        """

        # generator in same format as Roshambo2DataReaderh5.get_data()
        for dataset in self.datasets:
            yield dataset

    # def get_mol_group(self, fname, groupname):

    #     # return dataset corresponding  to the speficied filename. groupname is not used here.
    #     # kept for same format at Roshambo2DataReaderh5.get_mol_group() format

    #     idx = self.fname_dict[fname]
    #     dataset = self.datasets[idx]

    #     return dataset

    def __len__(self):
        return len(self.datasets)


class Roshambo2DataReaderRDKit:
    """ Roshambo2 RDKit data reader

        This class is a reader for a list of RDKit molecules

        This prepares and creates an Roshambo2Data object on initialization.
    """
    def __init__(self, rdkit_mols, color_generator, keep_original_coords=False, quiet=False, **kwargs):
        """

        Args:
            fnames (List): List of RDKit molecules
            color_generator: Color generator to be used by roshambo2.prepare.prepare_from_rdkitmols.
            keep_original_coords (bool):
            kwargs (dict, optional): extra keyword arguments are passed to roshambo2.prepare.prepare_from_rdkitmols()


        """
        #TODO: List of lists of mols, or option to split list into chunks.
        #TODO: can improve the way conformers are dealt with, should only need to get color features once per mol, not once per conf.
        #import is here to avoid circular dependency
        from roshambo2.prepare import prepare_from_rdkitmols

        dataset = prepare_from_rdkitmols(rdkit_mols, color_generator=color_generator, keep_original_coords=keep_original_coords, quiet=quiet,**kwargs)
        self.datasets = [dataset]
    
    def get_data(self, **kwargs):
        """Generator that returns the Roshambo2Data

            Note that currently there will be only 1 roshambo2 data object.
            This method is written to be in the same format as the other Roshambo2DataReaders.

            Yields:
                Roshambo2Data: The Roshambo2Data object
        """
        # generator in same format as Roshambo2DataReaderh5.get_data()
        for dataset in self.datasets:
            yield dataset
    
    def __len__(self):
        return len(self.datasets)
    





class Roshambo2Data:
    """Roshambo2Data object.

    Data format for storing many configurations in a way that optimizes memory access and disk read/write speed for 
    large datasets.

    - The names of the configurations are stored in a list of strings.
    - The coordinates are padded and stored in a numpy array of shape [N,L,4] where N is the number of confgurations, L
      is the length of the largest configuration [number of real atoms + number of dummy color atoms]. The 3rd dimension
      of size four is for the three x,y,z coordinates and forth entry is 0.0 for padded atoms, and 1.0 for real or dummy color
      atoms. In future it could be used to store the guassian width of atoms. But currently the "all carbon radii" approximation 
      is hardcoded. The arrays are always created such that padded atoms appear at the end.
    - If requested the SMILES strings of each configuration are stored in a list of length N. Additionally, the indexing orders that map the 
    - SMILES string to the coordinates are stored in a list of length N.
    - There are several book keeping arrays that record the number of non-padded atoms for each config, the number of real atoms for each config, and the type of each atom.
      The type is stored as an integer. Type 0 refers to a real atom. Types 1 and larger correspond to the color atom type.
      Note that there is no explicit type for a padded atom, this is because the book keeping arrays and the [:,:,4] 
      entries of the coordinates array encode this information.


    Attributes:
        f_names (List[str]): List of the configuration names. It is assumed by parts the program that different conformers of the same molecule are named as 'MOLNAME_X' where 'X' is an integer.
        f_x (np.array, dtype=float32): Numpy array of padded coordinates. Shape [N,L,4] where N is the number of configurations, L is the length of the largest configuration.
        f_n (np.array, dtype=int32): Numpy array that records the number of non-padded atoms for each config. Shape [N].
        f_n_real (np.array, dtype=int32): Numpy array that records the number of real atoms for each config. Shape [N]. 
        color (bool): True if this dataset is for molecules that have had color features assigned. False otherwise.
        f_smiles (List[str], Optional): List of smiles of each config. Will be None if keep_smiles=False was set during creation
        f_smiles_order (np.array, Optional): Numpy array of shape [N,L] of the smiles string order. Will be None if smiles is None.

    """
    #TODO: rename the variables to be more sensible and match the h5 variables.
    def __init__(self, input, key=None, keep_smiles=True, keep_original_coords=False):
        """ 

        Args:
            input (str | list): Can be a list of RDKit molecules or the filename of an Roshambo2 h5 file.
            key (str, optional): If input is an Roshambo2 h5 file this key value must be the name of one of the groups in the Roshambo2 h5 file.
            keep_smiles (bool, optional): If True the smiles are kept. Default is True.
            keep_original_coords (bool, optional): If True keep the original molecule coordinates. This should only be set to true for query molecules and is only necessary for query molecules. Default is False.

        """
        
        self.f_names = None  
        self.f_x = None      
        self.f_n = None      
        self.f_n_real = None
        self.f_types = None
        self.color = None
        self.f_smiles = None
        self.f_smiles_order = None

        # only used for query molecules
        self.original_coords_list = None
        self.original_color_coords_list = None

        
        if input is None:
            # create empty class
            pass
        elif isinstance(input, str): # load from roshambo2 h5
            self._create_from_file_group(input, key, load_smiles=keep_smiles)
        elif isinstance(input, list): # create from list of Roshambo2Mols
            self._create_from_mols(input, keep_smiles=keep_smiles, keep_original_coords=keep_original_coords)
        else:
            raise ValueError('Roshambo2Data must be created from an Roshambo2H5 file or a list of Roshambo2Mols')


    
    def _create_from_mols(self, mols, keep_smiles=True, keep_original_coords=False):

        if not isinstance(mols[0], Roshambo2Mol):
            raise ValueError('molecules must be Roshambo2Mol type')
        
        if mols[0].has_color:
            self.color=True
        else:
            self.color=False

        f_coords = [mol.coords.astype(np.float32)  for mol in mols]
        f_names  = [mol.name    for mol in mols]
        f_n_real = [mol.n_atoms for mol in mols]
        f_types  = [mol.types   for mol in mols]

        f_padded, f_n, f_types = _pad_with_types(f_coords, f_types) # outputs of this are np.arrays

        self.f_names  = f_names
        self.f_x      = f_padded.astype(np.float32)  # coordinates, padded to same length, shape [B, N, 4], [:,:,4] are 0 if it is padded entry
        self.f_n      = f_n       # number of non padded atoms, shape [B]
        self.f_n_real = np.array(f_n_real, dtype=int)  # number of non dummy atoms, shape [B]
        self.f_types  = f_types   # types of each atom (used for color dummy atoms, normal atoms are type 0), shape [B, N]


        if keep_smiles:
            self.f_smiles = [mol.smiles for mol in mols]

            # here we need to pad the ordering indexes
            B = self.f_x.shape[0]
            N = self.f_x.shape[1]

            self.f_smiles_order = np.zeros((B,N), dtype=np.int32)
            for i,mol in enumerate(mols):
                L = len(mol.smiles_order)
                self.f_smiles_order[i,:L] = np.array(mol.smiles_order, dtype=np.int32)
        
        else:
            self.f_smiles = None
            self.f_smiles_order = None

        if keep_original_coords:
            self.original_coords_list=[mol.original_coords for mol in mols]
            if self.color:
                self.original_color_coords_list=[mol.original_color_coords for mol in mols]
            



    def _create_from_file_group(self, fname, molkey, load_smiles=True):

        t1 = time.perf_counter()

        logger.debug(f"reading data from {fname}")
        
        with h5py.File(fname,'r') as f:

            group_name  = str(molkey)
            f_mols = f[molkey]
            self.f_names = f_mols['names'].asstr()[:].tolist()
            self.f_x = f_mols['r_padded'][:].astype(np.float32)
            # logger.debug(f'f_x dtype = {f_x.dtype}')
            self.f_n = f_mols['r_n'][:]
            self.f_types = f_mols['t_padded'][:]
            self.f_n_real = f_mols['r_n_real'][:]

            if 'color' in f_mols.attrs:
                self.color = True if f_mols.attrs['color'] == 'True' else False
            else:
                self.color = False

            if load_smiles:
                self.f_smiles = f_mols['smiles'].asstr()[:].tolist()
                self.f_smiles_order = f_mols['smiles_order'][:]

            else:
                self.f_smiles = None
                self.f_smiles_order = None


        t2 = time.perf_counter()
        logger.debug(f'Read group {molkey} which contains {self.f_x.shape[0]} molecules from {fname} in {t2-t1} seconds')
        



    def save_to_h5(self, fname, append=False, with_smiles=True):
        """Save the Roshambo2Data to an Roshambo2 H5 file.

        This will save the data to an Roshambo2 H5 file. In this format disk 
        read and write can be fast.

        The format is as follows:
        A group is created, if append is False or if it is the first group to be created in the file it will be called 
        molecules_0. 
        If append is true then the suffix will be incremented for each Roshambo2Data that gets written to the file and each
        group in the h5 file will correspond to an Roshambo2Data object. 

        Within the group the attributes of this class are saved as h5py datasets::

            with h5py.File(fname, 'w') as f:
                
                name = 'molecules_0'
                mol_group = f.create_group(name)
                mol_group.create_dataset('names', data = self.f_names, dtype=h5py.string_dtype(encoding='utf-8'))
                mol_group.create_dataset('r_padded', data = self.f_x, dtype=np.float32)
                mol_group.create_dataset('r_n', data = self.f_n)
                mol_group.create_dataset('t_padded', data = self.f_types)
                mol_group.create_dataset('r_n_real', data = self.f_n_real)
                if self.color:
                    mol_group.attrs['color'] = 'True'
                if with_smiles:
                    mol_group.create_dataset('smiles', data = self.f_smiles, dtype=h5py.string_dtype(encoding='utf-8'))
                    mol_group.create_dataset('smiles_order', data=self.f_smiles_order)


        Args:
            fname (str): The name of the Roshambo2H5 file.
            append (bool, optional): If True it will append this Roshambo2Data to an existing H5 file in a new group. Default is False.as_integer_ratio
            with_smiles (bool, optional): If True will save the smiles and smiles_order data. Default is True.
        """

        t1 = time.perf_counter()

        if append:
            with h5py.File(fname, 'a') as f:

                # these will not be in order created
                groupnames = [x for x in f.keys()]
                if len(groupnames) == 0:
                    Warning(f'Appending to empty h5 file {fname} check this is correct.')
                                
                # get the suffixes and find max value, we assume they are in format 'molecules_N' where N is an integer
                suffixes = [int(name.split('_')[1]) for name in groupnames]

                largest = max(suffixes, default=0)
                new_idx = largest+1
                name = f'molecules_{new_idx}'

                logger.info(f'creating {name} in {fname}')
                mol_group = f.create_group(name)
                mol_group.create_dataset('names', data = self.f_names, dtype=h5py.string_dtype(encoding='utf-8'))
                mol_group.create_dataset('r_padded', data = self.f_x, dtype=np.float32)
                mol_group.create_dataset('r_n', data = self.f_n)
                mol_group.create_dataset('t_padded', data = self.f_types)
                mol_group.create_dataset('r_n_real', data = self.f_n_real)
                if self.color:
                    mol_group.attrs['color'] = 'True'
                if with_smiles:
                    mol_group.create_dataset('smiles', data = self.f_smiles, dtype=h5py.string_dtype(encoding='utf-8'))
                    mol_group.create_dataset('smiles_order', data=self.f_smiles_order)



        else:
            with h5py.File(fname, 'w') as f:
                
                name = 'molecules_0'
                print(f'creating {name} in {fname}')
                mol_group = f.create_group(name)
                mol_group.create_dataset('names', data = self.f_names, dtype=h5py.string_dtype(encoding='utf-8'))
                mol_group.create_dataset('r_padded', data = self.f_x, dtype=np.float32)
                mol_group.create_dataset('r_n', data = self.f_n)
                mol_group.create_dataset('t_padded', data = self.f_types)
                mol_group.create_dataset('r_n_real', data = self.f_n_real)
                if self.color:
                    mol_group.attrs['color'] = 'True'
                if with_smiles:
                    mol_group.create_dataset('smiles', data = self.f_smiles, dtype=h5py.string_dtype(encoding='utf-8'))
                    mol_group.create_dataset('smiles_order', data=self.f_smiles_order)


        t2 = time.perf_counter()
        file_size = os.path.getsize(fname)/(1024*1024)
        n_mols = self.f_x.shape[0]
        logger.info(f'written {n_mols} molecules to group <{name}> in file <{fname}> in {t2-t1:.2f}s')


    def tofloat32(self):
        self.f_x = self.f_x.astype(np.float32)


    def __add__(self, other):
        """add method for Roshambo2Data objects.

            This will add two Roshambo2Data objects and pad arrays to the max length.
            e.g. inputs A and B with A.f_x.shape = (NA,LA,4) and B.f_x.shape=(NB,LB,4) which are added as C=A+B
            will result in C.f_x.shape=(NA+NB, max(LA,LB), 4). 

            Note due to the padding that must be done this method can be slow.

        """
        
        assert(self.color == other.color)

        # make a new empty one
        new = Roshambo2Data(None)
        
        M1 = len(self.f_names)
        M2 = len(other.f_names)

        # need to work out the max padding shape
        N1 = self.f_x.shape[1]
        N2 = other.f_x.shape[1]
        N = max(N1, N2)

        # re-pad and add 
        new.f_x = np.concatenate((_pad_BND(self.f_x,N), _pad_BND(other.f_x,N)), axis=0)
        new.f_types = np.concatenate((_pad_BN(self.f_types,N), _pad_BN(other.f_types,N)), axis=0)
        new.f_names = self.f_names + other.f_names # python list type
        
        new.f_n = np.concatenate((self.f_n, other.f_n))
        new.f_n_real = np.concatenate((self.f_n_real, other.f_n_real))

        if self.f_smiles is not None and other.f_smiles is not None:
            new.f_smiles = self.f_smiles + other.f_smiles # python list type

            # repad and add
            NS = max(self.f_smiles_order.shape[1], other.f_smiles_order.shape[1])
            new.f_smiles_order = np.concatenate((_pad_BN_int(self.f_smiles_order, NS), _pad_BN_int(other.f_smiles_order, NS)), axis=0) 

            assert((len(new.f_smiles)) == M1+M2)
            assert((len(new.f_smiles_order)) == M1+M2)
        
        else:
            new.f_smiles = None
            new.f_smiles_order = None
            

        assert(len(new.f_types) == M1+M2)
        assert(len(new.f_names) == M1+M2)

        new.color = self.color

        return new
    
    def __len__(self):
        """Return number of configurations."""
        return len(self.f_names)

    def get_memsize(self):
        """Return total size in bytes."""
        total_size = 0
        # we have lists and np arrays
        for key, value in self.__dict__.items():
            if value is not None:
                if isinstance(value, np.ndarray):
                    size = int(value.nbytes)
                elif isinstance(value, list):
                    size = int(np.sum([sys.getsizeof(_) for _ in value]))
                else:
                    size = int(sys.getsizeof(value))
            else:
                size = 0

            
            total_size += size
            logger.debug(f'Attribute name: {key}, size: {size/(1024*1024)} MB')

        return total_size
        
