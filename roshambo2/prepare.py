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
import h5py
from copy import copy, deepcopy
from tqdm import tqdm
from roshambo2.classes import Roshambo2Mol, Roshambo2Data, Roshambo2DataReaderh5
# from roshambo2.utils import _write_sdf
from roshambo2.pharmacophore import PharmacophoreGenerator
from rdkit import Chem
from rdkit.Chem import SDMolSupplier, SmilesMolSupplier
import numpy as np
import psutil
import gc
import logging
import multiprocessing
import math

logger = logging.getLogger(__name__)

Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)


def _process_mol_color(mol, color_generator, keep_smiles=True, remove_Hs_before_color_assignment=False):
    assert(color_generator is not None)    
    return(Roshambo2Mol(mol, color_generator=color_generator, keep_smiles=keep_smiles, 
                       remove_Hs_before_color_assignment=remove_Hs_before_color_assignment))
    

def _process_mol(mol, keep_smiles=True):

    return Roshambo2Mol(mol, keep_smiles=keep_smiles)


def _process_chunk(chunk, color, color_generator, keep_smiles, remove_Hs_before_color_assignment, quiet, id):

    N=len(chunk)
    for i in tqdm(range(N), disable=quiet, position=id, desc=f"preparing:cpu:{id}"):
        mol = chunk[i]
        if color_generator is not None:
            chunk[i]= _process_mol_color(
                    mol,
                    color_generator=color_generator,
                    keep_smiles=keep_smiles,
                    remove_Hs_before_color_assignment=remove_Hs_before_color_assignment,
                )
            
        else:
            chunk[i] = _process_mol(mol, keep_smiles=keep_smiles)
    
    return chunk


def _prepare_mols_inplace(mols, color=False, color_generator=None, quiet=False, keep_smiles=True,
    remove_Hs_before_color_assignment=False, n_cpus=None):

    if color is True:
        assert (color_generator is not None)

    N = len(mols)
    
    if n_cpus is None:
        n_cpus = multiprocessing.cpu_count()

    # Split the list into chunks
    chunk_size = math.ceil(N / n_cpus)
    chunks = [mols[i:i + chunk_size] for i in range(0, N, chunk_size)]

    args = [(chunk, color, color_generator, keep_smiles, remove_Hs_before_color_assignment, quiet, id) for id,chunk in enumerate(chunks)]

    with multiprocessing.Pool(n_cpus) as pool:
        results = list( pool.starmap(_process_chunk, args))

    
    mols[:] = [mol for chunk in results for mol in chunk]

    assert(len(mols) == N)
    return mols



def prepare_from_rdkitmols(mols, color=False, color_generator=None, keep_original_coords=False, remove_Hs_before_color_assignment=False, quiet=False, n_cpus=None):
    """Prepares a set of RDKit molecules for Roshambo2.

        Takes a list of RDKit molecules and returns an Roshambo2Data object

        Args:
            mols (List[rdkMol]): List of RDKit molecules.
            color (bool, Optional): If True color features will be assigned using the color_generator.
            color_generator (PharmacophoreGenerator, Optional). If passed this color_generator will be used. If None (the default) the Roshambo2 default color generator will be used. *Note the color_generator here must have the same features as the one used by the Roshambo2Server.*
            keep_original_coords (bool, optional): If True keep the original molecule coordinates. This should only be set to true for query molecules and is only necessary for query molecules. Default is False.
            remove_Hs_before_color_assignment (bool, optional): If True H atoms will be removed before the color feature 
                    assignment is done. If False the H atoms will be used to assign color features and them removed.
                    Note that H atoms are always removed for the shape calculations.
            n_cpus (int, Optional): Number of CPUs to use for molecule preparation and color assignment. Defaulof None will use all CPUs detected by the multiprocessing library.
        Returns:
            Roshambo2Data: An Roshambo2Data object.
    """

    if color and color_generator is None: # use the default
        color_generator = PharmacophoreGenerator()

    if color_generator is not None:
        color = True

    rdkitmols = []

    names_already_seen = set()
    
    total_confs = 0
    for mol in mols:

        # if we have multiple confs we split into N molecules with 1 conf
        n_confs = mol.GetNumConformers()
        total_confs+=n_confs
        
        if not mol.HasProp('_Name'):
            raise ValueError('RDKit molecules must have the "_Name" property set. '
                     'E.g., mol.SetProp("_Name", "name").')
        mol_name = mol.GetProp('_Name')

        if mol_name in names_already_seen:
            raise ValueError(f'RDKit molecules must have unique names. Molecule name {mol_name} is not unique!')
        else:
            names_already_seen.add(mol_name)

        confs=[]
        for x in mol.GetConformers():
            i = x.GetId()
            confi = Chem.Mol(mol, confId=i)
            confs.append(confi)

        for i,conf in enumerate(confs):

            # set each conf to have a different name
            mol_name = conf.GetProp('_Name')
            new_name = mol_name + f'_{i}'
            conf.SetProp('_Name', new_name)

            # store the rdkit molecule in binary format
            rdkitmols.append(conf.ToBinary())

    assert(len(rdkitmols) == total_confs)
    
    acemols = _prepare_mols_inplace(rdkitmols, color=color, color_generator=color_generator, quiet=quiet, 
                                            keep_smiles=True, 
                                            remove_Hs_before_color_assignment=remove_Hs_before_color_assignment, n_cpus=n_cpus)

    dataset = Roshambo2Data(acemols, keep_original_coords=keep_original_coords)
    return dataset



def prepare_from_file(infile, outname=None, color=False, color_generator=None, return_direct=False, quiet=False, 
                      conformers_have_unique_names=False, max_mols_per_group=1000000, keep_original_coords=False, 
                      remove_Hs_before_color_assignment=False, n_cpus=None):
    """Prepare Roshambo2Data from sdf file that contains 3d conformations.

        If outfile is provided it writes the prepared dataset to a h5 file.
        Or if return_direct=True it returns the Roshambo2Data object.

        Args:
            infile (str): Name of input SDF file.
            outname (str, Optional): Name of output Roshambo2 H5 file. If None and return_direct=True then the Roshambo2Data object will be returned and no file written.
            color (bool, Optional): If True color features will be assigned using the color_generator. Default is False.
            color_generator (PharmacophoreGenerator, Optional). If passed this color_generator will be used. If None (the default) the Roshambo2 default color generator will be used.
            return_direct (bool, Optional): If True then the Roshambo2Data object will be returned and no file written. Default is False.
            conformers_have_unique_names (bool, Optional): Flag indicating whether conformers have unique names. If set to True 
                then it will be assumed that conformers of the same molecule are named such that they differ by `_X` 
                where `X` is a integer. E.g. `ABC_0` and `ABC_1` are conformers of the same molecule. If False then the 
                Roshambo2 program will assume conformers have the same names and will assign the suffixes. Defaults to 
                False. 
            max_mols_per_group (int, Optional): The maximum number of configurations to store inside one Roshambo2Data object and in each Roshambo2 H5 group.
                Note that return direct only works when the total number of configs is smaller than this value. To have
                Multiple groups you must write to H5 files with this function. The Default is 1,000,000 (1 Million).
            keep_original_coords (bool, Optional): Passed to Roshambo2Data constructor. Only needed for preparation of query molecules by main Roshambo2 calculator.
            remove_Hs_before_color_assignment (bool, optional): If True H atoms will be removed before the color feature 
                assignment is done. If False the H atoms will be used to assign color features and them removed.
                Note that H atoms are always removed for the shape calculations.
        Returns:
            None | Roshambo2Data: An Roshambo2Data object if return_direct=True else it returns None.
    
    """
    
    # TODO: input validation

    if color_generator is not None: # assume that if a generator is given that color should be assigned
        color=True

    if color is True and color_generator is None: # use the default generator
        color_generator = PharmacophoreGenerator()

    keep_smiles = True 

    logger.info(f'input file = {infile}, out file = {outname}, assigning color features = {color}')
      
    NREADS = max_mols_per_group # number of mols to read before preparing
    
    # TODO: can check memory usage here and force a file write
    # process = psutil.Process()

    #total_mem = int(psutil.virtual_memory().total/(1024**2))
    #MEM_LIMIT = int(total_mem*0.8)
    #if not quiet:
    ##    print(f'Detected memory is {total_mem} MB. Setting limit of {MEM_LIMIT} MB (80%)')
    
    # store RDKit mols here
    mols = []
    dataset = None
    append = False

    # store Roshambo2Mols here
    processed_mols = []
    
    # store the names we have already read
    names_already_read = {}

    Supplier = SDMolSupplier

    with Supplier(fileName=infile, sanitize=True, removeHs=False) as reader:
        counter = 0
        total_N = len(reader)
        logger.info(f'Reading {infile} which contains {total_N} molecules')

        for mol in tqdm(reader, desc=f"Reading {infile}", disable=quiet):
            if mol is not None:

                # we need to be careful about the names here
                mol_name = mol.GetProp('_Name')

                # unless specified in the input parameters we assume that the input SDF files have molecules 
                # where each conformer has the same name.
                if not conformers_have_unique_names:
                    # if we have already read in this molecule we increment the suffix index
                    if mol_name in names_already_read:
                        
                        index = names_already_read[mol_name]+1

                        new_name = mol_name+'_'+str(index)

                        names_already_read[mol_name] = index

                        mol.SetProp('_Name', new_name)

                    # if it is the first time we see it we call it _0
                    else:
                        names_already_read[mol_name] = 0
                        new_name = mol_name+ '_0'
                        mol.SetProp('_Name', new_name)

                else:
                    if mol_name in names_already_read:
                        raise ValueError(f'Non unique names in {infile}')

                    names_already_read[mol_name] = 0

                # store the rdkit molecule in binary format
                mols.append(mol.ToBinary())
                counter+=1

                if len(mols) >= NREADS:

                    # too many mols read in, need to process
                    logger.debug(f'preparing subset of {len(mols)} molecules')        
                    
                    mols = _prepare_mols_inplace(mols, color=color, color_generator=color_generator, quiet=quiet, 
                                            keep_smiles=keep_smiles, 
                                            remove_Hs_before_color_assignment=remove_Hs_before_color_assignment, n_cpus=n_cpus)
                    
                    # the mols list now contains Roshambo2Mols, append to the list of processes mols
                    processed_mols += mols
                    mols.clear() # and clear the storage
                    logger.debug(f'len processed: {len(processed_mols)}')

                # now check if we need to save to file
                if len(processed_mols) >= max_mols_per_group:

                    if return_direct:
                        #TODO, should be able to hand multiple molecule groups here too
                        raise ValueError("Dataset size is too large, you must prepare the dataset in Roshambo2 H5 format.")

                    
                    dataset = Roshambo2Data(processed_mols[:max_mols_per_group], keep_smiles=True, keep_original_coords=keep_original_coords)
                    processed_mols = processed_mols[max_mols_per_group:]
                    
                    
                    dataset.save_to_h5(outname, append=append, with_smiles=True)

                    append = True # the first time we do not append                         

        # end of loop
        mols = _prepare_mols_inplace(mols, color=color, color_generator=color_generator, quiet=quiet, keep_smiles=keep_smiles, remove_Hs_before_color_assignment=remove_Hs_before_color_assignment, n_cpus=n_cpus)
        processed_mols += mols
        dataset = Roshambo2Data(processed_mols, keep_smiles=True, keep_original_coords=keep_original_coords)

        if return_direct:
            return dataset
        else:
            dataset.save_to_h5(outname, append=append, with_smiles=True )



def combine_datasets(inputs, outname, max_mols_per_group=1000000, verbosity=1):
    """Combines multiple Roshambo2 H5 files into one.

        Given a list of H5 files or a single H5 file with multiple groups it will output a single H5 file. It will load
        each group from the file(s) in turn. If the total number of molecules is less than max_mols_per_group it will 
        merge the groups. If it is higher it will not merge the groups. It will only keep at most two groups in memory.
        And will write previous groups to file.

        Args:
            inputs (str | List[str]): Name of an Roshambo2 H5 file or a list of names.
            outname (str): Name of the output Roshambo2 H5 file.
            max_mols_per_group (int, Optional): Limit that is used to decide if to merge Molecule groups. 
                If The currently loaded molecule group it larger than this limit it will be saved to the output file 
                before the next group is loaded and merged. Default is 1,000,000 (1 Million).


    """
    # TODO: give full control of group size


    data_reader = Roshambo2DataReaderh5(inputs)

    merged = None
    append = False
    for data in data_reader.get_data(read_smiles=True):
        if merged is None:
            merged = data
        else:
            logger.info(f'merging datasets of size ({len(merged)} and {len(data)})')
            merged += data

        current_len = len(merged)

        if current_len >= max_mols_per_group:
            merged.save_to_h5(outname, append=append)
            append = True
            merged = None

    if merged is not None:
        merged.save_to_h5(outname, append=append)


 