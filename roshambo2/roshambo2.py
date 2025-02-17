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
import warnings
from typing import List, Dict
import copy
import logging
import psutil
import gc
import tempfile
import shutil

import numpy as np
import pandas
from tqdm import tqdm
import h5py

from rdkit.Chem.rdchem import Mol as rdkitMol
from rdkit import Chem
from rdkit.Geometry import Point3D


from roshambo2.classes import *
from roshambo2.backends import CppShapeOverlay, CudaShapeOverlay
from roshambo2.pharmacophore import PharmacophoreGenerator, BasePharmacophoreGenerator
from roshambo2.utils import quaternion_to_rotation_matrix,  _write_sdf, _smiles_to_3d_noH, _get_prep_transform, _prepare_results, _merge_results, _append_results

# this is needed to make sure properties such as '_Name' are pickled
Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)

logger = logging.getLogger(__name__)


class Roshambo2():
    """ Roshambo2 class

        Main class for roshambo2 calculations
    """
    def __init__(self, queries, datasets, color=False, conformers_have_unique_names=False, data_mode="eager", verbosity=1, 
        color_generator=None, remove_Hs_before_color_assignment=False, n_cpus_prepare=None):
        """ 
        
        Example:

            Create an Roshambo2 calculator::
            
                from roshambo2 import Roshambo2
                roshambo2_calculator = Roshambo2(queries, datasets, color=False, conformers_have_unique_names=False, data_mode="eager", verbosity=1)
            
        
        Args: 
            queries (str or list): The name (or list of names) of the query file(s). Must be SDF with 3D coordinates, 
                Roshambo2 H5 format, or an RDKit molecule.
            datasets (str or list): The name (or list of names) of the query file(s). Must be SDF with 3D coordinates, 
                Roshambo2 H5 format, or can also be a list of RDKit molecules.
            color (bool): Flag indicating whether color information should be considered. If the input is SDF then color 
                features will be generated using RDKit. If the input is Roshambo2 H5 then a check will be done to ensure 
                the data has color assigned. Defaults to False. 
            conformers_have_unique_names (bool): Flag indicating whether conformers have unique names. If set to True 
                then it will be assumed that conformers of the same molecule are named such that they differ by `_X` 
                where `X` is a integer. E.g. `ABC_0` and `ABC_1` are conformers of the same molecule. If False then the 
                Roshambo2 program will assume conformers have the same names and will assign the suffixes. Defaults to 
                False. 
            data_mode (str): "eager" or "in_memory". A setting that controls how datasets are loaded and searched. 
                In "eager" mode each SDF file, or each H5 file group, will be loaded and then searched. After search the
                file/h5 group will be closed and then the next loaded and searched. This mode works well for fast disk,
                small memory, computers. In "in_memory" mode all data is loaded in RAM before searching. Defaults to 
                "eager".
            verbosity (int): 0, 1, or 2. Sets the verbosity for logging output. 0 is quiet, 1 gives some info, 
                2 gives lots of information.
            color_generator: subclass of roshambo2.pharmacophore.PharmacophoreGenerator.
            remove_Hs_before_color_assignment (bool, optional): If True H atoms will be removed before the color feature 
                    assignment is done. If False the H atoms will be used to assign color features and them removed.
                    Note that H atoms are always removed for the shape calculations.
            n_cpus_prepare (int, optional): Number of CPUs to use for multiprocessing the color assignment in the preparation
                    stage. Default is to use all detected cpus.
        """


        # input validation
        # query and datasets are done below
        if not isinstance(color, bool):
            raise ValueError("color argument must be True or False.")
        
        if not isinstance(conformers_have_unique_names, bool):
            raise ValueError("conformers_have_unique_names argument must be True or False.")

        if not  (isinstance(data_mode,str) and (data_mode == "eager" or data_mode == "in_memory")):
            raise ValueError("data_mode argument must be 'eager' or 'in_memory'.")

        if not (isinstance(verbosity, int) and (verbosity in [0,1,2])):
            raise ValueError("verbosity must be an integer with value in [0,1,2].")

        if not ((color_generator is None) or isinstance(color_generator, BasePharmacophoreGenerator)):
            raise ValueError("color_generator must be None or a PharmacophoreGenerator subclass")
        
        if not (isinstance(remove_Hs_before_color_assignment, bool)):
            raise ValueError("remove_Hs_before_color_assignment argument must be True or False.")
        
        if not ((n_cpus_prepare is None) or isinstance(n_cpus_prepare, int)):
            raise ValueError("n_cpus_prepare must be an integer.")


        


        self.verbosity = verbosity

        logging.basicConfig(
            level=[logging.WARNING, logging.INFO, logging.DEBUG][self.verbosity],
            # format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            # datefmt='%Y-%m-%d %H:%M:%S'
        )
        if self.verbosity == 0:
            self.quiet=True
        else:
            self.quiet=False

        # if user wants color but does not provide a generator we use the default one
        if color_generator is None:
            if color == True:
                self.color_generator = PharmacophoreGenerator()
            else: 
                self.color_generator = None
        else:
            self.color_generator = color_generator

        t1 = time.perf_counter()

        if isinstance(datasets, list):
            if all(isinstance(data, rdkitMol) for data in datasets):
                self.data_reader = Roshambo2DataReaderRDKit(datasets, color_generator=self.color_generator, remove_Hs_before_color_assignment=remove_Hs_before_color_assignment, quiet=self.quiet, n_cpus=n_cpus_prepare)
            elif all( dataset.endswith('.sdf') for dataset in datasets):
                self.data_reader = Roshambo2DataReaderSDF(datasets, color_generator=self.color_generator, conformers_have_unique_names=conformers_have_unique_names, remove_Hs_before_color_assignment=remove_Hs_before_color_assignment, quiet=self.quiet,n_cpus=n_cpus_prepare)
            elif all(dataset.endswith('.h5') for dataset in datasets):
                self.data_reader = Roshambo2DataReaderh5(datasets)

            else:
                raise ValueError('dataset input list must be all ".sdf" files, all ".h5" files, or all RDKit molecules')
                
        else:
            if datasets.endswith('.sdf'):
                self.data_reader = Roshambo2DataReaderSDF([datasets], color_generator=self.color_generator, conformers_have_unique_names=conformers_have_unique_names, remove_Hs_before_color_assignment=remove_Hs_before_color_assignment,quiet=self.quiet, n_cpus=n_cpus_prepare)
            elif datasets.endswith('.h5'):
                self.data_reader = Roshambo2DataReaderh5([datasets])
            else:
                raise ValueError('dataset file must be ".sdf" or ".h5"')


        if isinstance(queries, list):
            if all(isinstance(query, rdkitMol) for query in queries):
                self.query_data_reader = Roshambo2DataReaderRDKit(queries, color_generator=self.color_generator, keep_original_coords=True, remove_Hs_before_color_assignment=remove_Hs_before_color_assignment, quiet=self.quiet, n_cpus=1)
            elif all( query.endswith('.sdf') for query in queries):
                self.query_data_reader = Roshambo2DataReaderSDF(queries, color_generator=self.color_generator, conformers_have_unique_names=True, keep_original_coords=True, remove_Hs_before_color_assignment=remove_Hs_before_color_assignment,quiet=self.quiet, n_cpus=1)
            elif all(query.endswith('.h5') for query in queries):
                self.query_data_reader = Roshambo2DataReaderh5(queries)

            else:
                raise ValueError('query input list must be all ".sdf" files or all ".h5" files, or all RDKit molecules')
                
        else:
            if isinstance(queries, rdkitMol):
                self.query_data_reader = Roshambo2DataReaderRDKit([queries], color_generator=self.color_generator, keep_original_coords=True, remove_Hs_before_color_assignment=remove_Hs_before_color_assignment,quiet=self.quiet, n_cpus=1)
            elif queries.endswith('.sdf'):
                self.query_data_reader = Roshambo2DataReaderSDF([queries], color_generator=self.color_generator, conformers_have_unique_names=True, keep_original_coords=True, remove_Hs_before_color_assignment=remove_Hs_before_color_assignment,quiet=self.quiet, n_cpus=1)
            elif queries.endswith('.h5'):
                self.query_data_reader = Roshambo2DataReaderh5([queries])
           
            else:
                raise ValueError('query file must be ".sdf", ".h5", or an RDKit molecule')
        
      

        # merge queries into one dataset
        queries = [query for query in self.query_data_reader.get_data()]
        
        

        query = queries[0]
        for i in range(1,len(queries)):
            query+=queries[i]
        
        self.query_data = query

        # the queries must all have unique names for the ouput dict format to work correctly.
        seen  = set()
        for qname in self.query_data.f_names:
            if qname in seen:
                raise ValueError("query molecules must have unique names")
            else:
                seen.add(qname)

        
        self.dataset_list = []

        if data_mode == "in_memory":
            self.in_memory_mode = True
            self.eager_mode = False


            # preload all data into memory
            N = len(self.data_reader)

            total_size = 0

            logger.info(f"preloading {N} dataframes into memory")

            total_L = 0

            for dataset in tqdm(self.data_reader.get_data(read_smiles=True), total=N):
                total_size +=  dataset.get_memsize() / (1024**2)
                logger.debug(f"size of datasets =  {total_size} MB")
                total_L+=len(dataset)
                
                # # merge in mem
                # if len(self.dataset_list) > 0:
                #     # check if merge
                #     L = len(self.dataset_list[-1])
                #     if L < 1000000: # TODO renable in-memory merge
                #         logger.debug("merging in memory")
                #         self.dataset_list[-1]+=dataset
                #     else:
                #         logger.debug("adding dataset to list of datasets")
                #         self.dataset_list.append(dataset)
                # else:
                self.dataset_list.append(dataset)

                # check merge
                assert(total_L == sum([len(data) for data in self.dataset_list]))

        else:
            self.in_memory_mode = False
            self.eager_mode = True

        t2 = time.perf_counter()
        logger.info(f'Roshambo2 setup completed in {t2-t1}s')

            


    def compute(self, backend = 'cpp', start_mode=1, color_scores=False, optim_mode=None, combination_param=None, reduce_over_conformers=True,
                      write_scores = True, scores_csv_prefix = 'scores_query', max_results=10000, tanimoto_threshold=0.0, n_gpus=1, progress_variable=None,
                      keep_order=False, **kwargs):

        """ Run the shape overlap calculation for the query and dataset


        Examples:
            compute scores with default settings for shape only::

                scores = roshambo2_calculator.compute(self, backend = 'cuda', start_mode=1, color_scores=False, optim_mode='shape')


        Args:
            backend (str, optional): Specifies the backend for computation. Can be `'cuda'` or `'cpp'`. 
                Defaults to `'cpp'`.
            start_mode (int, optional): Specifies the starting mode for computation. Defaults to `1`.

                - `0`: Single start from COM and align to principal axis. (1 start)
                - `1`: Mode 0 plus 180-degree rotations around x, y, z axis. (4 starts)
                - `2`: Mode 0 plus -90, 90, and 180-degree rotations around x, y, z axis. (10 starts)
            
            color_scores (bool, optional): Flag indicating whether to compute color scores. Defaults to `False`.
            optim_mode (str, optional): Specifies the optimization mode. Can be `'shape'`, `'color'` or `'combination'`. 
                Defaults to `'shape'`.
            combination_param (float, optional): Specifies combination parameter to be used if `optim_mode='combination'`. 
                It is a float between `0.0` and `1.0`. Used to compute the `tanimoto_combination`, the objective function 
                to be optimized. `tanimoto_combination = (1-combination_param)*tanimoto_shape + combination_param*tanimoto_color`.
            reduce_over_conformers (bool, optional): Flag indicating whether to only keep the conformer with the highest score. 
                Defaults to `True`.
            write_scores (bool, optional): Flag indicating whether to write the scores to a CSV file. Defaults to `True`.
            scores_csv_prefix (str, optional): String used as the prefix of the scores CSV file. Files with scores for each query 
                are saved with name `<scores_csv_prefix>_<query_name>.csv`. Default is `'scores_query'`.
            max_results (int): Maximum number of results to keep. For large datasets, storing all results can cause 
                memory issues. Default is `10000`.
            tanimoto_threshold (float): Results will only be reported  if tanimoto_combination scores are above this threshold. Default is 0.0.
                Note that this does not override max_results, i.e. if there are more results with score above tanimoto_threshold but 
                max_results is set to a small number only max_results number will be returned. 
                Set max_results to a very high number to avoid this behavior. 
            n_gpus (int): Number of GPUs to use in `"cuda"` mode. Default is `1`.
            progress_variable (Any, optional): A variable used to track progress.
            keep_order (bool, optional): If True the results will be output in the same order as input. 
                To be used reduce_over_conformers must be false and max_results will be ignored. Default is False
        
        Keyword Args:
            optimizer_settings (dict): Settings passed to the optimizer, default is `{'lr_q':0.1, 'lr_t':0.1, 'steps':100}`.
        
        
        Returns:
                Dict[str, pd.DataFrame]: A dictionary of Pandas DataFrames for each query molecule, containing the sorted overlap scores of the best dataset molecules. Each DataFrame includes overlap volumes and self-overlap volumes to facilitate post-processing of other score types. The columns include:

                - `name`: The molecule name from the input SDF or H5. For multiple conformers, conformers are identified by a postfix `_N`, where `N` is an integer.
                - `tanimoto_combo_legacy`: The combined Tanimoto score of shape and color, computed as `tanimoto_color + tanimoto_shape`. Max value is `2.0`.
                - `tanimoto_combination`: The weighted sum of `tanimoto_shape` and `tanimoto_color`, based on `combination_param`. Max value is `1.0`.
                - `tanimoto_shape`: Shape tanimoto. 
                - `tanimoto_color`: Color tanimoto. 
                - `overlap_volume`: Shape overlap volume.
                - `overlap_color`: Color overlap volume.
                - `self_overlap_volume_fit`: Shape self-overlap of the dataset molecule.
                - `self_overlap_color_query`: Shape self-overlap of the query molecule.
                - `self_overlap_color_fit`: Color self-overlap of the dataset molecule.
                - `self_overlap_color_query`: Color self-overlap of the query molecule.
        

        """

        # input validation
        if not(isinstance(backend, str) and backend in ['cuda', 'cpp', '_pytorch']):
            raise ValueError("backend argument must be 'cpp' or 'cuda'.")
        
        if not(isinstance(start_mode, int) and start_mode in [0,1,2]):
            raise ValueError("start_mode argument must in an integer with value in [0,1,2]")

        if not(isinstance(color_scores, bool)):
            raise ValueError("color_scores argument must be True or False")
        
        if not((optim_mode is None) or (isinstance(optim_mode, str) and optim_mode in ["shape", "color", "combination"])):
            raise ValueError("optim_mode argument must be 'shape', 'color' or 'combination'")

        if not(isinstance(reduce_over_conformers, bool)):
            raise ValueError("reduce_over_conformers argument must be True or False")

        if not(isinstance(write_scores, bool)):
            raise ValueError("write_scores argument must be True or False")

        if not(isinstance(scores_csv_prefix, str)):
            raise ValueError("scores_csv_prefix argument must be a str")

        if not(isinstance(max_results, int)):
            raise ValueError("max_results argument must be an integer")
        
        if not(isinstance(n_gpus, int)):
            raise ValueError("n_gpus argument must be an integer")

        if not(isinstance(keep_order, bool)):
            raise ValueError("keep_order argument must be True or False")


        # make sure the options are compatible with each other

        assert (self.query_data is not None)

        if keep_order and reduce_over_conformers:
            logger.error("If 'keep_order=True' reduce_over_conformers must be False")
            assert(reduce_over_conformers is False)


        t1 = time.perf_counter()

        logger.info('Starting compute')
        
        if progress_variable is not None:
            progress_variable.value=0.0
            self.progress=0.0
        else:
            self.progress=0.0
        
        process = psutil.Process()

        self.color_scores = color_scores
        
        if optim_mode is not None:
            if optim_mode == 'shape':
                pass
            elif optim_mode == 'color' or optim_mode == 'combination':
                self.color_scores = True
            else:
                raise ValueError('optim_mode must be one of ["shape", "color", "combination"]')
        else:
            optim_mode = 'shape'
            
        
        if optim_mode == 'combination':
            if combination_param is None:
                warnings.warn('Using optim_mode combination but combination_param has not been set. Defaulting to 0.5')
                mixing = 0.5
            else:
                mixing = combination_param

        elif optim_mode == 'color':
            mixing = 1.0

        else:
            mixing = 0.0
            

        self.color = self.query_data.color

        if self.color_scores and not self.color:
            raise ValueError("Cannot compute color scores with input that was not prepared with color")

        if self.color:
            self.color_generator =  PharmacophoreGenerator()
        else:
            self.color_generator = None

        self.backend_mode = backend

        self.n_searched=0

        # list for each query to keep the top results
        all_results = [   {} for _ in range(len(self.query_data.f_names)) ]

        # loop over each batch of dataset molecules
        if self.eager_mode: # read from disk and compute straight away
            N = len(self.data_reader)
            self.datasource = self.data_reader.get_data(read_smiles=True)
        elif self.in_memory_mode: # all data is in memory
            N = len(self.dataset_list)
            self.datasource = self.dataset_list


        for i in tqdm(range(0,N), disable=(self.verbosity==0)):
        
            # load the next dataset

            if self.eager_mode: # use generator to load from disk
                data = next(self.datasource)
            elif self.in_memory_mode: # they are in a list
                data = self.datasource[i]
            
            
            logger.debug(f"Computing batch {i+1}/{N}")
            
            total_memuse_a = process.memory_info().rss / (1024**2)
            logger.debug(f'total memuse pre-compute: {total_memuse_a}')


            if self.color:
                assert(data.color and self.query_data.color)

            logger.debug('launching backend compute')
            results = self._compute_method(self.query_data, data, start_mode, mixing, n_gpus=n_gpus, **kwargs)
            logger.debug('completed backend compute.')
            logger.debug('merging results.')
            
            # we have the results for this batch
            # we cannot keep all of them because we will run out of memory
           
            # loop over each query
            for k, query_name in enumerate(self.query_data.f_names):
                
                
                if not keep_order:
                
                    # keep only the top ones
                    # sort them by the score.
                    # score is column 0

                    db_names = data.f_names[:]
                    
                    sort_idx = results[k][:,0].argsort()[::-1]

                    # sort the names
                    sorted_db_names = np.array(db_names)[sort_idx]
                    sorted_results = results[k][sort_idx]

                    # make a entry for this specific query+dataset
                    query_k_batch_results = {'names': sorted_db_names, 'scores': sorted_results, 'index': sort_idx, 'batch': np.full(len(sorted_db_names), i)}

                    # reduce over conformers
                    if reduce_over_conformers:
                        query_k_batch_results = self._reduce_over_conformers(query_k_batch_results)


                    # prepare in the correct format and take top N and apply tanimoto_threshold
                    prepared_results = _prepare_results(query_k_batch_results, data, max_results, tanimoto_threshold, logger)


                    # merge with current

                    if (len(all_results[k])) == 0:
                        # first batch, just add
                        all_results[k] = prepared_results

                    else:

                        all_results[k] =_merge_results(all_results[k], prepared_results, max_results, logger)




                else:

                    # dont sort

                    # score is column 0
                    db_names = data.f_names
                
                    # sort_idx = results[k][:,0].argsort()[::-1]
                    idx = np.arange(len(db_names))

                    # dont sort the names
                    unsorted_db_names = np.array(db_names)
                    unsorted_results = results[k]

                    # make a entry for this specific query+dataset
                    query_k_batch_results = {'names': unsorted_db_names, 'scores': unsorted_results, 'index': idx, 'batch': np.full(len(unsorted_db_names), i)}
                    
                    
                    # prepare in the correct format. Will not apply max_results or tanimoto threshold.
                    prepared_results = _prepare_results(query_k_batch_results, data, len(db_names), tanimoto_threshold, logger)


                    # merge with current
                    if (len(all_results[k])) == 0:
                        # first batch, just add
                        all_results[k] = prepared_results

                    else:

                        all_results[k] =_append_results(all_results[k], prepared_results, logger)
                   
                
            
            total_memuse_b = process.memory_info().rss / (1024**2)
            logger.debug(f'total memuse post-compute: {total_memuse_b}')

            if progress_variable is not None:
                progress_variable.value=(i+1)/N - 0.05
                self.progress=(i+1)/N - 0.05
            else:
                self.progress = (i+1)/N - 0.05

            logger.debug(f'progress: {self.progress*100}%')
        
        self.all_results = all_results


        # convert to pandas for nice output
        logger.debug(f'converting results to pandas')
        self.pandas_output_dict = {}
        for n, qname in enumerate(self.query_data.f_names):
            self.pandas_output_dict[qname] = self._convert_results_to_pandas(self.all_results[n])

        if write_scores:
            for query in self.pandas_output_dict:
                out_csv_name = f'{scores_csv_prefix}_{query}.csv'
                df = self.pandas_output_dict[query]
                logger.info(f'saving scores for query {query} to {out_csv_name}')
                df.to_csv(out_csv_name, index=False, sep='\t')


        total_memuse_end = process.memory_info().rss / (1024**2)
        logger.debug(f'total memuse end: {total_memuse_end}')

        self.progress = 1.0
        logger.debug(f'progress: {self.progress*100}%')
        t2 = time.perf_counter()
        logger.info(f'Roshambo2 compute completed in {np.round(t2-t1,3)} seconds. Searched {self.n_searched} configs.')
        return self.pandas_output_dict


    # def _compute_ref_tversky(self, scores, a=0.95):

    #     O_AB = scores['overlap_volume']
    #     O_AA = scores['self_overlap_volume_query']
    #     O_BB = scores['self_overlap_volume_fit']

    #     Tv = O_AB/(a*O_AA + (1-a)*O_BB)
    #     scores['ref_tversky_shape'] = Tv

    #     if self.color:

    #         O_AB = scores['overlap_color']
    #         O_AA = scores['self_overlap_color_query']
    #         O_BB = scores['self_overlap_color_fit']

    #         Tv_color = O_AB/(a*O_AA + (1-a)*O_BB)
    #         scores['ref_tversky_color'] = Tv_color

    def _compute_method(self, query_data, data, start_mode, mixing, n_gpus=1, **kwargs):


        if self.backend_mode == 'cpp':
            logger.debug('Initializing C++ backend')
            ShapeOverlay = CppShapeOverlay
        elif self.backend_mode == 'cuda':
            logger.debug('Initializing CUDA backend')
            ShapeOverlay = CudaShapeOverlay
        elif self.backend_mode == '_pytorch':
            from roshambo2.backends._pytorch_backend import PytorchShapeOverlay
            ShapeOverlay = PytorchShapeOverlay
        else:
            raise ValueError('backend must be one of [cpp, cuda]')

        #memuse_data = sum([data.get_memsize() / (1024**2) for data in data_list])
        #logger.debug(f'mem use dataset: {memuse_data}')


        logger.info('Initializing backend')

        backend = ShapeOverlay(query_data, data, start_mode, color_generator=self.color_generator, mixing=mixing, verbosity=self.verbosity, n_gpus=n_gpus, **kwargs)

        logger.info('Staring optimization')
        t1_opt = time.perf_counter()
        results = backend.optimize_overlap()

        #memuse_results = sum([results.nbytes / (1024**2) for results in results_list])

        #logger.debug(f'mem use volumes: {memuse_results}')
        
        t2_opt = time.perf_counter()

        n_scaned = len(data.f_x)

        self.n_searched+=n_scaned

        logger.info(f'completed optim of {len(query_data.f_x)} x {n_scaned} molecules in {t2_opt-t1_opt} seconds')

        return results



    def get_best_fit_structures(self, top_n=0,  get_color_pseudomols=False, feature_to_symbol_map=None):
        """Get a list of the best fit molecules in RDKit format

        Args:
            top_n (int, optional): Specifies the number of structures to save in order of score. Default is 0 (all structures will be saved).

        Returns:
            Dict[str, List[rdkitMol] ]: A dictionary of lists of RDKit molecules. There will be one list per query. The keys are the query names.

        """

        logger.debug('processing best fit structures')
        # coordinate info
        all_results=self.all_results

        mol_dict = {}

        if get_color_pseudomols:
            try:
                from moleculekit.molecule import Molecule as moleculekitMol
            except ModuleNotFoundError:
                logger.error("you must install moleculekit to output color pseudomolecules")
                raise ImportError

            color_mol_dict = {}
            
            if feature_to_symbol_map is None:
                feature_to_symbol = {F: f'X{i}' for F,i in self.color_generator.FEATURES_ENUM.items()} 
            else:
                feature_to_symbol = feature_to_symbol_map
            
            index_to_feature = self.color_generator.get_index_to_feature()
        
        else:
            color_mol_dict = None

        # loop over each query. we return a dict of lists of rdkit mols
        for n, query_results in zip(range(len(self.query_data)), all_results):

            mol_list = []
            color_mol_list= []

            qname = self.query_data.f_names[n]

            #prepped_query_sdf_fname = f"{prepped_query_sdf_prefix}_{qname}_prepped.sdf"
            
            # if self.query_data.f_binarymols is not None:
                
            #     # we have an rdkit mol so can easily write the sdf
            #     qmol = self.query_data.f_binarymols[n]
            #     mols = [ Chem.Mol(qmol) ]
            #     _write_sdf( prepped_query_sdf_fname, mols)
            
            # elif self.query_data.f_smiles is not None:
            # we have smiles so first need to make a rdkit mol then update the 3D coordinates

            qsmiles = self.query_data.f_smiles[n]
            qsmiles_order = self.query_data.f_smiles_order[n]


            # select the real atoms
            r = self.query_data.f_x[n,:self.query_data.f_n_real[n],:3]

            mol = _smiles_to_3d_noH(qsmiles, qsmiles_order, r)
            mols  = [ mol ]

            # store the reverse of the COM+PCA transformation we did before
            if self.query_data.original_coords_list is not None:
                orig_r = self.query_data.original_coords_list[n]

                Rq,tq = _get_prep_transform(orig_r,r)
            else:
                # cannot properly align
                Rq = None
                tq = None
                logger.warning("Query molecule data does not contain original coords. \
                Output confs will be aligned to the query in it COM+PCA.")

            
            #_write_sdf( prepped_query_sdf_fname, mols)
            
            # else:
            #     # we do not have enough info to build a proper molecule
            #     raise ValueError('data for query molecule does not contain enough info to create an sdf file')


            # loop over all the results for this query
            top_mols=[]
            scores=[]

            for i in range(len(query_results['names'])):

                # we rebuild from smiles

                smiles = query_results['molinfo'][i][0]
                smiles_order = query_results['molinfo'][i][1]

                transform = query_results['scores'][i, 9:20]
                
                # transform the coordinates
                q = transform[0:4]
                t = transform[4:7]
                qstartmode = transform[7:11]

                # print(q,t,qstartmode)
                M = quaternion_to_rotation_matrix(q)

                r = query_results['orig_coords'][i]

                # transformation by fit rotation and translation
                r = np.dot(r, M.T) + t

                # now we transform by the start mode
                Ms = quaternion_to_rotation_matrix(qstartmode)
                r = np.dot(r, Ms)

                # now we transform to original query
                r=np.dot(r,Rq.T) + tq


                if get_color_pseudomols: # do the same transformation with the color coordinates
                    c_r = query_results['orig_color_coords'][i]
                    c_r = np.dot(c_r, M.T) + t
                    c_r = np.dot(c_r, Ms)
                    c_r =np.dot(c_r,Rq.T) + tq
                    # build a psueodmol with moleculekit
                    color_mol = moleculekitMol().empty(c_r.shape[0])
                    color_mol.coords = c_r[...,None] # moleculekit needs [N,3,FRAME_NUM]
                    
                    color_mol.name = np.array(["CA"] * c_r.shape[0])
                    color_mol.record = np.array(["ATOM"] * c_r.shape[0])
                    color_mol.element = np.array([feature_to_symbol[index_to_feature[t]] for t in query_results['color_type'][i]])
                    color_mol.resid = np.arange(c_r.shape[0])
                    color_mol_list.append(color_mol)


                # build the 3d mol with the coordinates
                mol = _smiles_to_3d_noH(smiles, smiles_order, r)


                score = query_results['scores'][i,0]

                mol.SetProp('score', str(round(score,3)) )
                name = query_results['names'][i]
                
                props={}
                props['SMILES'] = smiles
                props['name'] = name
                props['tanimoto_combination'] = query_results['scores'][i,0]
                props['tanimoto_shape'] = query_results['scores'][i,1]
                props['tanimoto_color'] = query_results['scores'][i,2]
                props['overlap_volume'] = query_results['scores'][i,3]
                props['overlap_color']  = query_results['scores'][i,4]
                props['tanimoto_combo_legacy'] = props['tanimoto_shape']+props['tanimoto_color']
                props['query_name'] = qname
                props['query_SMILES'] = qsmiles


                mol.SetProp('_Name', name)
            
                for prop, value in props.items():
                    mol.SetProp(prop, str(value))

                mol.GetConformer().Set3D(True)

                mol_list.append(mol)

            if top_n==0: # get all
                pass
            else:
                # check it is not out of range
                if top_n > len(mol_list):
                    logger.warning(f"specified top_n of ({top_n}) is greater than number of returned molecules ({len(mol_list)}).")
                else:
                    mol_list=mol_list[:top_n]

                    assert(len(mol_list) == top_n)

                    if get_color_pseudomols:
                        color_mol_list = color_mol_list[:top_n]
                        assert(len(color_mol_list) == top_n)

            mol_dict[qname] = mol_list

            if get_color_pseudomols:
                color_mol_dict[qname] = color_mol_list

        if get_color_pseudomols:
            return mol_dict, color_mol_dict
        else:
            return mol_dict

    

    def write_best_fit_structures(self, top_n=0, hits_sdf_prefix='hits_for_query', append_query=False, write_color_pseudomols=False, feature_to_symbol_map=None):
        """Writes the best fit molecules to an SDF file.

        Args:
            top_n (int, optional): Specifies the number of structures to save in order of score. Default is 0 (all structures will be saved).
            hits_sdf_prefix (str optional): Prefix for hits sdf. SDF file of the top results will be created with name '<hits_sdf_prefix>_<query_name>.sdf'.
            append_query (bool optional): If True the query molecule will be written as the first molecule in the SDF file. Default is False.

        Returns:
            None

        """

        if write_color_pseudomols:
            try:
                from moleculekit.molecule import Molecule as moleculekitMol
            except ModuleNotFoundError:
                logger.error("you must install molecule kit to output color pseudomolecules")
                raise ImportError
            
            if feature_to_symbol_map is None:
                feature_to_symbol = {F: f'X{i}' for F,i in self.color_generator.FEATURES_ENUM.items()} 
            else:
                feature_to_symbol = feature_to_symbol_map

            index_to_feature = self.color_generator.get_index_to_feature()

            mol_dict, color_mol_dict = self.get_best_fit_structures(top_n=top_n, get_color_pseudomols=True, feature_to_symbol_map=feature_to_symbol)

        else:
            mol_dict = self.get_best_fit_structures(top_n=top_n, get_color_pseudomols=False)


        for n,qname in enumerate(mol_dict):

            top_mols = mol_dict[qname]

            if write_color_pseudomols:
                top_color_mols = color_mol_dict[qname]

            if append_query:

                # the mol_dict order should be the same as the order of query_data
                assert(qname == self.query_data.f_names[n])

                qsmiles = self.query_data.f_smiles[n]
                qsmiles_order = self.query_data.f_smiles_order[n]

                # we should have the original coordinates saved
                if self.query_data.original_coords_list is not None:
                    orig_r = self.query_data.original_coords_list[n]

                    qmol = _smiles_to_3d_noH(qsmiles, qsmiles_order, orig_r)

                    qmol.SetProp("Query", "True")
                    qmol.SetProp("_Name", qname)
                    qmol.GetConformer().Set3D(True)

                    top_mols = [qmol]+top_mols

                    if write_color_pseudomols:
                        c_r = self.query_data.original_color_coords_list[n]
                        
                        # build a psueodmol with moleculekit
                        color_mol = moleculekitMol().empty(c_r.shape[0])
                        color_mol.coords = c_r[...,None] # moleculekit needs [N,3,FRAME_NUM]
                        
                        color_mol.name = np.array(["CA"] * c_r.shape[0])
                        color_mol.record = np.array(["ATOM"] * c_r.shape[0])
                        color_mol.element = np.array([feature_to_symbol[index_to_feature[t]] for t in self.query_data.f_types[n][self.query_data.f_n_real[n]:self.query_data.f_n[n]]])
                        color_mol.resid = np.arange(c_r.shape[0])
                        
                        top_color_mols = [color_mol] + top_color_mols


                else:
                    logger.warning("Query molecule data does not contain original coords it will not be written to the output SDF.")

          
            # write SDF with RDkit
            _write_sdf(f'{hits_sdf_prefix}_{qname}.sdf', top_mols)

            if write_color_pseudomols:
                # write color pseudo molecules with Moleculekit
                fname = f'{hits_sdf_prefix}_{qname}_color_features.sdf'
                logger.info(f"writing color features to {fname}")

                with tempfile.TemporaryDirectory() as temp_dir:
                    temp_files = [] # work around b.c write does not seem to have append
                    for i, color_mol in enumerate(tqdm(top_color_mols)):
                        temp_file_path = os.path.join(temp_dir, f"tempfile_{i}.sdf")  # Define .sdf file path
                        temp_files.append(temp_file_path)
                        color_mol.write(temp_file_path)  # Write data to this .sdf file
                    with open(fname,'w') as outfile:
                        for temp_file in temp_files:
                            with open(temp_file,'r') as infile:
                                for line in infile.readlines():
                                    outfile.write(line)
                            
                            # need a newline between each mol
                            outfile.write("\n")
            
        

    def _reduce_over_conformers(self, results_dict):

        names_already_seen = set()

        N = len(results_dict['names'])

        # make a list of indexes to keep
        mask = np.full((N), True)

        for i in range(N):
            name = results_dict['names'][i]
            base_name = name.rsplit('_',1)[0]
            # print(name, base_name)
            
            if base_name in names_already_seen:
                mask[i] = False
            else:
                names_already_seen.add(base_name)

            #TODO: maybe better to keep full name?
            results_dict['names'][i] = base_name


        for key in results_dict:
            # print(key)
            results_dict[key] = results_dict[key][mask]

        return results_dict




    
    def _convert_results_to_pandas(self, results_dict):
        # convert results dict into roshambo2 1.0 pandas dataframe

        results_dict['scores'] = np.array(results_dict['scores'])

        data = { 'name':   results_dict['names'],
                 'smiles': [molinfo[0] for molinfo in results_dict['molinfo']], 
                 'tanimoto_combo_legacy': results_dict['scores'][:,1] + results_dict['scores'][:,2], # sum of shape and color so out of 2
                 'tanimoto_combination': results_dict['scores'][:,0],
                 'tanimoto_shape': results_dict['scores'][:,1],
                 'tanimoto_color': results_dict['scores'][:,2],
                 'overlap_volume': results_dict['scores'][:,3],
                 'overlap_color': results_dict['scores'][:,4],
                #  'ref_tversky_shape':
                #  'ref_tversky_color':
                 'self_overlap_volume_query': results_dict['scores'][:,5],
                 'self_overlap_volume_fit':   results_dict['scores'][:,6],
                 'self_overlap_color_query':  results_dict['scores'][:,7],
                 'self_overlap_color_fit':     results_dict['scores'][:,8],
        }

        # print(data)

        return pandas.DataFrame(data)




class Roshambo2ServerMode(Roshambo2):
    """Server mode version of Roshambo2

    The server mode is a subclass of the `Roshambo2` class.
    The key differences are that when you create it you only specify the dataset, not the query. It can load the entire dataset into memory. The query can then be passed as search time with the `search()` method.

    NOTE: This class was designed to be used by the flask sever application that can be found at '`scripts/server_app.py`'.

    Please look at that script before you use this class directly yourself.

    """
    def __init__(self, datasets, color=False, conformers_have_unique_names=False, verbosity=1, color_generator=None, remove_Hs_before_color_assignment=False,n_cpus_prepare=None, n_gpus=1):
        """

        Example:

            Create a server mode Roshambo2 calculator::
            
                roshambo2_server_calculator = Roshambo2ServerMode(datasets, color=False, conformers_have_unique_names=False, verbosity=1)
            
        
        Args: 
            datasets (str or list): The name (or list of names) of the query file(s). Must be SDF with 3D coordinates, 
                Roshambo2 H5 format, or can also be a list of RDKit molecules.
            color (bool): Flag indicating whether color information should be considered. If the input is SDF then color 
                features will be generated using RDKit. If the input is Roshambo2 H5 then a check will be done to ensure 
                the data has color assigned. Defaults to False. 
            conformers_have_unique_names (bool): Flag indicating whether conformers have unique names. If set to True 
                then it will be assumed that conformers of the same molecule are named such that they differ by `_X` 
                where `X` is a integer. E.g. `ABC_0` and `ABC_1` are conformers of the same molecule. If False then the 
                Roshambo2 program will assume conformers have the same names and will assign the suffixes. Defaults to 
                False.
            verbosity (int): 0, 1, or 2. Sets the verbosity for logging output. 0 is quiet, 1 gives some info, 
            color_generator: subclass of roshambo2.pharmacophore.PharmacophoreGenerator.
            remove_Hs_before_color_assignment (bool, optional): If True H atoms will be removed before the color feature 
                assignment is done. If False the H atoms will be used to assign color features and them removed.
                Note that H atoms are always removed for the shape calculations.
            n_cpus_prepare (int, optional): Number of CPUs to use for multiprocessing the color assignment in the preparation
                stage. Default is to use all detected cpus.
            n_gpus (int): number of gpus to use for searches
        """
        self.verbosity = verbosity

        assert(self.verbosity in [0,1,2])
        
        self.n_gpus=n_gpus

        logging.basicConfig(
            level=[logging.WARNING, logging.INFO, logging.DEBUG][self.verbosity],
            # format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            # datefmt='%Y-%m-%d %H:%M:%S'
        )

        # if user wants color but does not provide a generator we use the default one
        if color_generator is None:
            if color == True:
                self.color_generator = PharmacophoreGenerator()
            else: 
                self.color_generator = None
        else:
            self.color_generator = color_generator

        t1 = time.perf_counter()

        if isinstance(datasets, list):
            if all(isinstance(data, rdkitMol) for data in datasets):
                self.data_reader = Roshambo2DataReaderRDKit(datasets, color_generator=self.color_generator, remove_Hs_before_color_assignment=remove_Hs_before_color_assignment, n_cpus=n_cpus_prepare)
            elif all( dataset.endswith('.sdf') for dataset in datasets):
                self.data_reader = Roshambo2DataReaderSDF(datasets, color_generator=self.color_generator, conformers_have_unique_names=conformers_have_unique_names, remove_Hs_before_color_assignment=remove_Hs_before_color_assignment, n_cpus=n_cpus_prepare)  
            elif all(dataset.endswith('.h5') for dataset in datasets):
                self.data_reader = Roshambo2DataReaderh5(datasets)
            else:
                raise ValueError('dataset input list must be all ".sdf" files, all ".h5" files, or all RDKit molecules')
        else:
            if datasets.endswith('.sdf'):
                self.data_reader = Roshambo2DataReaderSDF([datasets], color_generator=self.color_generator, conformers_have_unique_names=conformers_have_unique_names, remove_Hs_before_color_assignment=remove_Hs_before_color_assignment, n_cpus=n_cpus_prepare)
            elif datasets.endswith('.h5'):
                self.data_reader = Roshambo2DataReaderh5([datasets])
            else:
                raise ValueError('dataset file must be ".sdf" or ".h5"')


    
        self.query_data = None
        
        self.dataset_list = []

        data_mode = "in_memory"

        if data_mode == "in_memory":
            self.in_memory_mode = True
            self.eager_mode = False


            # preload all data into memory
            N = len(self.data_reader)

            total_size = 0

            logger.info(f"preloading {N} dataframes into memory")

            total_L = 0

            #TODO: check if there is enough mem to read it all in
            for dataset in tqdm(self.data_reader.get_data(read_smiles=True), total=N):
                total_size +=  dataset.get_memsize() / (1024**2)
                logger.debug(f"size of datasets =  {total_size} MB")
                total_L+=len(dataset)
                
                # # merge in mem
                # if len(self.dataset_list) > 0:
                #     # check if merge
                #     L = len(self.dataset_list[-1])
                #     print(L)
                #     if L < 1000000: # TODO re-enable based on GPU size
                #         logger.debug("merging in memory")
                #         self.dataset_list[-1]+=dataset
                #     else:
                #         logger.debug("adding to list")
                #         self.dataset_list.append(dataset)
                #else:
                self.dataset_list.append(dataset)

                # check merge
                assert(total_L == sum([len(data) for data in self.dataset_list]))

        else:
            self.in_memory_mode = False
            self.eager_mode = True

        t2 = time.perf_counter()
        logger.info(f'Roshambo2 setup completed in {t2-t1}s')

    def search(self, query_data, options, progress_variable, get_structures):
        """Will search the loaded datasets with the given query molecule(s).

        Example:
            search the loaded datasets with a query::
                results = roshambo2_server_calculator.search(query_data, options)

        Args:
            query_data (Roshambo2Dataset): The query molecule(s) prepared and converted to an `Roshambo2Dataset`. 
            options (dict): A dictionary of key-value pairs corresponding to the arguments of the `Roshambo2.compute()` method.
            progress_variable: Used by server_app.
            get_structures (bool): If True the best fit molecules will also be returned as rdkit binary molecules.

        Returns:
            Tuple[Dict, Dict]: Two dictionaries, the first is the scores, the same as returned from .compute() method. The second is the best fit molecules as returned from .get_best_fit_structures() with the key difference that they are in RDKit binary format!
        
        """
        
        # prepare the query and then call compute

        self.query_data = query_data

        print(self.query_data)

        # TODO: more checks?
        options['write_scores']=False
        options['progress_variable']=progress_variable
        options['n_gpus']=self.n_gpus  # server sets this

        print(options)

        scores_dict = self.compute(**options)

        if get_structures:
            
            rdkit_mols_dict = self.get_best_fit_structures()

            # convert to rdkit binary
            binary_rdkit_mols_dict = {} 
            
            for qname in rdkit_mols_dict:
                binary_rdkit_mols_dict[qname] = [ mol.ToBinary() for mol in rdkit_mols_dict[qname]]

        else:
            binary_rdkit_mols_dict = None


        return scores_dict, binary_rdkit_mols_dict
    





