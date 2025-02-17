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

import copy
import gc
import numpy as np
from tqdm import tqdm
from scipy.spatial.transform import Rotation 
from rdkit.Chem import MultithreadedSDMolSupplier, SDMolSupplier, SDWriter, rdDepictor, MolFromSmiles
from rdkit.Geometry import Point3D

#TODO: not all of these functions are used

def _write_sdf(fname, mols):
    print(f'writing {fname}')
    with SDWriter(fname) as writer:
        for mol in tqdm(mols):
            writer.write(mol)

def _sdfReader(fname, removeHs=False, num_threads=1):
    mols = []
    print(f'reading {fname} with {num_threads} threads')

    # if num_threads == 1:
    with SDMolSupplier(fname, removeHs=removeHs, sanitize=True) as supplier:
        for mol in tqdm(supplier):
            if mol is not None:
                mols.append(mol)

    # this changes orders of things    
    # with MultithreadedSDMolSupplier(fname, removeHs=removeHs, sanitize=True, numWriterThreads=num_threads) as supplier:
    #     for mol in tqdm(supplier):
    #         if mol is not None:
    #             mols.append(mol.ToBinary())

    return mols

def _pad_with_types(mols, types, round_to_multiple=4):
    """
    turns a list of arrays into one array where the length
    is padded to the max. Adds a new column which has value zero if that
    row has been padded
    """

    max_len = max(a.shape[0] for a in mols)

    rounded_max = int(round_to_multiple*np.ceil(max_len/round_to_multiple))
    assert(rounded_max>=max_len)

    max_len = rounded_max

    padded = np.zeros((len(mols), max_len, 4), dtype=np.float32)
    ns = np.zeros(len(mols), dtype=int)
    padded_types = np.zeros((len(mols), max_len), dtype=int)

    # for i,(mol, type) in enumerate(tqdm(zip(mols, types), total=len(mols))):
    for i,(mol, type) in enumerate(zip(mols, types)):
        l = mol.shape[0]
        padded[i,:l,:3] = mol[:,:]
        padded[i,:l, 3] = 1.0
        ns[i]=l
        padded_types[i,:l] = type[:]

    return padded, ns, padded_types

def _pad_BND(arr, M):
    B,N,D = arr.shape
    delta = M-N
    
    assert (delta>=0)
    if delta == 0:
        return arr

    padding = [(0,0), (0,delta), (0,0)]

    padded_arr = np.pad(arr, padding, mode='constant', constant_values=0.0)
    assert padded_arr.shape == (B,M,D)

    return padded_arr

def _pad_BN(arr,M):
    B,N = arr.shape
    delta = M-N

    assert (delta>=0)
    if delta == 0:
        return arr

    padding = [(0,0), (0,delta)]

    padded_arr = np.pad(arr, padding, mode='constant', constant_values=0.0)
    assert padded_arr.shape == (B,M)

    return padded_arr

def _pad_BN_int(arr,M):
    B,N = arr.shape
    delta = M-N

    assert (delta>=0)
    if delta == 0:
        return arr

    padding = [(0,0), (0,delta)]

    padded_arr = np.pad(arr, padding, mode='constant', constant_values=0)
    assert padded_arr.shape == (B,M)

    return padded_arr




def _generate_start_modes(data, mode=1):
    
    # do start modes on the query molecule

    print('Generating starting modes')

    N = data.f_x.shape[0]

    #print(N)
    
    starts = []
    starts_n = []
    starts_n_real = []
    starts_type   = []

    start_indexes = [[] for _ in range(N)]
    j=0
    for i in tqdm(range(N),disable=(N<100)):
        x = data.f_x[i]
        xs = _create_start_coords(x, mode=mode)
        for r in xs:
            starts.append(r)
            starts_n.append(data.f_n[i])
            starts_n_real.append(data.f_n_real[i])
            starts_type.append(data.f_types[i])
            start_indexes[i].append(j)

            j+=1

    starts = np.array(starts)
    starts_n = np.array(starts_n, dtype=int)
    starts_n_real = np.array(starts_n_real, dtype=int)
    starts_type = np.array(starts_type, dtype=int)
    

    return starts, starts_n, starts_n_real, starts_type, start_indexes


def _get_start_mode_transform(k, mode):
    # TODO: this should be used by _create_start_coords to remove code duplication
    
    if k == 0 and mode == 0:
        R = np.eye(3)

    elif mode == 1:
        Rs = []
        Rs.append(np.eye(3))
        for axis in [[1,0,0], [0,1,0], [0,0,1]]:
            rotmat = Rotation.from_rotvec(180*np.array(axis), degrees=True).as_matrix()
            Rs.append(rotmat)

        R = Rs[k]

    elif mode == 2:
        Rs = []
        Rs.append(np.eye(3))
        for axis in [[1,0,0], [0,1,0], [0,0,1]]:
            for angle in [-90,90,180]:
                rotmat = Rotation.from_rotvec(angle*np.array(axis), degrees=True).as_matrix()
                Rs.append(rotmat)
        
        R = Rs[k]
    
    return R
        

def _create_start_coords(r, mode=1):
    # TODO: move this inside the c++/cuda code, no need to realize all the start modes, just need to best one.
    rs = []
    if mode == 0: # no nothing
        rs.append(r)
    elif mode == 1: # 3 rotations

        rs.append(r)
        for axis in [[1,0,0], [0,1,0], [0,0,1]]:
            rotmat = Rotation.from_rotvec(180*np.array(axis), degrees=True).as_matrix()
            
            r_i = copy.deepcopy(r)
            r_i[:,:3] = np.dot(r_i[:,:3], rotmat.T)

            rs.append(r_i)

    elif mode == 2: # 9 rotations
        rs.append(r)

        for axis in [[1,0,0], [0,1,0], [0,0,1]]:
            for angle in [-90,90,180]:
                rotmat = Rotation.from_rotvec(angle*np.array(axis), degrees=True).as_matrix()

                r_i = copy.deepcopy(r)
                r_i[:,:3] = np.dot(r_i[:,:3], rotmat.T)

                rs.append(r_i)

    return rs



def quaternion_to_rotation_matrix(quaternion):

    (w, x, y, z) = quaternion

    # Compute rotation matrix elements
    m00 = 1 - 2*y**2 - 2*z**2
    m01 = 2*x*y - 2*w*z
    m02 = 2*x*z + 2*w*y

    m10 = 2*x*y + 2*w*z
    m11 = 1 - 2*x**2 - 2*z**2
    m12 = 2*y*z - 2*w*x

    m20 = 2*x*z - 2*w*y
    m21 = 2*y*z + 2*w*x
    m22 = 1 - 2*x**2 - 2*y**2

    # Create the rotation matrix
    rotation_matrix = np.array([[m00, m01, m02],
                                   [m10, m11, m12],
                                   [m20, m21, m22]])

    return rotation_matrix





def rotation_matrix_to_quaternion(rotation_matrix):
    """
    Converts a 3x3 rotation matrix to a unit quaternion.
    """
    # Extract the elements of the rotation matrix
    m00, m01, m02 = rotation_matrix[0]
    m10, m11, m12 = rotation_matrix[1]
    m20, m21, m22 = rotation_matrix[2]

    # Calculate the components of the quaternion
    trace = m00 + m11 + m22

    if trace > 0:
        s = 0.5 / np.sqrt(trace + 1.0)
        w = 0.25 / s
        x = (m21 - m12) * s
        y = (m02 - m20) * s
        z = (m10 - m01) * s
    elif m00 > m11 and m00 > m22:
        s = 2.0 * np.sqrt(1.0 + m00 - m11 - m22)
        w = (m21 - m12) / s
        x = 0.25 * s
        y = (m01 + m10) / s
        z = (m02 + m20) / s
    elif m11 > m22:
        s = 2.0 * np.sqrt(1.0 + m11 - m00 - m22)
        w = (m02 - m20) / s
        x = (m01 + m10) / s
        y = 0.25 * s
        z = (m12 + m21) / s
    else:
        s = 2.0 * np.sqrt(1.0 + m22 - m00 - m11)
        w = (m10 - m01) / s
        x = (m02 + m20) / s
        y = (m12 + m21) / s
        z = 0.25 * s

    return np.array([w, x, y, z])


def _smiles_to_3d_noH(smiles: str, smiles_order: list, coords):

    # convert from smiles to molecule
    mol = MolFromSmiles(smiles)

    if mol is None:
        print("failed to convert {smiles} into a RDKit molecule")
        return None
    else:
        # convert to 2D
        rdDepictor.Compute2DCoords(mol)

        
        # assign the 3D coodinates we have specified
        # need to be careful about the order the smiles is in wrt to the coordinates we have
        # TODO: sdf file still says 2d, it should say 3d.
        conf = mol.GetConformer()
        for i in range(mol.GetNumAtoms()):

            j = smiles_order[i]
            x,y,z = (float(coord) for coord in coords[j,:])

            conf.SetAtomPosition(i, Point3D(x,y,z))

        mol.GetConformer().Set3D(True)

        return mol

def _get_prep_transform(A,B):

    ca = np.mean(A,axis=0)
    cb = np.mean(B,axis=0)

    A=A-ca
    B=B-cb

    rot,_ = Rotation.align_vectors(A,B)
    R=rot.as_matrix()

    t = ca-cb

    return R,t



def _prepare_results(results, data, max_results, tanimoto_threshold, logger):
    # assume already sorted

    logger.debug('preparing results')

    in_names  = results['names']
    in_scores = results['scores']
    in_index  = results['index']
    in_batch  = results['batch']
    NB = len(in_names)

    names = []
    scores = []
    index = []
    batch = []
    molinfo = []
    coords = []
    color_coords = []
    color_type = []

    for i in range(min(max_results,NB)):

        score = in_scores[i][0]
        if score >= tanimoto_threshold:

            names.append( copy.deepcopy(in_names[i]))
            scores.append(copy.deepcopy(in_scores[i]))
            index.append( copy.deepcopy(in_index[i]))
            batch.append( copy.deepcopy(in_batch[i]))
            molinfo.append( ( copy.deepcopy(data.f_smiles[ in_index[i] ][:]), copy.deepcopy(data.f_smiles_order[ in_index[i] ][:] )) )
            coords.append(  copy.deepcopy(data.f_x[in_index[i],:data.f_n_real[in_index[i]],:3][:]) )
            color_coords.append( copy.deepcopy(data.f_x[in_index[i],data.f_n_real[in_index[i]]:data.f_n[in_index[i]],:3][:]) )
            color_type.append( copy.deepcopy(data.f_types[in_index[i],data.f_n_real[in_index[i]]:data.f_n[in_index[i]]][:]) )


    return {'names': names, 'scores': scores, 'index': index, 'batch': batch, 'molinfo': molinfo, 'orig_coords': coords, 'orig_color_coords': color_coords, 'color_type': color_type}



def _merge_results(resultsA, resultsB, max_results, logger):

    names = []
    scores = []
    index = []
    batch = []
    molinfo = []
    coords = []
    color_coords = []
    color_type = []

    NA = len(resultsA['names'])
    NB = len(resultsB['names'])

    logger.debug(f"merging results of size {NA} {NB} into  {max((max_results, NA,NB))}")

    # loop over both lists
    i=0
    j=0
    k=0
    while i < NA and j < NB and k<max_results:
        if resultsA['scores'][i][0] >= resultsB['scores'][j][0]:
            #copy from A
            names.append(    resultsA['names'][i])
            scores.append(   resultsA['scores'][i])
            index.append(    resultsA['index'][i])
            batch.append(    resultsA['batch'][i])
            molinfo.append(  resultsA['molinfo'][i])
            coords.append(   resultsA['orig_coords'][i])
            color_coords.append(resultsA['orig_color_coords'][i])
            color_type.append(resultsA['color_type'][i])
                    

            i+=1
        else:
            # copy from B
            names.append(  resultsB['names'][j])
            scores.append( resultsB['scores'][j])
            index.append(  resultsB['index'][j])
            batch.append(  resultsB['batch'][j])
            molinfo.append(resultsB['molinfo'][j])
            coords.append( resultsB['orig_coords'][j])
            color_coords.append(resultsB['orig_color_coords'][j])
            color_type.append(resultsB['color_type'][j])
            j+=1

        k+=1


     # add any remaining
    while i < NA and k < max_results:
        names.append(resultsA['names'][i])
        scores.append(resultsA['scores'][i])
        index.append(resultsA['index'][i])
        batch.append(resultsA['batch'][i])
        molinfo.append(resultsA['molinfo'][i])
        coords.append(resultsA['orig_coords'][i])
        color_coords.append(resultsA['orig_color_coords'][i])
        color_type.append(resultsA['color_type'][i])

        i+=1
        k+=1

    while j < NB and k < max_results:
        names.append(  resultsB['names'][j])
        scores.append( resultsB['scores'][j])
        index.append(  resultsB['index'][j])
        batch.append(  resultsB['batch'][j])
        molinfo.append(resultsB['molinfo'][j])
        coords.append( resultsB['orig_coords'][j])
        color_coords.append(resultsB['orig_color_coords'][j])
        color_type.append(resultsB['color_type'][j])

        j+=1
        k+=1



    out = {'names': names, 'scores': scores, 'index': index, 'batch': batch, 'molinfo': molinfo, 'orig_coords': coords, 'orig_color_coords': color_coords, 'color_type': color_type}

    
    
    return  out



def _append_results(resultsA, resultsB, logger):

    names = []
    scores = []
    index = []
    batch = []
    molinfo = []
    coords = []
    color_coords = []
    color_type = []

    NA = len(resultsA['names'])
    NB = len(resultsB['names'])

    logger.debug(f"appending results of size {NB} to {NA}")

   
    for i in range(NA):
            #copy from A
            names.append(    resultsA['names'][i])
            scores.append(   resultsA['scores'][i])
            index.append(    resultsA['index'][i])
            batch.append(    resultsA['batch'][i])
            molinfo.append(  resultsA['molinfo'][i])
            coords.append(   resultsA['orig_coords'][i])
            color_coords.append(resultsA['orig_color_coords'][i])
            color_type.append(resultsA['color_type'][i])
                    

    for j in range(NB):
            # copy from B
            names.append(  resultsB['names'][j])
            scores.append( resultsB['scores'][j])
            index.append(  resultsB['index'][j])
            batch.append(  resultsB['batch'][j])
            molinfo.append(resultsB['molinfo'][j])
            coords.append( resultsB['orig_coords'][j])
            color_coords.append(resultsB['orig_color_coords'][j])
            color_type.append(resultsB['color_type'][j])



    out = {'names': names, 'scores': scores, 'index': index, 'batch': batch, 'molinfo': molinfo, 'orig_coords': coords, 'orig_color_coords': color_coords, 'color_type': color_type}

    
    
    return  out

