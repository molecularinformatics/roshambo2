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


# This code it used just for the testing framework. 
# It is not used in the main program
# It does not have the functionality to calculate color
# it only works with start mode 0


import torch
from rdkit.Geometry import Point3D
from tqdm import tqdm
import numpy as np


PI = 3.14159265358
KAPPA = 2.41798793102
CARBONRADII2 = 1.7**2
A = KAPPA/CARBONRADII2
CONSTANT = (PI/(2*A))**1.5



def quat_axis_angle(axis, angle):
    qw = torch.cos(0.5 * angle)
    qx = axis[0] * torch.sin(0.5 * angle)
    qy = axis[1] * torch.sin(0.5 * angle)
    qz = axis[2] * torch.sin(0.5 * angle)
    return torch.stack([qw, qx, qy, qz])




def overlap_volume(mol1, mol2, batch_weight_mat):

    a1 = A
    a2 = A

    dist_mat = torch.cdist(mol1, mol2)

    # weight mat masks the mol2 padded molecules
    Kij = torch.exp(-a1*a2*dist_mat**2/(a1+a2))*batch_weight_mat[:,None,:]

    Vij = 8*Kij*CONSTANT

    return torch.sum(Vij, dim=(1,2))


def overlap_volume_single(mol1, mol2):

    a1 = A
    a2 = A
    dist_mat = torch.cdist(mol1, mol2)

    Kij = torch.exp(-a1*a2*dist_mat**2/(a1+a2))

    Vij = 8*Kij*CONSTANT

    return torch.sum(Vij)


def loss_fn(q,t, ref, mat, batch_weight_mat):

    # q is quaternion
    # r is translation

    # make sure it is a unit quaternion
    norm = torch.linalg.norm(q, dim=1)
    q = q/norm[:,None]

    rot_mat = quaternion_to_rotation_matrix(q)

    rot_mat=torch.transpose(rot_mat, 1,2)
    
    rotated_mat = torch.matmul(mat, rot_mat)

    rotated_mat = rotated_mat + t[:,None,:]

    # compute volume
    vol = overlap_volume(ref, rotated_mat, batch_weight_mat)

    return vol


            

def optim(n_optim_steps, q, t, ref_coords, batch_padded_mat, batch_weight_mat):


    # get self volumes first
    vq = overlap_volume_single(ref_coords,ref_coords)
    vqs = torch.tensor([vq for i in range(batch_padded_mat.shape[0])], device=batch_padded_mat.device)
    vrs = []
    for i in range(batch_padded_mat.shape[0]):
        mol2 = batch_padded_mat[i,:,:3]
        wm = batch_weight_mat[i,:]
        idx = wm > 0
        mol2 = mol2[idx,:]
        vr = overlap_volume_single(mol2,mol2)
        vrs.append(vr)

    vrs = torch.tensor(vrs, device=batch_padded_mat.device)

    optimizer = torch.optim.Adagrad([{'params': q, 'lr': 0.1}, {'params': t, 'lr': 0.1}], )

    for ii in range(n_optim_steps):
      
      
      optimizer.zero_grad()
      
      # N overlaps
      overlaps = loss_fn(q,t, ref_coords, batch_padded_mat, batch_weight_mat)
      
      # normalise by self volumes
      losses = -overlaps/(vqs+vrs)

      # take mean so we can call autograd
      loss = torch.mean(losses)
      loss.backward()

      optimizer.step()
        
    with torch.no_grad():
        volumes = loss_fn(q,t, ref_coords, batch_padded_mat, batch_weight_mat)
        

    all_volumes = torch.stack((volumes,vqs,vrs,)).T
    
    return q.detach(), t.detach(), all_volumes.detach()




class PytorchShapeOverlay:
    def __init__(self, query_data, data, start_mode, color_generator=None, mixing=0.0, verbosity=None, n_gpus=None):

        # pytorch backend only works for simple case of shape overlap with single start pose
        assert(data.color is False and query_data.color is False)
        assert(len(query_data.f_names)==1)
        assert(start_mode == 0)
        
        self.data=data
        self.query_data = query_data

        self.device = 'cuda' if torch.cuda.is_available() else 'cpu'

        if self.device == 'cuda':
            print('PyTorch backend is using GPU')
        else:
            print('Warning: PyTorch backend unable to use GPU, it will use CPU which will be slow!')

        self.qs = torch.tensor(self.query_data.f_x, dtype=torch.float32, device=self.device)
        self.fs = torch.tensor(data.f_x, dtype=torch.float32, device=self.device)

        self.scores = np.zeros((len(self.qs), len(self.fs), 20))

    def optimize_overlap(self):

        N = self.fs.shape[0]
    
        # list of transforms
        q = torch.stack([ torch.tensor([1.0,0,0,0], dtype=torch.float32,device=self.device) for x in range(N)])
        t = torch.stack([ torch.tensor([0.0,0.0,0.0], dtype=torch.float32,device=self.device) for x in range(N) ])

        q.requires_grad_()
        t.requires_grad_()

        query = self.qs[0,:,:3]

        fits = self.fs[:,:,:3]
        bwm = self.fs[:,:,3]

        #print(query.shape, fits.shape, bwm.shape)
        q,t, overlap = optim(100, q, t, query, fits, bwm)

        print(overlap.shape)

        norm = torch.linalg.norm(q, dim=1)
        q = q/norm[:,None]

        # rot_mat = quaternion_to_rotation_matrix(q)

        # self.fs[istart:iend] = torch.matmul(self.fs[istart:iend], rot_mat) + t[:,None,:]

  
        # fill scores with the correct info
        # print(self.scores.shape)
        for i,(qi,ti,overlapi) in enumerate(zip(q,t,overlap)):
            # print(qi,ti,overlapi)
            tanimoto = overlapi[0]/(overlapi[1]+overlapi[2]-overlapi[0])
            self.scores[0,i,0] = tanimoto
            self.scores[0,i,1] = tanimoto
            self.scores[0,i,3] = overlapi[0]
            self.scores[0,i,5] = overlapi[1]
            self.scores[0,i,6] = overlapi[2]
            self.scores[0,i,9:13] = qi.cpu().numpy()
            self.scores[0,i,13:16] = ti.cpu().numpy()
            self.scores[0,i,16] = 1. # unit quaternion

        
        return  self.scores
        
        

    # def get_update_coords(self):

    #     N = self.fs.shape[0]

    #     coords = []

    #     for i in range(N):
    #         r = self.fs[i,:self.nfs[i],:].cpu().numpy()
    #         coords.append(r)

    #     return coords
            


@torch.jit.script
def quaternion_to_rotation_matrix(quaternion):
    """
    Convert a quaternion to a 3x3 rotation matrix.

    Parameters:
    - quaternion (torch.Tensor): Input quaternion with batch dimensions (shape: [..., 4])

    Returns:
    - torch.Tensor: 3x3 rotation matrix with batch dimensions (shape: [..., 3, 3])
    """
    w, x, y, z = quaternion.unbind(-1)

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
    rotation_matrix = torch.stack([torch.stack([m00, m01, m02], -1),
                                   torch.stack([m10, m11, m12], -1),
                                   torch.stack([m20, m21, m22], -1)], -2)

    return rotation_matrix
