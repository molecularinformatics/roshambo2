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


import numpy as np
from tqdm import tqdm
from _roshambo2_cuda import  optimize_overlap_color

class CudaShapeOverlay:
    """Wrapper around _roshambo2_cuda.optimize_overlap_color.

    """
    def __init__(self, query_data, data, start_mode, color_generator=None, mixing=None, verbosity=1, n_gpus=1, optimizer_settings={'lr_q':0.1, 'lr_t':0.1, 'steps':100}):
        
        self.data = data 
        self.query_data = query_data
        self.start_mode = start_mode
        self.verbosity = verbosity

        if color_generator is not None:
            self.color = True
            self.color_generator = color_generator
            self.mixing = mixing
        else:
            self.color = False
            self.color_generator = None
            self.mixing=0.0

        self.query_data.tofloat32()

        
        self.data.tofloat32() 

    
        # optimizer settings
        self.lr_q  = optimizer_settings['lr_q']
        self.lr_t  = optimizer_settings['lr_t']
        self.steps = optimizer_settings['steps']

        self.scores = np.zeros((len(self.query_data.f_x), len(self.data.f_names), 20), dtype=np.float32)  # scores and transforms

        self.n_gpus = n_gpus

    def optimize_overlap(self):

        if self.color:

            optimize_overlap_color(self.query_data.f_x, self.query_data.f_types, self.query_data.f_n_real, 
                                   self.data.f_x,       self.data.f_types,       self.data.f_n_real, 
                                   self.color_generator.interaction_matrix_r, 
                                   self.color_generator.interaction_matrix_p, 
                                   self.scores, True, self.mixing, self.lr_q, self.lr_t, self.steps, 
                                   self.start_mode, self.n_gpus, self.verbosity)

        else:

            dummy = np.zeros((1,1), dtype=np.float32)
            optimize_overlap_color(self.query_data.f_x, self.query_data.f_types, self.query_data.f_n_real, 
                                   self.data.f_x,       self.data.f_types,       self.data.f_n_real, 
                                    dummy, 
                                    dummy, 
                                    self.scores, False,0.0, self.lr_q, self.lr_t, self.steps, 
                                    self.start_mode, self.n_gpus, self.verbosity)
            

        return self.scores
    
                  
