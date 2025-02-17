// MIT License
// 
// Copyright (c) 2025 molecularinformatics  
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.


#include <iostream>
#include <stdexcept>
#include <string>

#define CUDA_CHECK_ERROR(err) \
    do { \
        cudaError_t cudaErr = (err); \
        if (cudaErr != cudaSuccess) { \
            throw std::runtime_error(std::string("CUDA error: ") + cudaGetErrorString(cudaErr) + \
                                     " at " + __FILE__ + ":" + std::to_string(__LINE__)); \
        } \
    } while (0)

// threads per block:
#define NTHREADS 256 //TODO: check best value



__device__ __constant__ float EPSILON = 1E-8f;

__device__ __constant__ float PI = 3.14159265358f;

__device__ __constant__ int start_mode_n[3] = {1, 4, 10};


// coordinate arrays are [N,D] (x,y,z, weight/could be used as radius)
__device__ __constant__ int D = 4;



void optimize_overlap_gpu(const float * molA, const int * molA_type, int molA_num_atoms,  int NmolA, 
                          const float * molBs,const int * molB_types, const int * molB_num_atoms, int NmolB, long num_molBs, 
                          const float * rmat, const float * pmat, int N_features, float * scores, bool optim_color,  float lr_q, float lr_t, int nsteps, float mixing_param, int start_mode, int device_id);

