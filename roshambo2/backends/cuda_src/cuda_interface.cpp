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



#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <iostream>
#include <stdio.h>
#include <cmath>
#include <stdexcept>
#include <pybind11/stl.h>
#include <string>
#include <cassert>
#include <cuda_runtime_api.h>
#include "cuda_functions.cuh"
#include <omp.h>

namespace py = pybind11;

enum loglevel {
  OFF,
  INFO,
  DEBUG
};



/// @brief optimization function. Entry point to c++ code from python (via pybind)
/// @param A 
/// @param AT 
/// @param AN 
/// @param B 
/// @param BT 
/// @param BN 
/// @param RMAT 
/// @param PMAT 
/// @param V 
/// @param optim_color 
/// @param mixing_param 
void optimize_overlap_color(py::array_t<float> A, py::array_t<int> AT, py::array_t<int> AN,
                            py::array_t<float> B, py::array_t<int> BT, py::array_t<int> BN, 
                            py::array_t<float> RMAT, py::array_t<float> PMAT, py::array_t<float> V,
                            bool optim_color, float mixing_param, float lr_q, float lr_t, int nsteps, 
                            int start_mode_method, int requested_gpus, int loglevel){
    
    const int D = 4;

    // here we query the number of gpus
    int deviceCount = 0;
    CUDA_CHECK_ERROR(cudaGetDeviceCount(&deviceCount));

    if (loglevel == INFO || loglevel == DEBUG){
        printf("roshambo2.cuda: Found %d cuda devices.\n", deviceCount);
    }

    if (deviceCount != requested_gpus){
        printf("roshambo2.cuda: Warning: %d GPUs requested but %d GPUs found on machine!. Check your settings.\n", requested_gpus, deviceCount);

        if (deviceCount < requested_gpus){
            printf("roshambo2.cuda: Error: %d GPUs requested but %d GPUs found on machine!. Check your settings.\n", requested_gpus, deviceCount);
            throw std::runtime_error("roshambo2.cuda: Error: "+std::to_string(requested_gpus)+" GPUs requested but "+std::to_string(deviceCount)+" GPUs found on machine!.");
        }
    }
    if (loglevel == INFO || loglevel == DEBUG){
        printf("roshambo2.cuda: Roshambo2 CUDA running on %d GPUS\n", requested_gpus);
    }

    if (loglevel == DEBUG){
        std::cout << "roshambo2.cuda: Optimizer settings = { lr_q:" <<lr_q<<" lr_t:"<<lr_t<<" steps: "<<nsteps<< "}"<<std::endl;
    }


    int start_index[requested_gpus];
    int chunk_size[requested_gpus];

    // indexing math
    long len_configs = B.unchecked<3>().shape(0);
    
    int base_len  = len_configs/requested_gpus;
    int remainder = len_configs%requested_gpus;

    int l=0;
    for( int k=0; k<requested_gpus; ++k){
        start_index[k]=l;
        
        if(k==0){
            chunk_size[k] = base_len+remainder;
        }
        else{
            chunk_size[k] = base_len;
        }
        l+=chunk_size[k];
    }
    
    if (loglevel == DEBUG && requested_gpus > 1){
        printf("roshambo2.cuda: Using chunks:\n");
        for( int k=0; k<requested_gpus; ++k){
        printf("roshambo2.cuda: chunk %d, start: %d, end: %d\n", k, start_index[k], start_index[k]+chunk_size[k]);
        }
    }

    int n_threads = requested_gpus;
 
    // make lists of pointers and values for each device
    float * molA_devices[n_threads];
    int * molA_type_devices[n_threads];
    float * molBs_devices[n_threads];
    int * molBs_type_devices[n_threads];
    int * molBs_num_atoms_devices[n_threads];
    long NmolBs[n_threads];
    long num_molBss[n_threads];
    float * rmat_devices[n_threads];
    float * pmat_devices[n_threads];
    float * scores_devices[n_threads];
    long scores_sizes[n_threads];
    long molA_sizes[n_threads];
    long molA_type_sizes[n_threads];

    // shared host pointers
    auto molAs  = A.unchecked<3>();
    auto molA_type = AT.unchecked<2>();
    auto molA_num_atoms = AN.unchecked<1>();
    auto rmat = RMAT.unchecked<2>();
    auto pmat = PMAT.unchecked<2>();
    const float * ptr_rmat = rmat.data(0,0);
    const float * ptr_pmat = pmat.data(0,0);
    int N_features = rmat.shape(0);
    int n_querys = molAs.shape(0);
    long NmolA = molAs.shape(1); // number of atoms in molA

    

    // Allocate data to all GPUs
    #pragma omp parallel for num_threads(n_threads)
    for(int thread_id=0; thread_id<n_threads; ++thread_id){

        int device_id = thread_id%deviceCount;
        int n = thread_id;
        int list_idx;

        # pragma omp critical
        {
        if (loglevel == DEBUG){
            printf("roshambo2.cuda: on thread %d on device %d\n", thread_id, device_id);
        }
        } // pragma
            
        // convert arrays them to the direct access versions
        auto molBs  = B.unchecked<3>();
        auto molB_type = BT.unchecked<2>();
        auto molB_num_atoms = BN.unchecked<1>(); // number of real atoms (non color, non dummy)
        auto scores = V.mutable_unchecked<3>();


        assert(scores.shape(2) == 20);

        // set the current device
        CUDA_CHECK_ERROR(cudaSetDevice(device_id));
        
        // Get the current CUDA device
        int currentDevice;
        cudaGetDevice(&currentDevice);
        assert(device_id == currentDevice);

        // Get the properties of the current CUDA device
        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp, currentDevice);


        // Print information about the current CUDA device
        if (loglevel == DEBUG){
            std::cout << "Current CUDA Device: " << deviceProp.name << std::endl;
            std::cout << "Compute Capability: " << deviceProp.major << "." << deviceProp.minor << std::endl;
            std::cout << "Total Global Memory: " << deviceProp.totalGlobalMem / (1024 * 1024) << " MB" << std::endl;
        }

        long malloced = 0;
        // first allocate molBs to device then loop over molAs (the query mols)

        // allocate molBs. these arrays correspond to f_x of the roshambo2data class

        // host pointer
        int    start_idx = start_index[n];
        const float *     ptr_molBs = molBs.data(start_idx,0,0);
        int    molB_shape0 = chunk_size[n];
  
        // device pointer
        float * molBs_device;
        long molBs_size = molB_shape0*molBs.shape(1)*molBs.shape(2)*sizeof(float);
        CUDA_CHECK_ERROR(cudaMalloc((void **)&molBs_device, molBs_size));
        CUDA_CHECK_ERROR(cudaMemcpy(molBs_device, ptr_molBs, molBs_size, cudaMemcpyHostToDevice));
        malloced += molBs_size;
        // store device pointer in the array
        molBs_devices[n] = molBs_device;

        // allocate molB_type. these arrays correspond to f_types of the roshambo2data class

        start_idx = start_index[n];
        const int *  ptr_molB_type = molB_type.data(start_idx,0);
        int molB_type_shape0 = chunk_size[n];
    
        int * molB_type_device;
        long molB_type_size = molB_type_shape0*molB_type.shape(1)*sizeof(int);
        CUDA_CHECK_ERROR(cudaMalloc((void **)&molB_type_device, molB_type_size));
        CUDA_CHECK_ERROR(cudaMemcpy(molB_type_device, ptr_molB_type, molB_type_size, cudaMemcpyHostToDevice));
        malloced += molB_type_size;
        molBs_type_devices[n] = molB_type_device;

        // allocate molB_num_atoms. these arrays correspond to f_n_real of the roshambo2data class

        start_idx = start_index[n];
        const int *  ptr_molB_num_atoms = molB_num_atoms.data(start_idx);
        int molB_num_atoms_shape0 = chunk_size[n];
    
        int * molB_num_atoms_device;
        long molB_num_atoms_size =  molB_num_atoms_shape0*sizeof(int);
        CUDA_CHECK_ERROR(cudaMalloc((void **)&molB_num_atoms_device, molB_num_atoms_size));
        CUDA_CHECK_ERROR(cudaMemcpy(molB_num_atoms_device, ptr_molB_num_atoms, molB_num_atoms_size, cudaMemcpyHostToDevice));
        malloced += molB_num_atoms_size;
        molBs_num_atoms_devices[n] = molB_num_atoms_device;

        // allocate molA to devices. Note here we just allocate the memory, we do not do the copy. That is done in the next loop

        float * molA_device;
        long molA_size = molAs.shape(1)*molAs.shape(2)*sizeof(float);
        CUDA_CHECK_ERROR(cudaMalloc((void **)&molA_device, molA_size));
        malloced += molA_size;
        molA_devices[n] = molA_device;
        molA_sizes[n] = molA_size;

        int * molA_type_device;
        long molA_type_size = molA_type.shape(1)*sizeof(int);
        CUDA_CHECK_ERROR(cudaMalloc((void **)&molA_type_device, molA_type_size));
        malloced += molA_type_size;
        molA_type_devices[n] = molA_type_device;
        molA_type_sizes[n] = molA_type_size;

        // allocate and copy the interaction maps
        float * rmat_device;
        long rmat_size = rmat.shape(0)*rmat.shape(1)*sizeof(float);
        CUDA_CHECK_ERROR(cudaMalloc((void **)&rmat_device, rmat_size));
        CUDA_CHECK_ERROR(cudaMemcpy(rmat_device, ptr_rmat, rmat_size, cudaMemcpyHostToDevice));
        malloced += rmat_size;
        rmat_devices[n] = rmat_device;

        float * pmat_device;
        long pmat_size = pmat.shape(0)*pmat.shape(1)*sizeof(float);
        CUDA_CHECK_ERROR(cudaMalloc((void **)&pmat_device, pmat_size));
        CUDA_CHECK_ERROR(cudaMemcpy(pmat_device, ptr_pmat, pmat_size, cudaMemcpyHostToDevice));
        malloced += pmat_size;
        pmat_devices[n] = pmat_device;

        // allocate the scores array
        int scores_shape1 = chunk_size[n];

        float * scores_device;
        long scores_size = scores_shape1*scores.shape(2)*sizeof(float);
        CUDA_CHECK_ERROR(cudaMalloc((void **)&scores_device, scores_size));
        malloced += scores_size;
        scores_devices[n] = scores_device;
        scores_sizes[n] = scores_size;

        // TODO: this naming is confusing 
        long NmolB = molBs.shape(1); // number of atoms in molB (padded to all be the max)
        long num_molBs = chunk_size[n];

        NmolBs[n] = NmolB;
        num_molBss[n] = num_molBs;
        long memory_malloced_in_mb = malloced / (1024 * 1024);

        if (loglevel == DEBUG){
            printf("roshambo2.cuda: Device %d, memory allocated = %d MB\n", currentDevice, memory_malloced_in_mb);
        }
    } // end for

    // loop over the query molecules
    for (py::ssize_t i = 0; i < molAs.shape(0); i++){

        // std::cout << "Optimizing for query molecule " << i << std::endl;
        if (loglevel == DEBUG){
            printf("roshambo2.cuda: query %d / %d\n", i+1, molAs.shape(0));
        }

        // loop over gpus
        // use OMP multithreading, one thread per gpu
        #pragma omp parallel for num_threads(n_threads)    
        for(int thread_id=0; thread_id<n_threads; ++thread_id){
            
            int n = thread_id;
            int device_id = n;
            
            CUDA_CHECK_ERROR(cudaSetDevice(device_id));

            //  host pointer to A
            const float * ptr_molA = molAs.data(i,0,0);


            // copy to device
            CUDA_CHECK_ERROR(cudaMemcpy(molA_devices[n], ptr_molA, molA_sizes[n], cudaMemcpyHostToDevice));

            const int * ptr_molA_type = molA_type.data(i,0);
            CUDA_CHECK_ERROR(cudaMemcpy(molA_type_devices[n], ptr_molA_type, molA_type_sizes[n], cudaMemcpyHostToDevice));

            int molA_num_atoms_i = molA_num_atoms(i);

            int   list_idx = 0; // list should be length 1, we split over GPUs
            int    idx = start_index[n];

            float * scores_host = V.mutable_data(i,start_index[n],0);
            CUDA_CHECK_ERROR(cudaMemcpy(scores_devices[n], scores_host, scores_sizes[n], cudaMemcpyHostToDevice)); // TODO: is this copy needed?
            
            if (loglevel == DEBUG){
                printf("roshambo2.cuda: launching kernel for thead %d on device %d\n", thread_id, device_id);
            }
            // loop over dataset molecules is done inside the kernel
            optimize_overlap_gpu(molA_devices[n], molA_type_devices[n], molA_num_atoms_i,      NmolA, 
                                molBs_devices[n],molBs_type_devices[n], molBs_num_atoms_devices[n], NmolBs[n], num_molBss[n], 
                                rmat_devices[n], pmat_devices[n], N_features, scores_devices[n], optim_color, lr_q, lr_t, nsteps, mixing_param, start_mode_method, device_id);

            if(loglevel == DEBUG){
                printf("roshambo2.cuda: completed kernel for thead %d on device %d\n", thread_id, device_id);
            }

            // copy data back
            CUDA_CHECK_ERROR(cudaMemcpy(scores_host, scores_devices[n], scores_sizes[n], cudaMemcpyDeviceToHost));
        } // for n
    } // for i
    

    // clean up
    for(int thread_id=0; thread_id<n_threads; ++thread_id){
        int n = thread_id;
        int device_id = n;

        CUDA_CHECK_ERROR(cudaSetDevice(device_id));
        cudaFree(molA_devices[n]);
        cudaFree(molBs_devices[n]);
        cudaFree(scores_devices[n]);
        cudaFree(molA_type_devices[n]);
        cudaFree(molBs_type_devices[n]);
        cudaFree(molBs_num_atoms_devices[n]);
        cudaFree(rmat_devices[n]);
        cudaFree(pmat_devices[n]);
    }
}


PYBIND11_MODULE(_roshambo2_cuda, m) { 
    m.def("optimize_overlap_color", &optimize_overlap_color, "computes overlap of ref mol A with fit mols B");
}