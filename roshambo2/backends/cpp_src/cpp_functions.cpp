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
#include <cmath>
#include <pybind11/stl.h>
#include <omp.h>
#include <cassert>

#include "cpp_helper_functions.h"

#define DTYPE float

namespace py = pybind11;

enum loglevel {
  OFF,
  INFO,
  DEBUG
};

////////////////////////////////////////////////////////////////////////////////
/// Constants
////////////////////////////////////////////////////////////////////////////////
constexpr int D = 4;
constexpr DTYPE PI = 3.14159265358;
constexpr DTYPE KAPPA = 2.41798793102;
constexpr DTYPE CARBONRADII2 = 1.7*1.7;
constexpr DTYPE A = KAPPA/CARBONRADII2;
const DTYPE CONSTANT = pow(PI/(2*A), 1.5);
constexpr DTYPE EPSILON = 1E-9;
constexpr std::array<int, 3> start_mode_n = {1,4,10};


////////////////////////////////////////////////////////////////////////////////
/// Math helper functions
////////////////////////////////////////////////////////////////////////////////

/// @brief convert quaternion to rotation matrix
/// @param q quaternion[4]
/// @return M matrix[3,3]
std::array<std::array<DTYPE,3>,3> quaternion_to_rotation_matrix(std::array<DTYPE,4> &q){

    std::array<std::array<DTYPE, 3>,3> M;

    // temp variables to make more readable
    auto w = q[0];
    auto x = q[1];
    auto y = q[2];
    auto z = q[3];

    // Compute rotation matrix elements
    M[0][0] = 1 - 2*y*y - 2*z*z;
    M[0][1] = 2*x*y - 2*w*z;
    M[0][2] = 2*x*z + 2*w*y;
    M[1][0] = 2*x*y + 2*w*z;
    M[1][1] = 1 - 2*x*x - 2*z*z;
    M[1][2] = 2*y*z - 2*w*x;
    M[2][0] = 2*x*z - 2*w*y;
    M[2][1] = 2*y*z + 2*w*x;
    M[2][2] = 1 - 2*x*x - 2*y*y;

    return M;
}

/// @brief convert quaternion to rotation matrix
/// @param q quaternion[4]
/// @return M matrix[3,3]
std::array<std::array<DTYPE,3>,3> quaternion_to_rotation_matrix(const double * q){

    std::array<std::array<DTYPE, 3>,3> M;

    // temp variables to make more readable
    auto w = q[0];
    auto x = q[1];
    auto y = q[2];
    auto z = q[3];

    // Compute rotation matrix elements
    M[0][0] = 1 - 2*y*y - 2*z*z;
    M[0][1] = 2*x*y - 2*w*z;
    M[0][2] = 2*x*z + 2*w*y;
    M[1][0] = 2*x*y + 2*w*z;
    M[1][1] = 1 - 2*x*x - 2*z*z;
    M[1][2] = 2*y*z - 2*w*x;
    M[2][0] = 2*x*z - 2*w*y;
    M[2][1] = 2*y*z + 2*w*x;
    M[2][2] = 1 - 2*x*x - 2*y*y;

    return M;
}




/// TODO: docs
std::array<DTYPE,4> axis_angle_to_quat(std::array<DTYPE,4> axis_angle){

    // axis must be normalized. We assume it is
    
    std::array<DTYPE,4> q;

    DTYPE hangle = axis_angle[3]*0.5;
    q[0] = cos(hangle);
    q[1] = axis_angle[0]*sin(hangle);
    q[2] = axis_angle[1]*sin(hangle);
    q[3] = axis_angle[2]*sin(hangle);

    return q;
}

////////////////////////////////////////////////////////////////////////////////
/// For debug
////////////////////////////////////////////////////////////////////////////////

template<typename T, size_t N>
void printArray(const std::array<T, N>& arr) {
    for (const auto& elem : arr) {
        std::cout << elem << " ";
    }
    std::cout << std::endl;
}


////////////////////////////////////////////////////////////////////////////////
/// Transform the query mol to the start mode
////////////////////////////////////////////////////////////////////////////////


/// @brief transform mol by start mode <method> index <index>
/// @param mol
/// @param transformed_mol
/// @param Nmol
/// @param method
/// @param index
/// @return transform (quaternion + vector displacement)
std::array<DTYPE, 7> start_mode_transform(const DTYPE * mol, DTYPE * transformed_mol, int Nmol, int method, int index){

    // store the transformation so we can use it later
    std::array<DTYPE, 7> out = {1,0,0,0,0,0,0}; 

    // initialise the transformed one to be the same as the original one
    for(int i=0; i<Nmol; i++){
        transformed_mol[D*i] = mol[D*i];
        transformed_mol[D*i+1] = mol[D*i+1];
        transformed_mol[D*i+2] = mol[D*i+2];
        transformed_mol[D*i+3] = mol[D*i+3];
    }

    // number of start mode configuration for the chosen method
    int method_n = start_mode_n[method];

    // the methods have a different number of start positions
    // they are listed in here as consts
    switch (method){

        case 0:
        {
            // method 0
            // do nothing
            break;
        }
        
        case 1:
        {
            // list the transformations
            const std::vector<std::array<DTYPE,4> > transforms_1 = {
                {1,0,0,0},
                {1,0,0,PI},
                {0,1,0,PI},
                {0,0,1,PI},
            };

            // get the one we need
            auto my_transform = transforms_1[index];

            // apply the transformation
            auto q = axis_angle_to_quat(my_transform);
            auto M = quaternion_to_rotation_matrix(q);

            for(int i=0; i<Nmol; i++){
                matvec3x3x3(M, &mol[D*i], &transformed_mol[D*i]);
            }

            out[0] = q[0];
            out[1] = q[1];
            out[2] = q[2];
            out[3] = q[3];

            break;
        }

        case 2:
        {
  
            // list the transformations
            const std::vector<std::array<DTYPE,4> > transforms_2 = {
                {1,0,0,0},
                {1,0,0,PI},
                {0,1,0,PI},
                {0,0,1,PI},
                {1,0,0,0.5*PI},
                {0,1,0,0.5*PI},
                {0,0,1,0.5*PI},
                {1,0,0,-0.5*PI},
                {0,1,0,-0.5*PI},
                {0,0,1,-0.5*PI},
            };


            // get the one we need
            auto my_transform = transforms_2[index];
            auto q = axis_angle_to_quat(my_transform);
            auto M = quaternion_to_rotation_matrix(q);

            // apply the transformation
            for(int i=0; i<Nmol; i++){
                matvec3x3x3(M, &mol[D*i], &transformed_mol[D*i]);
            }

            out[0] = q[0];
            out[1] = q[1];
            out[2] = q[2];
            out[3] = q[3];

            break;
        }
        default:
            break;

    } // switch

    return out;
}

////////////////////////////////////////////////////////////////////////////////
/// helper functions for self overlap
////////////////////////////////////////////////////////////////////////////////

DTYPE self_overlap_single(py::array_t<DTYPE> A){
    auto molA  = A.unchecked<2>();

    const DTYPE * ptr_molA = molA.data(0,0);
    int NmolA = molA.shape(0);
   
    auto v = volume(ptr_molA, NmolA, ptr_molA, NmolA);

    return v;

}

std::array<DTYPE,2> self_overlap_single_color(py::array_t<DTYPE> A, py::array_t<int> T, int N, py::array_t<DTYPE> RMAT, py::array_t<DTYPE> PMAT){
    auto molA  = A.unchecked<2>();
    auto molA_type = T.unchecked<1>();
    auto rmat = RMAT.unchecked<2>();
    auto pmat = PMAT.unchecked<2>();

    const DTYPE * ptr_molA = molA.data(0,0);

    int NmolA_real = N;
    int NmolA_color = molA.shape(0) - N;

    const int * ptr_molA_type = molA_type.data(0);
    const DTYPE * ptr_rmat = rmat.data(0,0);
    const DTYPE * ptr_pmat = pmat.data(0,0);
    int N_features = RMAT.shape(0);

    // v normally for real atoms:
    auto v = volume(ptr_molA, NmolA_real, ptr_molA, NmolA_real);


    // color needs a different function
    // offset the pointers to the correct locations
    auto vc = volume_color(&ptr_molA[NmolA_real*D], NmolA_color, &ptr_molA_type[NmolA_real], 
                           &ptr_molA[NmolA_real*D], NmolA_color, &ptr_molA_type[NmolA_real], 
                           ptr_rmat, ptr_pmat, N_features);

    // printf("single: %f %f\n", v, vc);
    return std::array<DTYPE,2>{v, vc};

}

/// @brief optimization function. Entry point to c++ code from python (via pybind)
/// @param A - the padded query coordinate data
/// @param AT - the padded query types
/// @param AN - the number of heavy atoms
/// @param B - the padded target coordinate data
/// @param BT - the padded target types
/// @param BN - the number of heavy atoms
/// @param RMAT - interaction matrix r - linearised square matrix for looking up r for features
/// @param PMAT - interaction matrix p - linearised square matrix for looking up p for features
/// @param V - the output scores
/// @param optim_color - whether to optimise on colors as well as shape
/// @param mixing_param - how to mix the 2 tanimoto values
/// @param lr_q - learning rate for optimising the quaternion
/// @param lr_t - learning rate for optimising the translation
/// @param nsteps - number of optimiser steps
/// @param legacy_combo_gradient - if true, weight the shape and colour
///        gradients by 1/(V_AA+V_BB) as previous releases did, instead of by
///        each channel's dT/dV. Only affects colour/combination modes.
/// @param start_mode_method - 0, 1, or 2 with increasingly more start points
/// @param loglevel - how wordy the output should be
void optimize_overlap_color(py::array_t<DTYPE> A, py::array_t<int> AT, py::array_t<int> AN,
                            py::array_t<DTYPE> B, py::array_t<int> BT, py::array_t<int> BN, 
                            py::array_t<DTYPE> RMAT, py::array_t<DTYPE> PMAT, py::array_t<DTYPE> V, 
                            bool optim_color, DTYPE mixing_param, DTYPE lr_q, DTYPE lr_t, int nsteps,
                            int start_mode_method, int loglevel, bool legacy_combo_gradient){

    // shared arrays/pointers
    auto molAs  = A.unchecked<3>();
    auto molA_type = AT.unchecked<2>();
    auto molA_num_atoms = AN.unchecked<1>();

    auto rmat = RMAT.unchecked<2>();
    auto pmat = PMAT.unchecked<2>();
    const DTYPE * ptr_rmat = rmat.data(0,0);
    const DTYPE * ptr_pmat = pmat.data(0,0);

    int N_features = rmat.shape(0);
    int n_querys = molAs.shape(0);

    

    auto molBs  = B.unchecked<3>();
    auto molB_type = BT.unchecked<2>();
    auto molB_num_atoms = BN.unchecked<1>();
    auto scores = V.mutable_unchecked<3>();

    assert(scores.shape(2) == 20);

    int numThreads = omp_get_max_threads();

    if(loglevel == DEBUG){
        std::cout << "roshambo2.c++: using: " << numThreads << " CPU threads" << std::endl;
        std::cout << "roshambo2.c++: Optimizer settings = { lr_q:" <<lr_q<<" lr_t:"<<lr_t<<" steps: "<<nsteps<< "}"<<std::endl;
    }

    // number of start mode configs for each query
    int n_starts = start_mode_n[start_mode_method];

    // loop over query molecules
    for (int k=0; k < n_querys; ++k){

        // get the query mol
        const DTYPE * ptr_molA_orig = molAs.data(k,0,0);
        int NmolA = molAs.shape(1);
        const int * ptr_molA_type = molA_type.data(k,0);

        // stored for the start mode transformed molA
        std::vector<DTYPE> transformed_molA(ptr_molA_orig, ptr_molA_orig+NmolA*D);
        DTYPE * ptr_molA = transformed_molA.data();


        // loop over start modes
        for (int start_index = 0; start_index < n_starts; ++start_index){
    
            // create a transformed copy of the query mol
            auto start_qr  = start_mode_transform(ptr_molA_orig, ptr_molA, NmolA, start_mode_method, start_index);

            // compute self overlap of A
            int NmolA_real = molA_num_atoms(k);
            int NmolA_color = NmolA -NmolA_real;
            auto self_overlap_A = volume(ptr_molA, NmolA_real, ptr_molA, NmolA_real);
            auto self_overlap_A_color = volume_color(&ptr_molA[NmolA_real*D], NmolA_color, &ptr_molA_type[NmolA_real], 
                                                    &ptr_molA[NmolA_real*D], NmolA_color, &ptr_molA_type[NmolA_real], 
                                                    ptr_rmat, ptr_pmat, N_features);

            // parallel loop over each dataset configuration
            {
            #pragma omp parallel for
            for (long j = 0; j < molBs.shape(0); j++){

                // make a copy of molA and molB to transform
                DTYPE molAT[NmolA*molAs.shape(2)];
                std::memcpy(molAT, ptr_molA, NmolA*molAs.shape(2)*sizeof(DTYPE));

                const DTYPE * ptr_molB = molBs.data(j,0,0);
                int NmolB = molBs.shape(1);
                const int * ptr_molB_type = molB_type.data(j,0);

                DTYPE molBT[NmolB*molBs.shape(2)];
                std::memcpy(molBT, ptr_molB, NmolB*molBs.shape(2)*sizeof(DTYPE));
                int NmolB_real = molB_num_atoms(j);
                int NmolB_color = NmolB -NmolB_real;

				DTYPE these_scores[20];
				single_conformer_optimiser(ptr_molA, molAT, ptr_molA_type,
                                           NmolA, NmolA_real, NmolA_color,
										   self_overlap_A, self_overlap_A_color,
                                           ptr_molB, molBT, ptr_molB_type,
										   NmolB, NmolB_real, NmolB_color, ptr_rmat, ptr_pmat, N_features,
										   optim_color, mixing_param, lr_q, lr_t, nsteps,
										   legacy_combo_gradient, these_scores);

                // Check the previous ones and keep the best. The first start
                // mode always writes, so that every row reports a real pose:
                // otherwise a pair whose every start mode scored <= 0 (or NaN,
                // which compares false) would keep the zero-initialised score
                // row, including a zero quaternion that is not a rotation.
                if (start_index == 0 || these_scores[0] > scores(k,j,0)){
                    scores(k,j,0) = these_scores[0]; // combination tanimoto of shape and color
                    scores(k,j,1) = these_scores[1]; // shape tanimoto
                    scores(k,j,2) = these_scores[2]; // color tanimoto
                    scores(k,j,3) = these_scores[3]; // volume shape
                    scores(k,j,4) = these_scores[4]; // volumes color
                    scores(k,j,5) = these_scores[5]; // self i
                    scores(k,j,6) = these_scores[6]; // self j
                    scores(k,j,7) = these_scores[7]; // self i color
                    scores(k,j,8) = these_scores[8]; // self j color
                    scores(k,j,9)  = these_scores[9]; // optimised quaternion value
                    scores(k,j,10) = these_scores[10]; // optimised quaternion value
                    scores(k,j,11) = these_scores[11]; // optimised quaternion value
                    scores(k,j,12) = these_scores[12]; // optimised quaternion value
                    scores(k,j,13) = these_scores[13]; // optimised translation value
                    scores(k,j,14) = these_scores[14]; // optimised translation value
                    scores(k,j,15) = these_scores[15]; // optimised translation value
                    scores(k,j,16) = start_qr[0];
                    scores(k,j,17) = start_qr[1];
                    scores(k,j,18) = start_qr[2];
                    scores(k,j,19) = start_qr[3];

                } // if

            } // for
            } // omp parallel
        } // for start_index
    } // for k
}

////////////////////////////////////////////////////////////////////////////////
/// wrapper functions for testing framework
////////////////////////////////////////////////////////////////////////////////


DTYPE test_overlap_single(py::array_t<DTYPE> A, py::array_t<DTYPE> B){
    auto molA = A.unchecked<2>();
    auto molB = B.unchecked<2>();

    const DTYPE * ptr_molA = molA.data(0,0);
    int NmolA = molA.shape(0);
    
    const DTYPE * ptr_molB = molB.data(0,0);
    int NmolB = molB.shape(0);

    auto v = volume(ptr_molA, NmolA, ptr_molB, NmolB);

    return v;

}


std::array<DTYPE,7> test_gradient(py::array_t<DTYPE> A, py::array_t<DTYPE> B){
    auto molA = A.unchecked<2>();
    auto molB = B.unchecked<2>();

    const DTYPE * ptr_molA = molA.data(0,0);
    int NmolA = molA.shape(0);
    

    const DTYPE * ptr_molB = molB.data(0,0);
    int NmolB = molB.shape(0);

    auto gq = get_gradient(ptr_molA, NmolA, ptr_molB, NmolB);

    return gq;

}


void test_overlap(py::array_t<DTYPE> A, py::array_t<DTYPE> B, py::array_t<DTYPE> V){
    auto molAs  = A.unchecked<3>(); // x must have ndim = 3; can be non-writeable
    auto molBs  = B.unchecked<3>(); // x must have ndim = 3; can be non-writeable
    auto molV = V.mutable_unchecked<2>();

    int numThreads = omp_get_max_threads();

    std::cout << "using: " << numThreads << " CPU threads" << std::endl;

    for (py::ssize_t i = 0; i < molAs.shape(0); i++){
        #pragma omp parallel for
        for (py::ssize_t j = 0; j < molBs.shape(0); j++){

            const DTYPE * ptr_molA = molAs.data(i,0,0);
            int NmolA = molAs.shape(1);
    

            const DTYPE * ptr_molB = molBs.data(j,0,0);
            int NmolB = molBs.shape(1);

            auto v = volume(ptr_molA, NmolA, ptr_molB, NmolB);

            molV(i,j) = v;
        }
    }
}







////////////////////////////////////////////////////////////////////////////////
/// Bindings for Python
////////////////////////////////////////////////////////////////////////////////
PYBIND11_MODULE(_roshambo2_cpp, m) { 
    m.def("optimize_overlap_color", &optimize_overlap_color, "computes overlap of ref mol A with fit mols B with color",
          py::arg("A"), py::arg("AT"), py::arg("AN"), py::arg("B"), py::arg("BT"), py::arg("BN"),
          py::arg("RMAT"), py::arg("PMAT"), py::arg("V"), py::arg("optim_color"), py::arg("mixing_param"),
          py::arg("lr_q"), py::arg("lr_t"), py::arg("nsteps"), py::arg("start_mode_method"),
          py::arg("loglevel"), py::arg("legacy_combo_gradient") = false);
    m.def("test_overlap", &test_overlap, "computes overlap of ref mol A with fit mols B");
    m.def("test_overlap_single", &test_overlap_single, "single overlap for testing");
    m.def("test_gradient", &test_gradient, "quaternion gradient of single overlap for testing");
}