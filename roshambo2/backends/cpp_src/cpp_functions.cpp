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
const int D = 4;
const DTYPE PI = 3.14159265358;
const DTYPE KAPPA = 2.41798793102;
const DTYPE CARBONRADII2 = 1.7*1.7;
const DTYPE A = KAPPA/CARBONRADII2;
const DTYPE CONSTANT = pow(PI/(2*A), 1.5);
const DTYPE EPSILON = 1E-9;
const std::array<int, 3> start_mode_n = {1,4,10};


////////////////////////////////////////////////////////////////////////////////
/// Math helper functions
////////////////////////////////////////////////////////////////////////////////

/// @brief 3x3 matrix x vector
/// @param mat matrix[3,3]
/// @param vec 3 vector[3]
/// @param result 3 vector[3]
void matvec3x3x3(DTYPE mat[][3], const DTYPE * vec, DTYPE * result){

    // only transforms the first 3 values in the array! The 4th value is not a 
    // coordinate so we do not transform it

    result[0] = 0.0;
    result[1] = 0.0;
    result[2] = 0.0;

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            result[i] += mat[i][j] * vec[j];
        }
    }
}

/// @brief 3x3 matrix x vector
/// @param mat matrix[3,3]
/// @param vec 3 vector[3]
/// @param result 3 vector[3]
void matvec3x3x3(const std::array<std::array<DTYPE,3>,3> &mat, const DTYPE * vec, DTYPE * result){

    // only transforms the first 3 values in the array! The 4th value is not a 
    // coordinate so we do not transform it

    result[0] = 0.0;
    result[1] = 0.0;
    result[2] = 0.0;

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            result[i] += mat[i][j] * vec[j];
        }
    }
}

/// @brief 3x3 matrix x vector
/// @param mat matrix[3,3]
/// @param vec 3 vector[3]
/// @return result 3 vector[3]
std::array<DTYPE,3> matvec3x3x3(const std::array<std::array<DTYPE,3>,3> &mat, const DTYPE * vec){

    // only transforms the first 3 values in the array! The 4th value is not a 
    // coordinate so we do not transform it
    std::array<DTYPE,3> result;
    result[0] = 0.0;
    result[1] = 0.0;
    result[2] = 0.0;

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            result[i] += mat[i][j] * vec[j];
        }
    }

    return result;
}




/// @brief convert quaternion to rotation matrix
/// @param q quaternion[4]
/// @param M matrix[3,3]
void quaternion_to_rotation_matrix(std::array<DTYPE,4> &q, DTYPE M[3][3]){

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
}


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




/// @brief translate molA to molAT by t
/// @param t vector[3]
/// @param molA molecular coordinates[N,3]
/// @param molAT molecular coordinates[N,3]
/// @param NmolA N
void translate(DTYPE t[3], const DTYPE * molA, DTYPE * molAT, int NmolA){

    // only transforms the first 3 values in the array! The 4th value is not a 
    // coordinate so we do not transform it

    for(int i=0;i<NmolA;i++){

        auto x = &molA[i*D];
        auto y = &molAT[i*D];

        y[0] = x[0] + t[0];
        y[1] = x[1] + t[1];
        y[2] = x[2] + t[2];
    }
}


/// @brief transform molB by mat
/// @param mat matrix[3,3]
/// @param molB molecular coordinates[N,3]
/// @param molBT molecular coordinates[N,3]
/// @param NmolB N
void transform(DTYPE mat[][3], const DTYPE * molB, DTYPE * molBT, int NmolB){

    // only transforms the first 3 values in the array! The 4th value is not a 
    // coordinate so we do not transform it

    for(int i=0;i<NmolB;i++){

        auto x = &molB[i*D];
        auto y = &molBT[i*D];

        matvec3x3x3(mat, x, y);
    }
}


/// @brief division of std::array by scalar
/// @tparam T 
/// @tparam N 
/// @param arr 
/// @param scalar 
/// @return 
template<typename T, size_t N>
std::array<T, N> operator/(const std::array<T, N>& arr, const T& scalar) {
    std::array<T, N> result;
    for (size_t i = 0; i < N; ++i) {
        result[i] = arr[i] / scalar;
    }
    return result;
}

/// @brief product of std::array by scalar
/// @tparam T 
/// @tparam N 
/// @param arr 
/// @param scalar 
/// @return 
template<typename T, size_t N>
std::array<T, N> operator*(const T& scalar, const std::array<T, N>& arr){
    std::array<T, N> result;
    for (size_t i = 0; i < N; ++i) {
        result[i] = scalar*arr[i]; 
    }
    return result;
}

/// @brief addition of std::arrays
/// @tparam T 
/// @tparam N 
/// @param arr1 
/// @param arr2
/// @return 
template<typename T, size_t N>
std::array<T, N> operator+(const std::array<T, N>& arr1, const std::array<T, N>& arr2) {
    std::array<T, N> result;
    for (size_t i = 0; i < N; ++i) {
        result[i] = arr1[i] + arr2[i];
    }
    return result;
}

/// @brief product of std::arrays
/// @tparam T 
/// @tparam N 
/// @param arr1 
/// @param arr2
/// @return 
template<typename T, size_t N>
std::array<T, N> operator*(const std::array<T, N>& arr1, const std::array<T, N>& arr2) {
    std::array<T, N> result;
    for (size_t i = 0; i < N; ++i) {
        result[i] = arr1[i]*arr2[i];
    }
    return result;
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



///////////////////////////////////////////////////////////////////////////////
/// Volume functions
///////////////////////////////////////////////////////////////////////////////


/// @brief shape overlap volume of molA and molB
/// @param molA 
/// @param NmolA 
/// @param molB 
/// @param NmolB 
/// @return volume 
DTYPE volume(const DTYPE * molA, int NmolA, const DTYPE * molB, int NmolB){

    DTYPE V = 0.0;

    for(int i=0; i < NmolA; i++){
        for( int j=0; j<NmolB; j++){
            
            DTYPE dx = molA[i*D]   - molB[j*D];
            DTYPE dy = molA[i*D+1] - molB[j*D+1];
            DTYPE dz = molA[i*D+2] - molB[j*D+2];

            DTYPE d2 = dx*dx + dy*dy + dz*dz;
        
            auto a1 = A;  // left easy to change to not doing all-carbon radii
            auto a2 = A;
            
            DTYPE wa = molA[i*D+3]; // wa,wb == zero means it is a padded atom 
            DTYPE wb = molB[j*D+3];

            DTYPE kij = exp(-a1*a2*d2/(a1+a2))*wa*wb;

            DTYPE vij = 8*kij*CONSTANT;

            V += vij;
        }
    }

    return V;
}

/// @brief color overlap volume of molA and molB
/// @param molA 
/// @param NmolA 
/// @param molA_type 
/// @param molB 
/// @param NmolB 
/// @param molB_type 
/// @param rmat 
/// @param pmat 
/// @param N_features 
/// @return volume
DTYPE volume_color(const DTYPE * molA, int NmolA, const int * molA_type, 
                   const DTYPE * molB, int NmolB, const int * molB_type, 
                   const DTYPE * rmat, const DTYPE * pmat, int N_features){

    DTYPE V = 0.0;

    for(int i=0; i < NmolA; i++){
        int ta = molA_type[i];
        if (ta==0) break; // padded atoms are at the end and have type==0

        for( int j=0; j<NmolB; j++){
            int tb = molB_type[j];
            if (tb==0) break;
            
            DTYPE dx = molA[i*D]   - molB[j*D];
            DTYPE dy = molA[i*D+1] - molB[j*D+1];
            DTYPE dz = molA[i*D+2] - molB[j*D+2];

            DTYPE d2 = dx*dx + dy*dy + dz*dz;
        
            auto a = rmat[ta*N_features+tb];
            auto p = pmat[ta*N_features+tb];
        
            
            DTYPE wa = molA[i*D+3];
            DTYPE wb = molB[j*D+3];

            DTYPE kij = exp(-a*a*d2/(a+a))*wa*wb;

            double constant = pow(PI/(2*a), 1.5);

            DTYPE vij = p*p*kij*constant;

            V += vij;
        }
    }

    return V;
}


////////////////////////////////////////////////////////////////////////////////
/// gradient functions
////////////////////////////////////////////////////////////////////////////////

/// @brief compute overlap gradients of molA and molB w.r.t molB
/// @param molA 
/// @param NmolA 
/// @param molB 
/// @param NmolB 
/// @return array[7] containing the quaternion gradients and the position gradients
std::array<DTYPE,7> get_gradient(const DTYPE * molA, int NmolA, const DTYPE * molB, int NmolB){

    std::array<DTYPE,7> grad = {0.0,0.0,0.0,0.0,0.0,0.0,0.0};

    for(int i=0; i < NmolA; i++){
        for( int j=0; j<NmolB; j++){
            
            DTYPE dx = molA[i*D]   - molB[j*D];
            DTYPE dy = molA[i*D+1] - molB[j*D+1];
            DTYPE dz = molA[i*D+2] - molB[j*D+2];

            DTYPE d2 = dx*dx + dy*dy + dz*dz;
        
            auto a1 = A; // left easy to change to not doing all-carbon radii
            auto a2 = A;

            DTYPE wa = molA[i*D+3];
            DTYPE wb = molB[j*D+3];

            DTYPE kij = exp(-a1*a2*d2/(a1+a2))*wa*wb;

            DTYPE vij = 8*kij*CONSTANT;

            DTYPE x = molB[j*D];
            DTYPE y = molB[j*D+1];
            DTYPE z = molB[j*D+2];

            DTYPE sks[3][3] = {
                {0,-2*z, 2*y},
                {2*z,0,-2*x},
                {-2*y,2*x,0}
            };

            DTYPE delta[3] = {dx, dy, dz};

            DTYPE mv[3] = {0.0,0.0,0.0};

            matvec3x3x3(sks, delta, mv);

            grad[1] += -2.0*(a1*a2)/(a1+a2)*vij*mv[0];
            grad[2] += -2.0*(a1*a2)/(a1+a2)*vij*mv[1];
            grad[3] += -2.0*(a1*a2)/(a1+a2)*vij*mv[2];

            // x,y,z
            grad[4] += -2.0*(a1*a2)/(a1+a2)*vij*delta[0];
            grad[5] += -2.0*(a1*a2)/(a1+a2)*vij*delta[1];
            grad[6] += -2.0*(a1*a2)/(a1+a2)*vij*delta[2];
        }
    }

    return grad;
}


/// @brief compute overlap color gradients of molA and molB w.r.t molB
/// @param molA 
/// @param NmolA 
/// @param molA_type 
/// @param molB 
/// @param NmolB 
/// @param molB_type 
/// @param rmat 
/// @param pmat 
/// @param N_features 
/// @return array[7] containing the quaternion gradients and the position gradients
std::array<DTYPE,7> get_gradient_color(const DTYPE * molA, int NmolA, const int * molA_type, 
                                       const DTYPE * molB, int NmolB, const int * molB_type, 
                                       const DTYPE * rmat, const DTYPE * pmat, int N_features){

    std::array<DTYPE,7> grad = {0.0,0.0,0.0,0.0,0.0,0.0,0.0};

    DTYPE V=0.0;

    for(int i=0; i < NmolA; i++){
        int ta = molA_type[i];
        if (ta==0) break;

        for( int j=0; j<NmolB; j++){
            int tb = molB_type[j];
            if (tb==0) break;
    
            DTYPE dx = molA[i*D]   - molB[j*D];
            DTYPE dy = molA[i*D+1] - molB[j*D+1];
            DTYPE dz = molA[i*D+2] - molB[j*D+2];

            DTYPE d2 = dx*dx + dy*dy + dz*dz;
        
            auto a = rmat[ta*N_features+tb];
            auto p = pmat[ta*N_features+tb];

            DTYPE wa = molA[i*D+3];
            DTYPE wb = molB[j*D+3];
            
            DTYPE kij = exp(-a*a*d2/(a+a))*wa*wb;

            double constant = pow(PI/(2*a), 1.5);

            DTYPE vij = p*p*kij*constant;

            DTYPE x = molB[j*D];
            DTYPE y = molB[j*D+1];
            DTYPE z = molB[j*D+2];

            DTYPE sks[3][3] = {
                {0,-2*z, 2*y},
                {2*z,0,-2*x},
                {-2*y,2*x,0}
            };

            DTYPE delta[3] = {dx, dy, dz};

            DTYPE mv[3] = {0.0,0.0,0.0};

            matvec3x3x3(sks, delta, mv);

            grad[1] += -2.0*(a*a)/(a+a)*vij*mv[0];
            grad[2] += -2.0*(a*a)/(a+a)*vij*mv[1];
            grad[3] += -2.0*(a*a)/(a+a)*vij*mv[2];

            // x,y,z
            grad[4] += -2.0*(a*a)/(a+a)*vij*delta[0];
            grad[5] += -2.0*(a*a)/(a+a)*vij*delta[1];
            grad[6] += -2.0*(a*a)/(a+a)*vij*delta[2];
        }
    }
    return grad;
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



////////////////////////////////////////////////////////////////////////////////
/// Optimization functions
////////////////////////////////////////////////////////////////////////////////


// Adagrad optimization step
void adagrad_step(std::array<DTYPE,4> &q, std::array<DTYPE,3> &t,  std::array<DTYPE,7> g, 
                  std::array<DTYPE,7> &cache, DTYPE lr_q, DTYPE lr_t) {

    cache =  cache + g*g;

    q[0] -= lr_q*g[0]/(sqrt(cache[0])+EPSILON);
    q[1] -= lr_q*g[1]/(sqrt(cache[1])+EPSILON);
    q[2] -= lr_q*g[2]/(sqrt(cache[2])+EPSILON);
    q[3] -= lr_q*g[3]/(sqrt(cache[3])+EPSILON);

    t[0] -= lr_t*g[4]/(sqrt(cache[4])+EPSILON);
    t[1] -= lr_t*g[5]/(sqrt(cache[5])+EPSILON);
    t[2] -= lr_t*g[6]/(sqrt(cache[6])+EPSILON);
}





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
void optimize_overlap_color(py::array_t<DTYPE> A, py::array_t<int> AT, py::array_t<int> AN,
                            py::array_t<DTYPE> B, py::array_t<int> BT, py::array_t<int> BN, 
                            py::array_t<DTYPE> RMAT, py::array_t<DTYPE> PMAT, py::array_t<DTYPE> V, 
                            bool optim_color, DTYPE mixing_param, DTYPE lr_q, DTYPE lr_t, int nsteps,
                            int start_mode_method, int loglevel){


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


    DTYPE lr = lr_q;
    DTYPE lrt = lr_t;


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


                std::array<DTYPE,4> q = {1.0,0.0,0.0,0.0}; // initial q
                std::array<DTYPE,3> t = {0.0,0.0,0.0}; // initial t

                const DTYPE * ptr_molB = molBs.data(j,0,0);
                int NmolB = molBs.shape(1);
                const int * ptr_molB_type = molB_type.data(j,0);

                // make a copy of molA and molB to transform
                DTYPE molBT[NmolB*molBs.shape(2)];
                std::memcpy(molBT, ptr_molB, NmolB*molBs.shape(2)*sizeof(DTYPE));

                DTYPE molAT[NmolA*molAs.shape(2)];
                std::memcpy(molAT, ptr_molA, NmolA*molAs.shape(2)*sizeof(DTYPE));


                DTYPE M[3][3] = {{0.0,0.0,0.0},
                                {0.0,0.0,0.0},
                                {0.0,0.0,0.0}};



                int NmolB_real = molB_num_atoms(j);
                int NmolB_color = NmolB -NmolB_real;

                auto self_overlap_B = volume(molBT, NmolB_real, molBT, NmolB_real);
                auto self_overlap_B_color = volume_color(&molBT[NmolB_real*D],    NmolB_color, &ptr_molB_type[NmolB_real], 
                                                    &molBT[NmolB_real*D],    NmolB_color, &ptr_molB_type[NmolB_real], 
                                                    ptr_rmat, ptr_pmat, N_features);


                std::array<DTYPE,7> cache = {0.0,0.0,0.0,0.0,0.0,0.0,0.0};
                
                // optimization loop
                for(int m=0; m<nsteps; m++){


                    // important: we must compute gradients in the reference frame of molB


                    // 1. first we rotate molB

                    quaternion_to_rotation_matrix(q, M);

                    transform(M, ptr_molB, molBT, NmolB);


                    // 2. translate molA by -t
                    DTYPE mt[3] = {-t[0],-t[1],-t[2]};
                    translate(mt, ptr_molA, molAT, NmolA);



                    auto g = get_gradient(molAT, NmolA_real, molBT, NmolB_real);
                    g = g/(self_overlap_A+self_overlap_B); // normalize


                    if(optim_color){
                        auto g_c = get_gradient_color(&molAT[NmolA_real*D],    NmolA_color, &ptr_molA_type[NmolA_real], 
                                                &molBT[NmolB_real*D],    NmolB_color, &ptr_molB_type[NmolB_real], 
                                                ptr_rmat, ptr_pmat, N_features);
                    
                        // normalize and combine
                        g_c = g_c/(self_overlap_A_color+self_overlap_B_color);
                        auto g_combo = (1-mixing_param)*g + mixing_param*g_c;
                        adagrad_step(q, t, g_combo, cache, lr_q, lr_t);

                    }else{
                        adagrad_step(q, t, g, cache, lr_q, lr_t);
                        
                    }

                    //normalize q so that it is a unit quaternion
                    DTYPE magq = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] * q[3]*q[3]);

                    q = q/magq;
                }

                // get volume with transformed B
                quaternion_to_rotation_matrix(q, M);
                transform(M, ptr_molB, molBT, NmolB);
                translate(t.data(), molBT, molBT, NmolB);

                
                auto vol = volume(ptr_molA, NmolA_real, molBT, NmolB_real);
                // printf("v: %f\n", vol);


                // printf("here1\n");
                // printf("N,N: %d, %d\n", NmolA, NmolB);
                // printf("N,N: %d, %d\n", NmolA_real, NmolB_real);
                // printf("N,N: %d, %d\n", NmolA_color, NmolB_color);


                auto vc = volume_color(&ptr_molA[NmolA_real*D], NmolA_color, &ptr_molA_type[NmolA_real], 
                                    &molBT[NmolB_real*D],    NmolB_color, &ptr_molB_type[NmolB_real], 
                                    ptr_rmat, ptr_pmat, N_features);
                // printf("here2\n");
                // printf("v color: %f\n", vc);


                // compute tanimoto scores
                auto ts = vol / (self_overlap_A + self_overlap_B - vol);
                
                DTYPE tc = 0.0;
                if(optim_color){
                    tc = vc / (self_overlap_A_color + self_overlap_B_color - vc);
                }
                
                // we use the mixing param to weight the tanimotos
                // this is the objective function we have optimized
                auto total = (ts *(1-mixing_param) + tc*mixing_param);

                // printf("scores: %f %f %f\n", ts, tc, total);

            
                // check the previous ones and keep the best
                if (total > scores(k,j,0)){
                    scores(k,j,0) = total; // combination tanimoto of shape and color
                    scores(k,j,1) = ts; // shape tanimoto
                    scores(k,j,2) = tc; // color tanimoto
                    scores(k,j,3) = vol; // volume shape
                    scores(k,j,4) = vc; // volumes color
                    scores(k,j,5) = self_overlap_A; // self i
                    scores(k,j,6) = self_overlap_B; // self j
                    scores(k,j,7) = self_overlap_A_color; // self i color
                    scores(k,j,8) = self_overlap_B_color; // self j color
                    scores(k,j,9)  = q[0];
                    scores(k,j,10) = q[1];
                    scores(k,j,11) = q[2];
                    scores(k,j,12) = q[3];
                    scores(k,j,13) = t[0];
                    scores(k,j,14) = t[1];
                    scores(k,j,15) = t[2];
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
    m.def("optimize_overlap_color", &optimize_overlap_color, "computes overlap of ref mol A with fit mols B with color");
    m.def("test_overlap", &test_overlap, "computes overlap of ref mol A with fit mols B");
    m.def("test_overlap_single", &test_overlap_single, "single overlap for testing");
    m.def("test_gradient", &test_gradient, "quaternion gradient of single overlap for testing");
}