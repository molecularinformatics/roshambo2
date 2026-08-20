//
// Created by dave on 04/12/2025.
//

#ifndef _ROSHAMBO2_CPP_HELPER_FUNCTIONS_H
#define _ROSHAMBO2_CPP_HELPER_FUNCTIONS_H

using DTYPE = float;

////////////////////////////////////////////////////////////////////////////////
/// Math helper functions
////////////////////////////////////////////////////////////////////////////////

/// @brief division of std::array by scalar
/// @tparam T
/// @tparam N
/// @param arr
/// @param scalar
/// @return
template <typename T, size_t N>
std::array<T, N> operator/(const std::array<T, N> &arr, const T &scalar);

/// @brief product of std::array by scalar
/// @tparam T
/// @tparam N
/// @param arr
/// @param scalar
/// @return
template <typename T, size_t N>
std::array<T, N> operator*(const T &scalar, const std::array<T, N> &arr);

/// @brief addition of std::arrays
/// @tparam T
/// @tparam N
/// @param arr1
/// @param arr2
/// @return
template <typename T, size_t N>
std::array<T, N> operator+(const std::array<T, N> &arr1,
                           const std::array<T, N> &arr2);

/// @brief product of std::arrays
/// @tparam T
/// @tparam N
/// @param arr1
/// @param arr2
/// @return
template <typename T, size_t N>
std::array<T, N> operator*(const std::array<T, N> &arr1,
                           const std::array<T, N> &arr2);

/// @brief 3x3 matrix x vector
/// @param mat matrix[3,3]
/// @param vec 3 vector[3]
/// @param result 3 vector[3]
void matvec3x3x3(DTYPE mat[][3], const DTYPE *vec, DTYPE *result);

/// @brief 3x3 matrix x vector
/// @param mat matrix[3,3]
/// @param vec 3 vector[3]
/// @return result 3 vector[3]
std::array<DTYPE, 3> matvec3x3x3(const std::array<std::array<DTYPE, 3>, 3> &mat,
                                 const DTYPE *vec);

/// @brief 3x3 matrix x vector
/// @param mat matrix[3,3]
/// @param vec 3 vector[3]
/// @param result 3 vector[3]
void matvec3x3x3(const std::array<std::array<DTYPE, 3>, 3> &mat,
                 const DTYPE *vec, DTYPE *result);

/// @brief convert quaternion to rotation matrix
/// @param q quaternion[4]
/// @param M matrix[3,3]
void quaternion_to_rotation_matrix(std::array<DTYPE, 4> &q, DTYPE M[3][3]);

/// @brief transform molB by mat
/// @param mat matrix[3,3]
/// @param molB molecular coordinates[N,3]
/// @param molBT molecular coordinates[N,3]
/// @param NmolB N
void transform(DTYPE mat[][3], const DTYPE *molB, DTYPE *molBT, int NmolB);

/// @brief translate molA to molAT by t
/// @param t vector[3]
/// @param molA molecular coordinates[N,3]
/// @param molAT molecular coordinates[N,3]
/// @param NmolA N
void translate(DTYPE t[3], const DTYPE *molA, DTYPE *molAT, int NmolA);

///////////////////////////////////////////////////////////////////////////////
/// Volume functions
///////////////////////////////////////////////////////////////////////////////

/// @brief shape overlap volume of molA and molB
/// @param molA
/// @param NmolA
/// @param molB
/// @param NmolB
/// @return volume
DTYPE volume(const DTYPE *molA, int NmolA, const DTYPE *molB, int NmolB);

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
DTYPE volume_color(const DTYPE *molA, int NmolA, const int *molA_type,
                   const DTYPE *molB, int NmolB, const int *molB_type,
                   const DTYPE *rmat, const DTYPE *pmat, int N_features);

////////////////////////////////////////////////////////////////////////////////
/// gradient functions
////////////////////////////////////////////////////////////////////////////////

/// @brief compute overlap gradients of molA and molB w.r.t molB
/// @param molA
/// @param NmolA
/// @param molB
/// @param NmolB
/// @param volume_out - if not null, receives the overlap volume at this pose
/// @return array[7] containing the quaternion gradients and the position
/// gradients
std::array<DTYPE, 7> get_gradient(const DTYPE *molA, int NmolA,
                                  const DTYPE *molB, int NmolB,
                                  DTYPE *volume_out = nullptr);

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
/// @param volume_out - if not null, receives the colour overlap volume at this
/// pose
/// @return array[7] containing the quaternion gradients and the position
/// gradients
std::array<DTYPE, 7> get_gradient_color(const DTYPE *molA, int NmolA,
                                        const int *molA_type, const DTYPE *molB,
                                        int NmolB, const int *molB_type,
                                        const DTYPE *rmat, const DTYPE *pmat,
                                        int N_features,
                                        DTYPE *volume_out = nullptr);

////////////////////////////////////////////////////////////////////////////////
/// Optimization functions
////////////////////////////////////////////////////////////////////////////////

/// @brief per-channel weights for mixing the shape and colour gradients
///
/// The reported score is a blend of two Tanimotos, so its gradient weights each
/// channel by dT/dV = (V_AA + V_BB) / (V_AA + V_BB - V)^2 at the current pose.
/// @param legacy - if true, use the previous weight 1/(V_AA + V_BB) instead
void combo_gradient_weights(DTYPE self_a, DTYPE self_b, DTYPE v,
                            DTYPE self_a_color, DTYPE self_b_color,
                            DTYPE v_color, bool legacy, DTYPE *w_shape,
                            DTYPE *w_color);

// Adagrad optimization step
void adagrad_step(std::array<DTYPE, 4> &q, std::array<DTYPE, 3> &t,
                  std::array<DTYPE, 7> g, std::array<DTYPE, 7> &cache,
                  DTYPE lr_q, DTYPE lr_t);

/// @brief Do the overlay for a single conformer of A against a single conformer
/// of B.  The output in scores is the rotation and translation that moves B
/// to optimise its score with A.
/// @param molA - the query molecule as 1D array of 4 * NmolA entries. Each
/// block of 4 is the coords and w parameter
/// @param molAT - the working copy of A
/// @param molA_type - the features types for molecule A
/// @param NmolA - the number of atoms and features in A
/// @param NmolA_real - the number of atoms in A
/// @param NmolA_color - the number of features in A
/// @param self_overlap_A - overlap volume of A with itself
/// @param self_overlap_A_color - color overlap of A with itself
/// @param molB - the target molecule as 1D array of 4 * NmolB entries. Each
/// block of 4 is the coords and w parameter
/// @param molBT - the working copy of B - the final overlaid coordinates will
/// be left here
/// @param molB_type - the features types for molecule B
/// @param NmolB - the number of atoms and features in B
/// @param NmolB_real - the number of atoms in B
/// @param NmolB_color - the number of features in B
/// @param rmat - interaction matrix r - linearised square matrix for looking up
/// r for features
/// @param pmat - interaction matrix p - linearised square matrix for looking up
/// p for features
/// @param N_features - the number of feature types for looking up values in
/// rmat and pmat
/// @param optim_color - whether to optimise on colors as well as shape
/// @param mixing_param - how to mix the 2 tanimoto values
/// @param lr_q - scale factor for optimising the quaternion (?)
/// @param lr_t - scale factor for optimising the translation (?)
/// @param nsteps - number of optimiser steps
/// @return scores - the output scores and transformation to reproduce the
/// overlay - an array of size 20. Only the first 16 are used here. They are:
/// 0 - the combo score
/// 1 - the shape tanimoto
/// 2 - the color tanimoto
/// 3 - the overlap volume
/// 4 - the color overlap volume
/// 5 - the volume of A
/// 6 - the volume of B
/// 7 - the color volume of A
/// 8 - the color volume of B
/// 9-12 - the quaternion to rotate B onto A. Applied first.
/// 13-15 - the translation to move B onto A. Applied second.
/// 16-19 - returned as zeros.
void single_conformer_optimiser(
    const DTYPE *molA, DTYPE *molAT, const int *molA_type, int NmolA,
    int NmolA_real, int NmolA_color, DTYPE self_overlap_A,
    DTYPE self_overlap_A_color, const DTYPE *molB, DTYPE *molBT,
    const int *molB_type, int NmolB, int NmolB_real, int NmolB_color,
    const DTYPE *rmat, const DTYPE *pmat, int N_features, bool optim_color,
    DTYPE mixing_param, DTYPE lr_q, DTYPE lr_t, int nsteps,
    bool legacy_combo_gradient, DTYPE *scores);

#endif  //_ROSHAMBO2_CPP_HELPER_FUNCTIONS_H