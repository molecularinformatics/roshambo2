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
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include <array>
#include <cmath>
#include <iostream>

using DTYPE = float;

////////////////////////////////////////////////////////////////////////////////
/// Constants
////////////////////////////////////////////////////////////////////////////////
constexpr int D = 4;
constexpr DTYPE PI = 3.14159265358;
constexpr DTYPE KAPPA = 2.41798793102;
constexpr DTYPE CARBONRADII2 = 1.7 * 1.7;
constexpr DTYPE A = KAPPA / CARBONRADII2;
const DTYPE CONSTANT = pow(PI / (2 * A), 1.5);
constexpr DTYPE EPSILON = 1E-9;

////////////////////////////////////////////////////////////////////////////////
/// For debug
////////////////////////////////////////////////////////////////////////////////

template <typename T, size_t N>
void printArray(const std::array<T, N> &arr) {
  for (const auto &elem : arr) {
    std::cout << elem << " ";
  }
  std::cout << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Math helper functions
////////////////////////////////////////////////////////////////////////////////

template <typename T, size_t N>
std::array<T, N> operator/(const std::array<T, N> &arr, const T &scalar) {
  std::array<T, N> result;
  for (size_t i = 0; i < N; ++i) {
    result[i] = arr[i] / scalar;
  }
  return result;
}

template <typename T, size_t N>
std::array<T, N> operator*(const T &scalar, const std::array<T, N> &arr) {
  std::array<T, N> result;
  for (size_t i = 0; i < N; ++i) {
    result[i] = scalar * arr[i];
  }
  return result;
}

template <typename T, size_t N>
std::array<T, N> operator+(const std::array<T, N> &arr1,
                           const std::array<T, N> &arr2) {
  std::array<T, N> result;
  for (size_t i = 0; i < N; ++i) {
    result[i] = arr1[i] + arr2[i];
  }
  return result;
}

template <typename T, size_t N>
std::array<T, N> operator*(const std::array<T, N> &arr1,
                           const std::array<T, N> &arr2) {
  std::array<T, N> result;
  for (size_t i = 0; i < N; ++i) {
    result[i] = arr1[i] * arr2[i];
  }
  return result;
}

void matvec3x3x3(DTYPE mat[][3], const DTYPE *vec, DTYPE *result) {
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

std::array<DTYPE, 3> matvec3x3x3(const std::array<std::array<DTYPE, 3>, 3> &mat,
                                 const DTYPE *vec) {
  // only transforms the first 3 values in the array! The 4th value is not a
  // coordinate so we do not transform it
  std::array<DTYPE, 3> result;
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

void matvec3x3x3(const std::array<std::array<DTYPE, 3>, 3> &mat,
                 const DTYPE *vec, DTYPE *result) {
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

void quaternion_to_rotation_matrix(std::array<DTYPE, 4> &q, DTYPE M[3][3]) {
  // temp variables to make more readable
  auto w = q[0];
  auto x = q[1];
  auto y = q[2];
  auto z = q[3];

  // Compute rotation matrix elements
  M[0][0] = 1 - 2 * y * y - 2 * z * z;
  M[0][1] = 2 * x * y - 2 * w * z;
  M[0][2] = 2 * x * z + 2 * w * y;
  M[1][0] = 2 * x * y + 2 * w * z;
  M[1][1] = 1 - 2 * x * x - 2 * z * z;
  M[1][2] = 2 * y * z - 2 * w * x;
  M[2][0] = 2 * x * z - 2 * w * y;
  M[2][1] = 2 * y * z + 2 * w * x;
  M[2][2] = 1 - 2 * x * x - 2 * y * y;
}

void transform(DTYPE mat[][3], const DTYPE *molB, DTYPE *molBT, int NmolB) {
  // only transforms the first 3 values in the array! The 4th value is not a
  // coordinate so we do not transform it

  for (int i = 0; i < NmolB; i++) {
    auto x = &molB[i * D];
    auto y = &molBT[i * D];

    matvec3x3x3(mat, x, y);
  }
}

void translate(DTYPE t[3], const DTYPE *molA, DTYPE *molAT, int NmolA) {
  // only transforms the first 3 values in the array! The 4th value is not a
  // coordinate so we do not transform it

  for (int i = 0; i < NmolA; i++) {
    auto x = &molA[i * D];
    auto y = &molAT[i * D];

    y[0] = x[0] + t[0];
    y[1] = x[1] + t[1];
    y[2] = x[2] + t[2];
  }
}

///////////////////////////////////////////////////////////////////////////////
/// Volume functions
///////////////////////////////////////////////////////////////////////////////

/// @brief Calculate the atom volume overlap between mols A and B
/// @param molA - 1D array of 4 * NmolA entries. Each block of 4 is the coords
/// and w parameter
/// @param NmolA - number of atoms in A
/// @param molB - 1D array of 4 * NmolB entries. Each block of 4 is the coords
/// and w parameter
/// @param NmolA - number of atoms in B
DTYPE volume(const DTYPE *molA, int NmolA, const DTYPE *molB, int NmolB) {
  DTYPE V = 0.0;

  for (int i = 0; i < NmolA; i++) {
    for (int j = 0; j < NmolB; j++) {
      DTYPE dx = molA[i * D] - molB[j * D];
      DTYPE dy = molA[i * D + 1] - molB[j * D + 1];
      DTYPE dz = molA[i * D + 2] - molB[j * D + 2];

      DTYPE d2 = dx * dx + dy * dy + dz * dz;

      auto a1 = A;  // left easy to change to not doing all-carbon radii
      auto a2 = A;

      DTYPE wa = molA[i * D + 3];  // wa,wb == zero means it is a padded atom
      DTYPE wb = molB[j * D + 3];

      DTYPE kij = exp(-a1 * a2 * d2 / (a1 + a2)) * wa * wb;

      DTYPE vij = 8 * kij * CONSTANT;

      V += vij;
    }
  }

  return V;
}

/// @brief Calculate the feature/color volume overlap between mols A and B.
/// @param molA - 1D array of 4 * NmolA entries. Each block of 4 is the coords
/// and w parameter.
/// @param NmolA - number of features in A
/// @param molA_type - types of features
/// @param molB - 1D array of 4 * NmolB entries. Each block of 4 is the coords
/// and w parameter
/// @param NmolA - number of features in B
/// @param molB_type - types of features
/// @param rmat - interaction matrix r - linearised square matrix for looking up
/// r for features
/// @param pmat - interaction matrix p - linearised square matrix for looking up
/// p for features
/// @param N_features - the number of feature types for looking up values in
/// rmat and pmat
DTYPE volume_color(const DTYPE *molA, int NmolA, const int *molA_type,
                   const DTYPE *molB, int NmolB, const int *molB_type,
                   const DTYPE *rmat, const DTYPE *pmat, int N_features) {
  DTYPE V = 0.0;

  for (int i = 0; i < NmolA; i++) {
    int ta = molA_type[i];
    if (ta == 0) break;  // padded atoms are at the end and have type==0
    for (int j = 0; j < NmolB; j++) {
      int tb = molB_type[j];
      if (tb == 0) break;

      DTYPE dx = molA[i * D] - molB[j * D];
      DTYPE dy = molA[i * D + 1] - molB[j * D + 1];
      DTYPE dz = molA[i * D + 2] - molB[j * D + 2];

      DTYPE d2 = dx * dx + dy * dy + dz * dz;

      auto a = rmat[ta * N_features + tb];
      auto p = pmat[ta * N_features + tb];

      // A zero interaction width means this feature pair does not interact,
      // not that the Gaussian has zero width: its contribution is zero and
      // evaluating it would divide by zero.
      if (a == 0) continue;

      DTYPE wa = molA[i * D + 3];
      DTYPE wb = molB[j * D + 3];

      DTYPE kij = exp(-a * a * d2 / (a + a)) * wa * wb;

      double constant = pow(PI / (2 * a), 1.5);

      DTYPE vij = p * p * kij * constant;

      V += vij;
    }
  }

  return V;
}

////////////////////////////////////////////////////////////////////////////////
/// gradient functions
////////////////////////////////////////////////////////////////////////////////

std::array<DTYPE, 7> get_gradient(const DTYPE *molA, int NmolA,
                                  const DTYPE *molB, int NmolB) {
  std::array<DTYPE, 7> grad = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  for (int i = 0; i < NmolA; i++) {
    for (int j = 0; j < NmolB; j++) {
      DTYPE dx = molA[i * D] - molB[j * D];
      DTYPE dy = molA[i * D + 1] - molB[j * D + 1];
      DTYPE dz = molA[i * D + 2] - molB[j * D + 2];

      DTYPE d2 = dx * dx + dy * dy + dz * dz;

      auto a1 = A;  // left easy to change to not doing all-carbon radii
      auto a2 = A;

      DTYPE wa = molA[i * D + 3];
      DTYPE wb = molB[j * D + 3];

      DTYPE kij = exp(-a1 * a2 * d2 / (a1 + a2)) * wa * wb;

      DTYPE vij = 8 * kij * CONSTANT;

      DTYPE x = molB[j * D];
      DTYPE y = molB[j * D + 1];
      DTYPE z = molB[j * D + 2];

      DTYPE sks[3][3] = {
          {0, -2 * z, 2 * y}, {2 * z, 0, -2 * x}, {-2 * y, 2 * x, 0}};

      DTYPE delta[3] = {dx, dy, dz};

      DTYPE mv[3] = {0.0, 0.0, 0.0};

      matvec3x3x3(sks, delta, mv);

      grad[1] += -2.0 * (a1 * a2) / (a1 + a2) * vij * mv[0];
      grad[2] += -2.0 * (a1 * a2) / (a1 + a2) * vij * mv[1];
      grad[3] += -2.0 * (a1 * a2) / (a1 + a2) * vij * mv[2];

      // x,y,z
      grad[4] += -2.0 * (a1 * a2) / (a1 + a2) * vij * delta[0];
      grad[5] += -2.0 * (a1 * a2) / (a1 + a2) * vij * delta[1];
      grad[6] += -2.0 * (a1 * a2) / (a1 + a2) * vij * delta[2];
    }
  }

  return grad;
}

////////////////////////////////////////////////////////////////////////////////
/// Optimization functions
////////////////////////////////////////////////////////////////////////////////

void adagrad_step(std::array<DTYPE, 4> &q, std::array<DTYPE, 3> &t,
                  std::array<DTYPE, 7> g, std::array<DTYPE, 7> &cache,
                  DTYPE lr_q, DTYPE lr_t) {
  cache = cache + g * g;

  q[0] -= lr_q * g[0] / (sqrt(cache[0]) + EPSILON);
  q[1] -= lr_q * g[1] / (sqrt(cache[1]) + EPSILON);
  q[2] -= lr_q * g[2] / (sqrt(cache[2]) + EPSILON);
  q[3] -= lr_q * g[3] / (sqrt(cache[3]) + EPSILON);

  t[0] -= lr_t * g[4] / (sqrt(cache[4]) + EPSILON);
  t[1] -= lr_t * g[5] / (sqrt(cache[5]) + EPSILON);
  t[2] -= lr_t * g[6] / (sqrt(cache[6]) + EPSILON);
}

std::array<DTYPE, 7> get_gradient_color(const DTYPE *molA, int NmolA,
                                        const int *molA_type, const DTYPE *molB,
                                        int NmolB, const int *molB_type,
                                        const DTYPE *rmat, const DTYPE *pmat,
                                        int N_features) {
  std::array<DTYPE, 7> grad = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  for (int i = 0; i < NmolA; i++) {
    int ta = molA_type[i];
    if (ta == 0) break;

    for (int j = 0; j < NmolB; j++) {
      int tb = molB_type[j];
      if (tb == 0) break;

      DTYPE dx = molA[i * D] - molB[j * D];
      DTYPE dy = molA[i * D + 1] - molB[j * D + 1];
      DTYPE dz = molA[i * D + 2] - molB[j * D + 2];

      DTYPE d2 = dx * dx + dy * dy + dz * dz;

      auto a = rmat[ta * N_features + tb];
      auto p = pmat[ta * N_features + tb];

      // see volume_color()
      if (a == 0) continue;

      DTYPE wa = molA[i * D + 3];
      DTYPE wb = molB[j * D + 3];

      DTYPE kij = exp(-a * a * d2 / (a + a)) * wa * wb;

      double constant = pow(PI / (2 * a), 1.5);

      DTYPE vij = p * p * kij * constant;

      DTYPE x = molB[j * D];
      DTYPE y = molB[j * D + 1];
      DTYPE z = molB[j * D + 2];

      DTYPE sks[3][3] = {
          {0, -2 * z, 2 * y}, {2 * z, 0, -2 * x}, {-2 * y, 2 * x, 0}};

      DTYPE delta[3] = {dx, dy, dz};

      DTYPE mv[3] = {0.0, 0.0, 0.0};

      matvec3x3x3(sks, delta, mv);

      grad[1] += -2.0 * (a * a) / (a + a) * vij * mv[0];
      grad[2] += -2.0 * (a * a) / (a + a) * vij * mv[1];
      grad[3] += -2.0 * (a * a) / (a + a) * vij * mv[2];

      // x,y,z
      grad[4] += -2.0 * (a * a) / (a + a) * vij * delta[0];
      grad[5] += -2.0 * (a * a) / (a + a) * vij * delta[1];
      grad[6] += -2.0 * (a * a) / (a + a) * vij * delta[2];
    }
  }
  return grad;
}

void single_conformer_optimiser(
    const DTYPE *molA, DTYPE *molAT, const int *molA_type, int NmolA,
    int NmolA_real, int NmolA_color, DTYPE self_overlap_A,
    DTYPE self_overlap_A_color, const DTYPE *molB, DTYPE *molBT,
    const int *molB_type, int NmolB, int NmolB_real, int NmolB_color,
    const DTYPE *rmat, const DTYPE *pmat, int N_features, bool optim_color,
    DTYPE mixing_param, DTYPE lr_q, DTYPE lr_t, int nsteps, DTYPE *scores) {
  std::array<DTYPE, 4> q = {1.0, 0.0, 0.0, 0.0};  // initial q
  std::array<DTYPE, 3> t = {0.0, 0.0, 0.0};       // initial t

  DTYPE M[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

  auto self_overlap_B = volume(molBT, NmolB_real, molBT, NmolB_real);
  // For the color volume, offset the coords to start at the color features.
  // With colour scoring off there are no colour features: the range holds only
  // padding and rmat/pmat are a 1x1 placeholder, so indexing them is both out
  // of bounds and numerically undefined. The colour volume is zero.
  DTYPE self_overlap_B_color = 0.0;
  if (optim_color) {
    self_overlap_B_color =
        volume_color(&molBT[NmolB_real * D], NmolB_color, &molB_type[NmolB_real],
                     &molBT[NmolB_real * D], NmolB_color, &molB_type[NmolB_real],
                     rmat, pmat, N_features);
  }
  std::array<DTYPE, 7> cache = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  // optimization loop
  for (int m = 0; m < nsteps; m++) {
    // important: we must compute gradients in the reference frame of molB

    // 1. first we rotate molB
    quaternion_to_rotation_matrix(q, M);

    // transform molB, putting the result in molBT.
    transform(M, molB, molBT, NmolB);

    // 2. translate molA by -t, result into molAT.
    DTYPE mt[3] = {-t[0], -t[1], -t[2]};
    translate(mt, molA, molAT, NmolA);

    auto g = get_gradient(molAT, NmolA_real, molBT, NmolB_real);
    g = g / (self_overlap_A + self_overlap_B);  // normalize

    if (optim_color) {
      // For the color gradients, offset the coords and types to start at the
      // color features.
      auto g_c = get_gradient_color(
          &molAT[NmolA_real * D], NmolA_color, &molA_type[NmolA_real],
          &molBT[NmolB_real * D], NmolB_color, &molB_type[NmolB_real], rmat,
          pmat, N_features);
      // normalize and combine
      g_c = g_c / (self_overlap_A_color + self_overlap_B_color);
      auto g_combo = (1 - mixing_param) * g + mixing_param * g_c;
      adagrad_step(q, t, g_combo, cache, lr_q, lr_t);
    } else {
      adagrad_step(q, t, g, cache, lr_q, lr_t);
    }
    // normalize q so that it is a unit quaternion
    DTYPE magq = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] * q[3] * q[3]);

    q = q / magq;
  }
  // get volume with transformed B
  quaternion_to_rotation_matrix(q, M);
  transform(M, molB, molBT, NmolB);
  translate(t.data(), molBT, molBT, NmolB);

  auto vol = volume(molA, NmolA_real, molBT, NmolB_real);

  DTYPE vc = 0.0;
  if (optim_color) {
    vc = volume_color(&molA[NmolA_real * D], NmolA_color, &molA_type[NmolA_real],
                      &molBT[NmolB_real * D], NmolB_color, &molB_type[NmolB_real],
                      rmat, pmat, N_features);
  }

  // compute tanimoto scores
  auto ts = vol / (self_overlap_A + self_overlap_B - vol);

  DTYPE tc = 0.0;
  if (optim_color) {
    tc = vc / (self_overlap_A_color + self_overlap_B_color - vc);
  }

  // we use the mixing param to weight the tanimotos
  // this is the objective function we have optimized
  auto total = (ts * (1 - mixing_param) + tc * mixing_param);

  scores[0] = total;                 // combination tanimoto of shape and color
  scores[1] = ts;                    // shape tanimoto
  scores[2] = tc;                    // color tanimoto
  scores[3] = vol;                   // volume shape
  scores[4] = vc;                    // volumes color
  scores[5] = self_overlap_A;        // self i
  scores[6] = self_overlap_B;        // self j
  scores[7] = self_overlap_A_color;  // self i color
  scores[8] = self_overlap_B_color;  // self j color
  scores[9] = q[0];
  scores[10] = q[1];
  scores[11] = q[2];
  scores[12] = q[3];
  scores[13] = t[0];
  scores[14] = t[1];
  scores[15] = t[2];
  scores[16] = 0.0;
  scores[17] = 0.0;
  scores[18] = 0.0;
  scores[19] = 0.0;
}
