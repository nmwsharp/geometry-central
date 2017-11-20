#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

// === Various helper functions and sanity checks which are useful for linear
// algebra code

// Verify that a matrix has finite entries, error out if not. Does nothing if
// NDEBUG is defined.
template <typename T>
void checkFinite(const Eigen::SparseMatrix<T>& m);

template <typename T, int R, int C>
void checkFinite(const Eigen::Matrix<T, R, C>& m);

// Verify that a sparse matrix is symmetric (hermitian), error out if not. Does
// nothing if NDEBUG is defined.
void checkSymmetric(const Eigen::SparseMatrix<double>& m);
void checkSymmetric(const Eigen::SparseMatrix<std::complex<double>>& m);

#include <linear_algebra_utilities.ipp>
