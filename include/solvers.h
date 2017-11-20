#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <utilities.h>

// === Solvers ===

// Note: Since we only expect these functions to be called on a small handful of
// types (double, complex, etc), they are
// explicitly instantiated in the .cpp file.

// Compute the eigenvector corresponding to the smallest eigenvalue for a
// symmetric positive definite matrix using
// inverse power iterations
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> smallestEigenvectorPositiveDefinite(
    Eigen::SparseMatrix<T, Eigen::ColMajor> &energyMatrix,
    Eigen::SparseMatrix<T> &massMatrix, unsigned int nIterations = 50);

// Solve sparse general Ax=b
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> solve(
    Eigen::SparseMatrix<T, Eigen::ColMajor> &A,
    Eigen::Matrix<T, Eigen::Dynamic, 1> &b);

// Solve sparse Ax=b where A is square
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> solveSquare(
    Eigen::SparseMatrix<T, Eigen::ColMajor> &A,
    Eigen::Matrix<T, Eigen::Dynamic, 1> &b);

// Solve sparse Ax=b where A is symmetric positive definite
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> solvePositiveDefinite(
    Eigen::SparseMatrix<T, Eigen::ColMajor> &A,
    Eigen::Matrix<T, Eigen::Dynamic, 1> &b);

// === Eigen solver variants
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> solve_EigenQR(
    Eigen::SparseMatrix<T, Eigen::ColMajor> &A,
    Eigen::Matrix<T, Eigen::Dynamic, 1> &b);

// === Suitesparse solver variants
// (some of which are just Eigen wrappers)
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> solve_SPQR(
    Eigen::SparseMatrix<T, Eigen::ColMajor> &A,
    Eigen::Matrix<T, Eigen::Dynamic, 1> &b);
