#pragma once

#include "geometrycentral/utilities.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <iostream>

// === Various helper functions and sanity checks which are useful for linear algebra code

// Convenience typedef
namespace {
template <typename T> using Vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;
}

namespace geometrycentral {


// ==== Simple utilities

template <typename T>
Eigen::SparseMatrix<T> identityMatrix(size_t N);

// Shift the the diagonal of a matrix by a constant offset
template <typename T>
void shiftDiagonal(Eigen::SparseMatrix<T>& m, T shiftAmount = 1e-4);


// ==== Sanity checks


// Verify that a matrix has finite entries, error out if not. Does nothing if NDEBUG is defined.
template <typename T>
void checkFinite(const Eigen::SparseMatrix<T>& m);

template <typename T, int R, int C>
void checkFinite(const Eigen::Matrix<T, R, C>& m);

// Verify that a sparse matrix is symmetric (hermitian), error out if not. Does nothing if NDEBUG is defined.
template <typename T>
void checkHermitian(const Eigen::SparseMatrix<T>& m);


// ==== Permutations and blocking

// Build a permutation matrix
// template <typename T>
// Eigen::SparseMatrix<T> buildPermutationMatrix(const Vector<size_t>& p);

// Apply a permutation to the rows and columns of a matrix
// template <typename T>
// void permute(const Eigen::SparseMatrix<T>& m, const Vector<size_t>& p);


// Block-decompose a square sparse matrix with interleaved index sets A and B
template <typename T>
struct BlockDecompositionResult {
  // The index of each element of A (resp. B) in the original system
  Vector<size_t> origIndsA;
  Vector<size_t> origIndsB;
 
  // Index of each orignal element in the new system  (either in the A system or B)
  Vector<size_t> newInds;
  Vector<bool> isA;

  // The four "block" matrices
  Eigen::SparseMatrix<T> AA;
  Eigen::SparseMatrix<T> AB;
  Eigen::SparseMatrix<T> BA;
  Eigen::SparseMatrix<T> BB;
};
template <typename T>
BlockDecompositionResult<T> blockDecomposeSquare(const Eigen::SparseMatrix<T>& m, const Vector<bool>& Aset, bool buildBuildBside=true);

// Apply a decomposition to a vector
template <typename T>
void decomposeVector(BlockDecompositionResult<T>& decomp, const Vector<T>& vec, Vector<T>& vecAOut, Vector<T>& vecBOut);

#include "geometrycentral/linear_algebra_utilities.ipp"

} // namespace geometrycentral
