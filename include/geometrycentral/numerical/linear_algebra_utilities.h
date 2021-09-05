#pragma once

#include "geometrycentral/numerical/linear_algebra_types.h"
#include "geometrycentral/utilities/utilities.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/StdVector>

#include <complex>
#include <fstream> // ofsteam, ifstream
#include <iomanip> // setprecision
#include <iostream>

// === Various helper functions and sanity checks which are useful for linear algebra code


namespace geometrycentral {

// ==== Simple utilities

template <typename T>
SparseMatrix<T> identityMatrix(size_t N);

// Shift the the diagonal of a matrix by a constant offset
template <typename T>
void shiftDiagonal(SparseMatrix<T>& m, T shiftAmount = 1e-4);

template <typename T>
SparseMatrix<T> verticalStack(const std::vector<SparseMatrix<T>, Eigen::aligned_allocator<SparseMatrix<T>>>& mats);

template <typename T>
SparseMatrix<T> horizontalStack(const std::vector<SparseMatrix<T>, Eigen::aligned_allocator<SparseMatrix<T>>>& mats);

// Blow up an NxM complex system to a 2N x2M real system.
SparseMatrix<double> complexToReal(const SparseMatrix<std::complex<double>>& m);
Vector<double> complexToReal(const Vector<std::complex<double>>& v);
Vector<std::complex<double>> realToComplex(const Vector<double>& v);

template <typename T>
std::vector<std::vector<T>> unpackMatrixToStdVector(const DenseMatrix<T>& mat);

// ==== Sanity checks


// Verify that a matrix has finite entries, error out if not
template <typename T>
void checkFinite(const SparseMatrix<T>& m);

template <typename T, int R, int C>
void checkFinite(const Eigen::Matrix<T, R, C>& m);

// Verify that a sparse matrix is symmetric, error out if not.
template <typename T>
void checkSymmetric(const SparseMatrix<T>& m, double absoluteEPS = -1.);

// Verify that a sparse matrix is hermitian, error out if not.  For real matrices, coincides with checkSymmetric
template <typename T>
void checkHermitian(const SparseMatrix<T>& m, double absoluteEPS = -1.);


// ==== Permutations and blocking

// Build a permutation matrix
// template <typename T>
// SparseMatrix<T> buildPermutationMatrix(const Vector<size_t>& p);

// Apply a permutation to the rows and columns of a matrix
// template <typename T>
// void permute(const SparseMatrix<T>& m, const Vector<size_t>& p);


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
  SparseMatrix<T> AA;
  SparseMatrix<T> AB;
  SparseMatrix<T> BA;
  SparseMatrix<T> BB;
};
template <typename T>
BlockDecompositionResult<T> blockDecomposeSquare(const SparseMatrix<T>& m, const Vector<bool>& Aset,
                                                 bool buildBuildBside = true);

// Apply a decomposition to a vector
template <typename T>
void decomposeVector(BlockDecompositionResult<T>& decomp, const Vector<T>& vec, Vector<T>& vecAOut, Vector<T>& vecBOut);
template <typename T>
Vector<T> reassembleVector(BlockDecompositionResult<T>& decomp, const Vector<T>& vecA, const Vector<T>& vecB);

// === IO

// WARNING: this follows matlab convention and thus is 1-indexed
template <typename T>
void saveSparseMatrix(std::string filename, SparseMatrix<T>& matrix);
template <typename T>
void saveSparseMatrix(std::ostream& out, SparseMatrix<T>& matrix);

template <typename T>
void saveDenseMatrix(std::string filename, DenseMatrix<T>& matrix);
template <typename T>
void saveDenseMatrix(std::ostream& out, DenseMatrix<T>& matrix);

template <typename T>
SparseMatrix<T> loadSparseMatrix(std::string filename);
template <typename T>
SparseMatrix<T> loadSparseMatrix(std::istream& in);

template <typename T>
DenseMatrix<T> loadDenseMatrix(std::string filename);
template <typename T>
DenseMatrix<T> loadDenseMatrix(std::istream& in);


#include "geometrycentral/numerical/linear_algebra_utilities.ipp"

} // namespace geometrycentral
