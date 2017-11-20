#include "sparse_matrix.h"

#include <iostream>

namespace GC {

template <class T>
LDLFactorization<T>::LDLFactorization(SparseMatrix<T>& A_)
    : A(A_), validSymbolic(false), validNumeric(false) {}

template <class T>
LDLFactorization<T>::~LDLFactorization(void) {}

template <class T>
void LDLFactorization<T>::clear(void) {
  validSymbolic = false;
  validNumeric = false;
}

template <class T>
void LDLFactorization<T>::clearNumeric(void) {
  validNumeric = false;
}

template <class T>
void LDLFactorization<T>::solve(DenseVector<T>& x, const DenseVector<T>& b) {
  if (!validSymbolic) {
    build();
  } else if (!validNumeric) {
    update();
  }

  Eigen::Matrix<T, Eigen::Dynamic, 1> result = solver.solve(b.toEigen());
  x = result;
}

template <class T>
void LDLFactorization<T>::build(void) {
  clear();

  Eigen::SparseMatrix<T> Aeigen = A.toEigen();
  Aeigen.makeCompressed();

  solver.analyzePattern(Aeigen);
  if (solver.info() != Eigen::Success) {
    throw std::runtime_error("Matrix symbolic factorization failed.");
  }
  validSymbolic = true;

  solver.factorize(Aeigen);
  if (solver.info() != Eigen::Success) {
    throw std::runtime_error("Matrix numeric factorization failed.");
  }
  validNumeric = true;
}

template <class T>
void LDLFactorization<T>::update(void) {
  Eigen::SparseMatrix<T> Aeigen = A.toEigen();
  Aeigen.makeCompressed();

  solver.factorize(Aeigen);
  if (solver.info() != Eigen::Success) {
    throw std::runtime_error("Matrix numeric factorization failed.");
  }
  validNumeric = true;
}
}