#include "geometrycentral/sparse_matrix.h"

#include <iostream>

namespace geometrycentral {
template <class T>
LLFactorization<T>::LLFactorization(SparseMatrix<T>& A_)
    : A(A_), validSymbolic(false), validNumeric(false) {}

template <class T>
LLFactorization<T>::~LLFactorization(void) {}

template <class T>
void LLFactorization<T>::clear(void) {
  validSymbolic = false;
  validNumeric = false;
}

template <class T>
void LLFactorization<T>::clearNumeric(void) {
  validNumeric = false;
}

template <class T>
void LLFactorization<T>::solve(DenseVector<T>& x, const DenseVector<T>& b) {
  if (!validSymbolic) {
    build();
  } else if (!validNumeric) {
    update();
  }

  Eigen::Matrix<T, Eigen::Dynamic, 1> result = solver.solve(b.toEigen());
  x = result;
}

template <class T>
void LLFactorization<T>::build(void) {
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
void LLFactorization<T>::update(void) {
  Eigen::SparseMatrix<T> Aeigen = A.toEigen();
  Aeigen.makeCompressed();

  solver.factorize(Aeigen);
  if (solver.info() != Eigen::Success) {
    throw std::runtime_error("Matrix numeric factorization failed.");
  }

  validNumeric = true;
}
}