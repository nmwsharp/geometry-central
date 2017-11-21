#include "sparse_matrix.h"

#include <iostream>

namespace geometrycentral {

template <class T>
LUFactorization<T>::LUFactorization(SparseMatrix<T>& A_)
    : A(A_), validSymbolic(false), validNumeric(false) {}

template <class T>
LUFactorization<T>::~LUFactorization(void) {}

template <class T>
void LUFactorization<T>::clear(void) {
  validSymbolic = false;
  validNumeric = false;
}

template <class T>
void LUFactorization<T>::clearNumeric(void) {
  validNumeric = false;
}

template <class T>
void LUFactorization<T>::solve(DenseVector<T>& x, const DenseVector<T>& b) {
  if (!validSymbolic) {
    build();
  } else if (!validNumeric) {
    update();
  }

  Eigen::Matrix<T, Eigen::Dynamic, 1> result = solver.solve(b.toEigen());
  x = result;
}

template <class T>
void LUFactorization<T>::build(void) {
  clear();

  Eigen::SparseMatrix<T> Aeigen = A.toEigen();
  Aeigen.makeCompressed();

  solver.analyzePattern(Aeigen);
  if (solver.info() != Eigen::Success) {
    std::string message = solver.lastErrorMessage();
    throw std::runtime_error("Matrix symbolic factorization failed. Message: " +
                             message);
  }
  validSymbolic = true;

  solver.factorize(Aeigen);
  if (solver.info() != Eigen::Success) {
    std::string message = solver.lastErrorMessage();
    throw std::runtime_error("Matrix numeric factorization failed. Message: " +
                             message);
  }
  validNumeric = true;
}

template <class T>
void LUFactorization<T>::update(void) {
  Eigen::SparseMatrix<T> Aeigen = A.toEigen();
  Aeigen.makeCompressed();

  solver.factorize(Aeigen);
  if (solver.info() != Eigen::Success) {
    std::string message = solver.lastErrorMessage();
    throw std::runtime_error("Matrix numeric factorization failed. Message: " +
                             message);
  }
  validNumeric = true;
}
}
