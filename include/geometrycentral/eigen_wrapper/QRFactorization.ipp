#include "geometrycentral/sparse_matrix.h"

#include <iostream>

namespace geometrycentral {

template <class T>
QRFactorization<T>::QRFactorization(SparseMatrix<T>& A_)
    : A(A_), validSymbolic(false), validNumeric(false) {}

template <class T>
QRFactorization<T>::~QRFactorization(void) {}

template <class T>
void QRFactorization<T>::clear(void) {
  validSymbolic = false;
  validNumeric = false;
}

template <class T>
void QRFactorization<T>::clearNumeric(void) {
  validNumeric = false;
}

template <class T>
void QRFactorization<T>::solve(DenseVector<T>& x, const DenseVector<T>& b) {
  if (!validSymbolic) {
    build();
  } else if (!validNumeric) {
    update();
  }

  Eigen::Matrix<T, Eigen::Dynamic, 1> result = solver.solve(b.toEigen());
  x = result;
}

template <class T>
size_t QRFactorization<T>::rank() {
  if (!validSymbolic) {
    build();
  } else if (!validNumeric) {
    update();
  }

  return solver.rank();
}

template <class T>
void QRFactorization<T>::build(void) {
  clear();

  Eigen::SparseMatrix<T> Aeigen = A.toEigen();
  Aeigen.makeCompressed();

  // NOTE: Eigen library always seems to throw errors when using incrementally
  // like this...
  //       Might be a library issue or a me issue, but either way updates are
  //       disabled with
  //       Eigen SparseQR for now.
  // solver.analyzePattern(Aeigen);
  // if(solver.info() != Eigen::Success) {
  //    std::string message = solver.lastErrorMessage();
  //    throw std::runtime_error("Matrix symbolic factorization failed. Message:
  //    " + message);
  // }
  // validSymbolic = true;

  // solver.factorize(Aeigen);
  // if(solver.info() != Eigen::Success) {
  //    std::string message = solver.lastErrorMessage();
  //    throw std::runtime_error("Matrix numeric factorization failed. Message:
  //    " + message);
  // }
  // validNumeric = true;

  solver.compute(Aeigen);
  if (solver.info() != Eigen::Success) {
    std::string message = solver.lastErrorMessage();
    throw std::runtime_error("Matrix factorization failed. Message: " +
                             message);
  }

  validSymbolic = true;
  validNumeric = true;
}

template <class T>
void QRFactorization<T>::update(void) {
  build();

  // Eigen::SparseMatrix<T> Aeigen = A.toEigen();
  // Aeigen.makeCompressed();

  // solver.factorize(Aeigen);
  // if(solver.info() != Eigen::Success) {
  //    std::string message = solver.lastErrorMessage();
  //    throw std::runtime_error("Matrix numeric factorization failed. Message:
  //    " + message);
  // }
  // validNumeric = true;
}
}
