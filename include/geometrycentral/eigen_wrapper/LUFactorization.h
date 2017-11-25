#pragma once

#include <Eigen/SparseLU>

#include "geometrycentral/utilities.h"

namespace geometrycentral {
template <class T>
class SparseMatrix;

template <class T>
class LUFactorization {
 public:
  LUFactorization(SparseMatrix<T>& A);
  ~LUFactorization(void);

  void clear(void);
  // clears both the symbolic and numeric factorization -- should
  // be called following any change to the nonzero pattern of the
  // corresponding matrix

  void clearNumeric(void);
  // clears only the numeric factorization -- should be called
  // following any change to the values of nonzero entries in
  // the corresponding matrix

  void solve(DenseVector<T>& x, const DenseVector<T>& b);
  // solve using this factorization, building the factorization
  // first if needed

 protected:
  void build(void);
  void update(void);

  SparseMatrix<T>& A;
  Eigen::SparseLU<Eigen::SparseMatrix<T>> solver;

  bool validSymbolic;
  bool validNumeric;
};

template <>
void LUFactorization<Complex>::update(void);
}

#include "geometrycentral/eigen_wrapper/LUFactorization.ipp"