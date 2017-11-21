#pragma once

#include <Eigen/SparseCholesky>

#include "utilities.h"

namespace geometrycentral {
template <class T>
class SparseMatrix;

template <class T>
class LLFactorization {
 public:
  LLFactorization(SparseMatrix<T>& A);
  ~LLFactorization(void);

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
  Eigen::SimplicialLLT<Eigen::SparseMatrix<T>> solver;

  bool validSymbolic;
  bool validNumeric;
};
}

#include "eigen_wrapper/LLFactorization.ipp"