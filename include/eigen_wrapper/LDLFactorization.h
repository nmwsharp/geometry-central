#pragma once

#include <Eigen/SparseCholesky>

#include "utilities.h"

namespace GC {
template <class T>
class SparseMatrix;

template <class T>
class LDLFactorization {
 public:
  LDLFactorization(SparseMatrix<T>& A);
  ~LDLFactorization(void);

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
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<T>> solver;

  bool validSymbolic;
  bool validNumeric;
};
}

#include "eigen_wrapper/LDLFactorization.ipp"