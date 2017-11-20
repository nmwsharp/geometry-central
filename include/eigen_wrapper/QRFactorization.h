#pragma once

#include <Eigen/SparseQR>

#include "utilities.h"
#include "quaternion.h"

namespace GC {
template <class T>
class SparseMatrix;

template <class T>
class QRFactorization {
 public:
  QRFactorization(SparseMatrix<T>& A);
  ~QRFactorization(void);

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

  size_t rank();
  // Returns the rank using the QR factorization

 protected:
  void build(void);
  void update(void);

  SparseMatrix<T>& A;
  Eigen::SparseQR<Eigen::SparseMatrix<T, Eigen::ColMajor>,
                  Eigen::COLAMDOrdering<int>> solver;

  bool validSymbolic;
  bool validNumeric;
};
}

#include "eigen_wrapper/QRFactorization.ipp"