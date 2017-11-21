#pragma once

#include "eigen_wrapper/LLFactorization.h"
#include "eigen_wrapper/LDLFactorization.h"
#include "eigen_wrapper/LUFactorization.h"
#include "eigen_wrapper/QRFactorization.h"

namespace geometrycentral {
template <class T>
class SparseMatrix;

template <class T>
class SparseFactorization {
 public:
  SparseFactorization(SparseMatrix<T>& A);

  void clear(void);
  void clearNumeric(void);

  LLFactorization<T> ll;
  LDLFactorization<T> ldl;
  LUFactorization<T> lu;
  QRFactorization<T> qr;
};
}

#include "sparse_matrix.h"
#include "eigen_wrapper/sparse_factorization.ipp"