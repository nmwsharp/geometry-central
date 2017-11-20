#pragma once

#include "suitesparse_wrapper/LLFactorization.h"
#include "suitesparse_wrapper/LDLFactorization.h"
#include "suitesparse_wrapper/LUFactorization.h"
#include "suitesparse_wrapper/QRFactorization.h"

namespace GC {
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
#include "suitesparse_wrapper/sparse_factorization.ipp"