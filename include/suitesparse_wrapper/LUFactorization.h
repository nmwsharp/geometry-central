#pragma once

#include <cholmod.h>
#include <umfpack.h>

#include "utilities.h"

namespace GC {
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

  void* get(void);
  // returns pointer to UMFPACK numeric factorization, building
  // (or rebuilding) the symbolic and numeric parts as necessary

 protected:
  void build(void);
  void update(void);

  SparseMatrix<T>& A;
  void* symbolic;
  void* numeric;

  bool validSymbolic;
  bool validNumeric;
};

template <>
void LUFactorization<Complex>::update(void);
}

#include "suitesparse_wrapper/LUFactorization.ipp"