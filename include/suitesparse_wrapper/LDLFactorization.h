#pragma once

#include <cholmod.h>
#include "utilities.h"

namespace geometrycentral {
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

  cholmod_factor* get(void);
  // returns pointer to CHOLMOD factorization, building (or
  // rebuilding) the symbolic and numeric parts as necessary

  void makePositive(void);
  // converts to a factorization of a nearby positive-definite
  // matrix by flipping the sign of negative diagonal entries

 protected:
  void build(void);
  void update(void);

  SparseMatrix<T>& A;
  cholmod_factor* factor;

  bool validSymbolic;
  bool validNumeric;
};
}

#include "suitesparse_wrapper/LDLFactorization.ipp"