#pragma once

#include <cholmod.h>

#include <SuiteSparseQR.hpp>

#include "quaternion.h"
#include "utilities.h"

namespace geometrycentral {
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

  // Type system craziness that maps to the proper entry type.
  // For every type T except quaternion, the entry type should be T. For
  // T=Quaternion
  // it should be double. This implements that mapping.
  typedef typename std::conditional<std::is_same<T, Quaternion>::value, double,
                                    T>::type AUTO_ENTRYTYPE;

  SuiteSparseQR_factorization<typename QRFactorization<T>::AUTO_ENTRYTYPE>* get(
      void);
  // returns pointer to SuiteSparseQR factorization, building
  // (or rebuilding) the symbolic and numeric parts as necessary

 protected:
  void build(void);
  void update(void);

  SparseMatrix<T>& A;
  SuiteSparseQR_factorization<typename QRFactorization<T>::AUTO_ENTRYTYPE>*
      factor;

  bool validSymbolic;
  bool validNumeric;
};
}

#include "suitesparse_wrapper/QRFactorization.ipp"