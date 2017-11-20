#include "sparse_matrix.h"

// SuiteSparse's spqr.hpp rudely disables debugging and does not reenable it.
// Circumvent that here.
#ifndef NDEBUG
#define NDEBUG_WASNT_DEFINED
#endif
#include <spqr.hpp>
#ifdef NDEBUG_WASNT_DEFINED
#undef NDEBUG
#endif

#include <iostream>

namespace GC {
extern LinearContext context;

template <class T>
QRFactorization<T>::QRFactorization(SparseMatrix<T>& A_)
    : A(A_), factor(NULL), validSymbolic(false), validNumeric(false) {}

template <class T>
QRFactorization<T>::~QRFactorization(void) {
  clear();
}

template <class T>
void QRFactorization<T>::clear(void) {
  if (factor) {
    SuiteSparseQR_free<typename QRFactorization<T>::AUTO_ENTRYTYPE>(&factor,
                                                                    context);
    factor = NULL;
  }

  validSymbolic = false;
  validNumeric = false;
}

template <class T>
void QRFactorization<T>::clearNumeric(void) {
  validNumeric = false;
}

template <class T>
SuiteSparseQR_factorization<typename QRFactorization<T>::AUTO_ENTRYTYPE>*
QRFactorization<T>::get(void) {
  if (!validSymbolic) {
    build();
  } else if (!validNumeric) {
    update();
  }

  return factor;
}

template <class T>
void QRFactorization<T>::build(void) {
  clear();

  SparseMatrix<T> AT = A.transpose();
  cholmod_sparse* Ac = AT.to_cholmod();

  //  ordering options:
  //       0 or 3: fixed
  //       1: natural (only look for singletons)
  //       2: colamd after finding singletons
  //       4: CHOLMOD fter finding singletons
  //       5: amd(A'*A) after finding singletons
  //       6: metis(A'*A) after finding singletons
  //       7: SuiteSparseQR default (selects COLAMD, AMD, or METIS)
  const int ordering = 7;
  const int allow_tol = true;

  START_TIMING(sym)
  factor = SuiteSparseQR_symbolic<typename QRFactorization<T>::AUTO_ENTRYTYPE>(
      ordering, allow_tol, Ac, context);
  std::cout << "[qr] symbolic: " << pretty_time(FINISH_TIMING(sym))
            << std::endl;

  validSymbolic = true;

  update();
}

template <class T>
void QRFactorization<T>::update(void) {
  SparseMatrix<T> AT = A.transpose();
  cholmod_sparse* Ac = AT.to_cholmod();

  // only accept singletons above tol.  If tol <= -2, then
  // use the default tolerance
  const double tol = -2.;

  START_TIMING(num)
  SuiteSparseQR_numeric<typename QRFactorization<T>::AUTO_ENTRYTYPE>(
      tol, Ac, factor, context);
  std::cout << "[qr] numeric: " << pretty_time(FINISH_TIMING(num)) << std::endl;

  validNumeric = true;
}
}
