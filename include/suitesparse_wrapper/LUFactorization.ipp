#include "sparse_matrix.h"

#include <iostream>

namespace GC {
extern LinearContext context;

template <class T>
LUFactorization<T>::LUFactorization(SparseMatrix<T>& A_)
    : A(A_),
      symbolic(NULL),
      numeric(NULL),
      validSymbolic(false),
      validNumeric(false) {}

template <class T>
LUFactorization<T>::~LUFactorization(void) {
  clear();
}

template <class T>
void LUFactorization<T>::clear(void) {
  if (symbolic) {
    umfpack_dl_free_symbolic(&symbolic);
    symbolic = NULL;
  }

  if (numeric) {
    umfpack_dl_free_numeric(&numeric);
    numeric = NULL;
  }
}

template <class T>
void LUFactorization<T>::clearNumeric(void) {
  validNumeric = false;
}

template <class T>
void* LUFactorization<T>::get(void) {
  if (!validSymbolic) {
    build();
  } else if (!validNumeric) {
    update();
  }

  return numeric;
}

template <class T>
void LUFactorization<T>::build(void) {
  clear();

  cholmod_sparse* Ac = A.to_cholmod();
  int n = Ac->nrow;
  SuiteSparse_long* Ap = (SuiteSparse_long*)Ac->p;
  SuiteSparse_long* Ai = (SuiteSparse_long*)Ac->i;
  double* Ax = (double*)Ac->x;

  START_TIMING(sym)
  umfpack_dl_symbolic(n, n, Ap, Ai, Ax, &symbolic, NULL, NULL);
  std::cout << "[lu] symbolic: " << pretty_time(FINISH_TIMING(sym))
            << std::endl;

  validSymbolic = true;

  update();
}

template <class T>
void LUFactorization<T>::update(void) {
  cholmod_sparse* Ac = A.to_cholmod();
  SuiteSparse_long* Ap = (SuiteSparse_long*)Ac->p;
  SuiteSparse_long* Ai = (SuiteSparse_long*)Ac->i;
  double* Ax = (double*)Ac->x;

  if (numeric) {
    umfpack_dl_free_numeric(&numeric);
    numeric = NULL;
  }

  START_TIMING(num)
  umfpack_dl_numeric(Ap, Ai, Ax, symbolic, &numeric, NULL, NULL);
  std::cout << "[lu] numeric: " << pretty_time(FINISH_TIMING(num)) << std::endl;

  validNumeric = true;
}
}
