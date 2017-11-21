#include "sparse_matrix.h"

#include <algorithm>
#include <iostream>

namespace geometrycentral {
extern LinearContext context;

template <class T>
LDLFactorization<T>::LDLFactorization(SparseMatrix<T>& A_)
    : A(A_), factor(NULL), validSymbolic(false), validNumeric(false) {}

template <class T>
LDLFactorization<T>::~LDLFactorization(void) {
  clear();
}

template <class T>
void LDLFactorization<T>::clear(void) {
  if (factor) {
    cholmod_l_free_factor(&factor, context);
    factor = NULL;
  }

  validSymbolic = false;
  validNumeric = false;
}

template <class T>
void LDLFactorization<T>::clearNumeric(void) {
  validNumeric = false;
}

template <class T>
cholmod_factor* LDLFactorization<T>::get(void) {
  if (!validSymbolic) {
    build();
  } else if (!validNumeric) {
    update();
  }

  return factor;
}

template <class T>
void LDLFactorization<T>::build(void) {
  clear();

  context.setSimplicial();

  cholmod_sparse* Ac = A.to_cholmod();
  Ac->stype = 1;

  factor = cholmod_l_analyze(Ac, context);

  validSymbolic = true;

  update();
}

template <class T>
void LDLFactorization<T>::update(void) {
  cholmod_sparse* Ac = A.to_cholmod();
  Ac->stype = 1;

  cholmod_l_factorize(Ac, factor, context);

  validNumeric = true;
}

template <class T>
void LDLFactorization<T>::makePositive(void) {
  size_t n = factor->n;
  SuiteSparse_long* p =
      (SuiteSparse_long*)factor->p;  // [0..ncol], the column pointers
  SuiteSparse_long* r =
      (SuiteSparse_long*)factor->i;  // [0..nzmax-1], the row indices
  double* x = (double*)factor->x;    // [0..nzmax-1], the numerical values
  SuiteSparse_long* nz = (SuiteSparse_long*)factor->nz;

  for (size_t i = 0; i < n; i++) {
    for (size_t l = p[i]; l < p[i] + nz[i]; l++) {
      size_t j = r[l];

      if (i == j) {
        const double delta = 1e-10;
        x[l] = fabs(x[l]);
        x[l] = std::max(x[l], delta);
      }
    }
  }
}
}