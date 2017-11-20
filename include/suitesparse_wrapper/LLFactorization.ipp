#include "sparse_matrix.h"

#include <iostream>

namespace GC {
extern LinearContext context;

template <class T>
LLFactorization<T>::LLFactorization(SparseMatrix<T>& A_)
    : A(A_), factor(NULL), validSymbolic(false), validNumeric(false) {}

template <class T>
LLFactorization<T>::~LLFactorization(void) {
  clear();
}

template <class T>
void LLFactorization<T>::clear(void) {
  if (factor) {
    cholmod_l_free_factor(&factor, context);
    factor = NULL;
  }

  validSymbolic = false;
  validNumeric = false;
}

template <class T>
void LLFactorization<T>::clearNumeric(void) {
  validNumeric = false;
}

template <class T>
cholmod_factor* LLFactorization<T>::get(void) {
  if (!validSymbolic) {
    build();
  } else if (!validNumeric) {
    update();
  }

  return factor;
}

template <class T>
void LLFactorization<T>::build(void) {
  clear();

  context.setSupernodal();

  cholmod_sparse* Ac = A.to_cholmod();
  Ac->stype = 1;

  factor = cholmod_l_analyze(Ac, context);

  validSymbolic = true;

  update();
}

template <class T>
void LLFactorization<T>::update(void) {
  cholmod_sparse* Ac = A.to_cholmod();
  Ac->stype = 1;

  cholmod_l_factorize(Ac, factor, context);

  validNumeric = true;
}
}