// -----------------------------------------------------------------------------
// libDDG -- LinearContext.h
// -----------------------------------------------------------------------------
//
// LinearContext is the global solver context needed to interface with the
// SuiteSparse library.  It is essentially a wrapper around cholmod_common.  A
// single static instance of LinearContext is declared in LinearContext.cpp and
// is shared by all instances of DenseMatrix, SparseMatrix, and LinearSystem.
// In other words, you shouldn't have to instantiate LinearContext yourself
// unless you're doing something really fancy!
//

#pragma once

#ifdef HAVE_SUITESPARSE
#include <cholmod.h>

namespace GC {
class LinearContext {
 public:
  LinearContext(void);
  // constructor

  ~LinearContext(void);
  // destructor

  void setSimplicial(void);
  void setSupernodal(void);
  // set mode for Cholesky factorization

  operator cholmod_common*(void);
  // allows LinearContext to be treated as a cholmod_common*

 protected:
  cholmod_common context;
};
}
#endif
