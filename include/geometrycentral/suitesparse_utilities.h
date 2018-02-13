#pragma once

#include "Eigen/Sparse"

// Suitesparse includes, as needed
#ifdef HAVE_SUITESPARSE
#include <cholmod.h>
#else
#error "USING SUITESPARSE FEATURES, BUT DON'T HAVE SUITESPARSE"
#endif

namespace geometrycentral {


enum class SType { UNSYMMETRIC = 0, SYMMETRIC};

// Suitesparse context class
class CholmodContext {
public:
  // constructor
  CholmodContext(void);

  // destructor
  ~CholmodContext(void);

  // set mode for Cholesky factorization
  void setSimplicial(void);
  void setSupernodal(void);

  // allows CholmodContext to be treated as a cholmod_common*
  operator cholmod_common*(void);

protected:
  cholmod_common context;
};

// === Conversion functions
// Caller is responsible for deallocating

// Convert a sparse matrix. Always returns a "full" matrix (stype=0, ie not triangular)
template <typename T>
cholmod_sparse* toCholmod(Eigen::SparseMatrix<T, Eigen::ColMajor>& A, CholmodContext& context, SType stype=SType::UNSYMMETRIC);

// Convert a vector
template <typename T>
cholmod_dense* toCholmod(const Eigen::Matrix<T, Eigen::Dynamic, 1>& v, CholmodContext& context);

// Convert a vector
template <typename T>
void toEigen(cholmod_dense* cVec, CholmodContext& context, Eigen::Matrix<T, Eigen::Dynamic, 1>& xOut);

} // namespace geometrycentral
