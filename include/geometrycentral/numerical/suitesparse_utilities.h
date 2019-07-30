#pragma once

#include "geometrycentral/numerical/linear_algebra_utilities.h"

#include "Eigen/Sparse"

// Suitesparse includes, as needed
#ifdef GC_HAVE_SUITESPARSE
#include "geometrycentral/numerical/suitesparse_utilities.h"
#include <SuiteSparseQR.hpp>
#include <cholmod.h>
#else
#error "USING SUITESPARSE FEATURES, BUT DON'T HAVE SUITESPARSE"
#endif


namespace geometrycentral {


enum class SType { UNSYMMETRIC = 0, SYMMETRIC };

// Suitesparse context class
class CholmodContext {
public:
  // constructor
  CholmodContext();

  // destructor
  ~CholmodContext();

  // set mode for Cholesky factorization
  void setSimplicial();
  void setSupernodal();

  // set LL vs LDL mode
  void setLL();
  void setLDL();

  // allows CholmodContext to be treated as a cholmod_common*
  operator cholmod_common*();

  cholmod_common context;
};


  // Type helper. This type is 'double' if T == 'float', and T otherwise
template<typename T>
struct SOLVER_ENTRYTYPE {
  typedef typename std::conditional<std::is_same<T, float>::value, double, T>::type type;
};

// === Conversion functions
// Caller is responsible for deallocating

// Convert a sparse matrix. Always returns a "full" matrix (stype=0, ie not triangular)
template <typename T>
cholmod_sparse* toCholmod(Eigen::SparseMatrix<T, Eigen::ColMajor>& A, CholmodContext& context,
                          SType stype = SType::UNSYMMETRIC);

// Convert a vector
template <typename T>
cholmod_dense* toCholmod(const Eigen::Matrix<T, Eigen::Dynamic, 1>& v, CholmodContext& context);

// Convert a vector
template <typename T>
void toEigen(cholmod_dense* cVec, CholmodContext& context, Eigen::Matrix<T, Eigen::Dynamic, 1>& xOut);

} // namespace geometrycentral
