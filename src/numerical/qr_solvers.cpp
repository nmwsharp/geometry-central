#include "geometrycentral/numerical/linear_solvers.h"

#include "geometrycentral/numerical/linear_algebra_utilities.h"

#ifdef GC_HAVE_SUITESPARSE
#include "geometrycentral/numerical/suitesparse_utilities.h"
#include <SuiteSparseQR.hpp>
#include <cholmod.h>
#endif

namespace geometrycentral {

template <typename T>
struct QRSolverInternals {
  // Implementation-specific quantities
#ifdef GC_HAVE_SUITESPARSE
  CholmodContext context;
  cholmod_sparse* cMat = nullptr;
  cholmod_sparse* cMatTrans = nullptr;
  SuiteSparseQR_factorization<typename SOLVER_ENTRYTYPE<T>::type>* factorization = nullptr;
  double zero_tolerance = -2; // (use default)
#else
  Eigen::SparseQR<SparseMatrix<T>, Eigen::COLAMDOrdering<int>> solver;
#endif
};

template <typename T>
Solver<T>::~Solver() {
#ifdef GC_HAVE_SUITESPARSE
  if (internals->cMat != nullptr) {
    cholmod_l_free_sparse(&internals->cMat, internals->context);
    internals->cMat = nullptr;
  }
  if (internals->cMatTrans != nullptr) {
    cholmod_l_free_sparse(&internals->cMatTrans, internals->context);
    internals->cMatTrans = nullptr;
  }
  if (internals->factorization != nullptr) {
    SuiteSparseQR_free(&(internals->factorization), internals->context);
  }
#endif
}

template <typename T>
Solver<T>::Solver(SparseMatrix<T>& mat) : LinearSolver<T>(mat), internals(new QRSolverInternals<T>()) {

  // Is the system underdetermined?
  if (this->nRows < this->nCols) {
    underdetermined = true;
  } else {
    underdetermined = false;
  }


// Check some sanity
#ifndef GC_NLINALG_DEBUG
  checkFinite(mat);
#endif

  mat.makeCompressed();

// Suitesparse version
#ifdef GC_HAVE_SUITESPARSE

  // Convert suitesparse format
  // Either use A or A^T, depending on whether the system underdetermined
  if (internals->cMat != nullptr) {
    cholmod_l_free_sparse(&(internals->cMat), internals->context);
  }
  internals->cMat = toCholmod(mat, internals->context);
  if (underdetermined) {
    if (internals->cMatTrans != nullptr) {
      cholmod_l_free_sparse(&(internals->cMatTrans), internals->context);
    }
    internals->cMatTrans = cholmod_l_transpose(internals->cMat, 2, internals->context);

    // free the non-transposed matrix, its not needed
    if (internals->cMat != nullptr) {
      cholmod_l_free_sparse(&(internals->cMat), internals->context);
    }
  }

  // Factor

  //  ordering options:
  //       0 or 3: fixed
  //       1: natural (only look for singletons)
  //       2: colamd after finding singletons
  //       4: CHOLMOD fter finding singletons
  //       5: amd(A'*A) after finding singletons
  //       6: metis(A'*A) after finding singletons
  //       7: SuiteSparseQR default (selects COLAMD, AMD, or METIS)
  const int ordering = 7;

  if (internals->factorization != nullptr) {
    SuiteSparseQR_free(&(internals->factorization), internals->context);
  }

  if (underdetermined) {
    internals->factorization = SuiteSparseQR_factorize<typename SOLVER_ENTRYTYPE<T>::type>(
        ordering, internals->zero_tolerance, internals->cMatTrans, internals->context);
  } else {
    internals->factorization = SuiteSparseQR_factorize<typename SOLVER_ENTRYTYPE<T>::type>(
        ordering, internals->zero_tolerance, internals->cMat, internals->context);
  }

  if (internals->factorization == nullptr) {
    throw std::logic_error("Factorization failed");
  }

// Eigen version
#else
  if (underdetermined) {
    throw std::logic_error("Eigen's sparse QR solver doesn't like underdetermined systems");
  }

  internals->solver.setPivotThreshold(0.);
  internals->solver.compute(mat);
  if (internals->solver.info() != Eigen::Success) {
    std::cerr << "Solver factorization error: " << internals->solver.info() << std::endl;
    throw std::invalid_argument("Solver factorization failed");
  }
#endif
};

template <typename T>
Vector<T> Solver<T>::solve(const Vector<T>& rhs) {
  Vector<T> out;
  solve(out, rhs);
  return out;
}

template <typename T>
void Solver<T>::solve(Vector<T>& x, const Vector<T>& rhs) {

  // Check some sanity
  if ((size_t)rhs.rows() != this->nRows) {
    throw std::logic_error("Vector is not the right length");
  }
#ifndef GC_NLINALG_DEBUG
  checkFinite(rhs);
#endif

// Suitesparse version
#ifdef GC_HAVE_SUITESPARSE

  // Convert input to suitesparse format
  cholmod_dense* inVec = toCholmod(rhs, internals->context);
  cholmod_dense* outVec;

  // Solve
  // outVec = SuiteSparseQR<typename Solver<T>::SOLVER_ENTRYTYPE>(cMat, inVec, internals->context);
  // Note that the solve strategy is different for underdetermined systems
  if (underdetermined) {

    // solve y = R^-T b
    cholmod_dense* y = SuiteSparseQR_solve<typename SOLVER_ENTRYTYPE<T>::type>(
        SPQR_RTX_EQUALS_B, internals->factorization, inVec, internals->context);

    // compute x = Q*y
    outVec = SuiteSparseQR_qmult<typename SOLVER_ENTRYTYPE<T>::type>(SPQR_QX, internals->factorization, y,
                                                                     internals->context);
    cholmod_l_free_dense(&y, internals->context);

  } else {

    // compute y = Q^T b
    cholmod_dense* y = SuiteSparseQR_qmult<typename SOLVER_ENTRYTYPE<T>::type>(SPQR_QTX, internals->factorization,
                                                                               inVec, internals->context);

    // solve x = R^-1 y
    // TODO what is this E doing here?
    outVec = SuiteSparseQR_solve<typename SOLVER_ENTRYTYPE<T>::type>(SPQR_RETX_EQUALS_B, internals->factorization, y,
                                                                     internals->context);

    cholmod_l_free_dense(&y, internals->context);
  }

  // Convert back
  toEigen(outVec, internals->context, x);

  // Free
  cholmod_l_free_dense(&outVec, internals->context);
  cholmod_l_free_dense(&inVec, internals->context);

// Eigen version
#else
  // Solve
  x = internals->solver.solve(rhs);
  if (internals->solver.info() != Eigen::Success) {
    std::cerr << "Solver error: " << internals->solver.info() << std::endl;
    // std::cerr << "Solver says: " << solver.lastErrorMessage() << std::endl;
    throw std::invalid_argument("Solve failed");
  }
#endif
}

template <typename T>
Vector<T> solve(SparseMatrix<T>& A, const Vector<T>& rhs) {
  Solver<T> s(A);
  return s.solve(rhs);
}


template <typename T>
size_t Solver<T>::rank() {
#ifdef GC_HAVE_SUITESPARSE
  return internals->factorization->rank;
#else
  return internals->solver.rank();
#endif
}

// Explicit instantiations
template class Solver<double>;
template class Solver<float>;
template class Solver<std::complex<double>>;

template Vector<float> solve(SparseMatrix<float>& A, const Vector<float>& rhs);
template Vector<double> solve(SparseMatrix<double>& A, const Vector<double>& rhs);
template Vector<std::complex<double>> solve(SparseMatrix<std::complex<double>>& A,
                                            const Vector<std::complex<double>>& rhs);

} // namespace geometrycentral
