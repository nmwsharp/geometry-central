#include "geometrycentral/linear_solvers.h"

#include "geometrycentral/linear_algebra_utilities.h"

#ifdef HAVE_SUITESPARSE
#include "geometrycentral/suitesparse_utilities.h"
#include <SuiteSparseQR.hpp>
#include <cholmod.h>
#endif


using namespace Eigen;

namespace geometrycentral {

template <typename T>
Solver<T>::~Solver() {
#ifdef HAVE_SUITESPARSE
  if (cMat != nullptr) {
    cholmod_l_free_sparse(&cMat, context);
    cMat = nullptr;
  }
  if (cMatTrans != nullptr) {
    cholmod_l_free_sparse(&cMatTrans, context);
    cMatTrans = nullptr;
  }
  if (factorization != nullptr) {
    SuiteSparseQR_free(&factorization, context);
  }
#endif
}

template <typename T>
void Solver<T>::prepare() {

  size_t Nrows = this->mat.rows();
  size_t Ncols = this->mat.cols();

// Check some sanity
#ifndef GC_NLINALG_DEBUG
  checkFinite(this->mat);
#endif

  this->mat.makeCompressed();

// Suitesparse version
#ifdef HAVE_SUITESPARSE

  // Is the system underdetermined?
  if (Nrows < Ncols) {
    underdetermined = true;
     throw std::logic_error("is not well tested, be careful");
  } else {
    underdetermined = false;
  }

  // Convert suitesparse format
  // Either use A or A^T, depending on whether the system underdetermined
  if (cMat != nullptr) {
    cholmod_l_free_sparse(&cMat, context);
  }
  cMat = toCholmod(this->mat, context);
  if (underdetermined) {
    if (cMatTrans != nullptr) {
      cholmod_l_free_sparse(&cMatTrans, context);
    }
    cMatTrans = cholmod_l_transpose(cMat, 2, context);
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

  if (factorization != nullptr) {
    SuiteSparseQR_free(&factorization, context);
  }

  if (underdetermined) {
    factorization =
        SuiteSparseQR_factorize<typename Solver<T>::SOLVER_ENTRYTYPE>(ordering, zero_tolerance, cMatTrans, context);
  } else {
    factorization =
        SuiteSparseQR_factorize<typename Solver<T>::SOLVER_ENTRYTYPE>(ordering, zero_tolerance, cMat, context);
  }

  if (factorization == nullptr) {
    throw std::logic_error("Factorization failed");
  }

// Eigen version
#else
  solver.compute(this->mat);
  if (solver.info() != Eigen::Success) {
    std::cerr << "Solver factorization error: " << solver.info() << std::endl;
    throw std::invalid_argument("Solver factorization failed");
  }
  std::cout << "Eigen done factoring" << std::endl;
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

  size_t N = this->mat.rows();

// Check some sanity
#ifndef GC_NLINALG_DEBUG
  if ((size_t)rhs.rows() != N) {
    throw std::logic_error("Vector is not the right length");
  }
  checkFinite(rhs);
#endif

// Suitesparse version
#ifdef HAVE_SUITESPARSE

  // Convert input to suitesparse format
  cholmod_dense* inVec = toCholmod(rhs, context);
  cholmod_dense* outVec;

  // Solve
  // outVec = SuiteSparseQR<typename Solver<T>::SOLVER_ENTRYTYPE>(cMat, inVec, context);
  // Note that the solve strategy is different for underdetermined systems
  if (underdetermined) {

    // solve y = R^-T b
    cholmod_dense* y =
        SuiteSparseQR_solve<typename Solver<T>::SOLVER_ENTRYTYPE>(SPQR_RTX_EQUALS_B, factorization, inVec, context);

    // compute x = Q*y
    outVec = SuiteSparseQR_qmult<typename Solver<T>::SOLVER_ENTRYTYPE>(SPQR_QX, factorization, y, context);
    cholmod_l_free_dense(&y, context);

  } else {

    // compute y = Q^T b
    cholmod_dense* y =
        SuiteSparseQR_qmult<typename Solver<T>::SOLVER_ENTRYTYPE>(SPQR_QTX, factorization, inVec, context);

    // solve x = R^-1 y
    // TODO what is this E doing here?
    outVec = SuiteSparseQR_solve<typename Solver<T>::SOLVER_ENTRYTYPE>(SPQR_RETX_EQUALS_B, factorization, y, context);

    cholmod_l_free_dense(&y, context);
  }

  // Convert back
  toEigen(outVec, context, x);

  // Free
  cholmod_l_free_dense(&outVec, context);
  cholmod_l_free_dense(&inVec, context);

// Eigen version
#else
  // Solve
  x = solver.solve(rhs);
  if (solver.info() != Eigen::Success) {
    std::cerr << "Solver error: " << solver.info() << std::endl;
    // std::cerr << "Solver says: " << solver.lastErrorMessage() << std::endl;
    throw std::invalid_argument("Solve failed");
  }
#endif

// Compute residual to spot bad solves
#ifndef GC_NLINALG_DEBUG
  Matrix<T, Dynamic, 1> residual = this->mat * x - rhs;
  double residualNorm = residual.norm();
  double relativeResidualNorm = residualNorm / rhs.norm();
  std::cout << "  -- Residual norm: " << residualNorm << "   relative residual norm: " << relativeResidualNorm
            << std::endl;
#endif
}

template <typename T>
Vector<T> solve(const Eigen::SparseMatrix<T>& A, const Vector<T>& rhs) {
  Solver<T> s(A);
  return s.solve(rhs);
}


template <typename T>
size_t Solver<T>::rank() {
#ifdef HAVE_SUITESPARSE
  return factorization->rank;
#else
  return solver.rank();
#endif
}

// Explicit instantiations
template class Solver<double>;
template class Solver<float>;
template class Solver<Complex>;

template Vector<float> solve(const Eigen::SparseMatrix<float>& A, const Vector<float>& rhs);
template Vector<double> solve(const Eigen::SparseMatrix<double>& A, const Vector<double>& rhs);
template Vector<Complex> solve(const Eigen::SparseMatrix<Complex>& A, const Vector<Complex>& rhs);

} // namespace geometrycentral
