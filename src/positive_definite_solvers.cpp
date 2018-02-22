#include "geometrycentral/linear_solvers.h"

#include "geometrycentral/linear_algebra_utilities.h"

#ifdef HAVE_SUITESPARSE
#include "geometrycentral/suitesparse_utilities.h"
#endif

using namespace Eigen;
using std::cout;
using std::endl;

namespace geometrycentral {

template <typename T>
PositiveDefiniteSolver<T>::~PositiveDefiniteSolver() {
#ifdef HAVE_SUITESPARSE
  if (cMat != nullptr) {
    cholmod_l_free_sparse(&cMat, context);
    cMat = nullptr;
  }
  if (factorization != nullptr) {
    cholmod_l_free_factor(&factorization, context);
  }
#endif
}

template <typename T>
void PositiveDefiniteSolver<T>::prepare() {

  size_t N = this->mat.rows();

// Check some sanity
#ifndef GC_NLINALG_DEBUG
  if ((size_t)this->mat.cols() != N) {
    throw std::logic_error("Matrix must be square");
  }
  checkFinite(this->mat);
  checkHermitian(this->mat);
#endif

  this->mat.makeCompressed();

  // Suitesparse version
#ifdef HAVE_SUITESPARSE

  // Convert suitesparse format
  if (cMat != nullptr) {
    cholmod_l_free_sparse(&cMat, context);
  }
  cMat = toCholmod(this->mat, context, SType::SYMMETRIC);

  // Factor
  factorization = cholmod_l_analyze(cMat, context);
  cholmod_l_factorize(cMat, factorization, context);

  // Eigen version
#else
  solver.compute(this->mat);
  if (solver.info() != Eigen::Success) {
    std::cerr << "Solver factorization error: " << solver.info() << std::endl;
    throw std::invalid_argument("Solver factorization failed");
  }
#endif
};

template <typename T>
Vector<T> PositiveDefiniteSolver<T>::solve(const Vector<T>& rhs) {
  Vector<T> out;
  solve(out, rhs);
  return out;
}

template <typename T>
void PositiveDefiniteSolver<T>::solve(Vector<T>& x, const Vector<T>& rhs) {

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

  // Solve
  cholmod_dense* outVec = cholmod_l_solve(CHOLMOD_A, factorization, inVec, context);

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
Vector<T> solvePositiveDefinite(const Eigen::SparseMatrix<T>& A, const Vector<T>& rhs) {
  PositiveDefiniteSolver<T> s(A);
  return s.solve(rhs);
}


// Explicit instantiations
template class PositiveDefiniteSolver<double>;
template class PositiveDefiniteSolver<float>;
template class PositiveDefiniteSolver<Complex>;

template Vector<float> solvePositiveDefinite<float>(const Eigen::SparseMatrix<float>& A, const Vector<float>& rhs);
template Vector<double> solvePositiveDefinite<double>(const Eigen::SparseMatrix<double>& A, const Vector<double>& rhs);
template Vector<Complex> solvePositiveDefinite<Complex>(const Eigen::SparseMatrix<Complex>& A, const Vector<Complex>& rhs);


} // namespace geometrycentral
