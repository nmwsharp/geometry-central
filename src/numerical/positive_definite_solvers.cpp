#include "geometrycentral/numerical/linear_solvers.h"

#include "geometrycentral/numerical/linear_algebra_utilities.h"

#ifdef GC_HAVE_SUITESPARSE
#include "geometrycentral/numerical/suitesparse_utilities.h"
#endif

using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;

namespace geometrycentral {

template <typename T>
struct PSDSolverInternals {
#ifdef GC_HAVE_SUITESPARSE
  CholmodContext context;
  cholmod_sparse* cMat = nullptr;
  cholmod_factor* factorization = nullptr;
#else
  Eigen::SimplicialLDLT<SparseMatrix<T>> solver;
#endif
};

template <typename T>
PositiveDefiniteSolver<T>::~PositiveDefiniteSolver() {
#ifdef GC_HAVE_SUITESPARSE
  if (internals->cMat != nullptr) {
    cholmod_l_free_sparse(&internals->cMat, internals->context);
    internals->cMat = nullptr;
  }
  if (internals->factorization != nullptr) {
    cholmod_l_free_factor(&internals->factorization, internals->context);
  }
#endif
}

template <typename T>
PositiveDefiniteSolver<T>::PositiveDefiniteSolver(SparseMatrix<T>& mat)
    : LinearSolver<T>(mat), internals(new PSDSolverInternals<T>()) {


  // Check some sanity
  if (this->nRows != this->nCols) {
    throw std::logic_error("Matrix must be square");
  }
  size_t N = this->nRows;
#ifndef GC_NLINALG_DEBUG
  checkFinite(mat);
  checkHermitian(mat);
#endif

  mat.makeCompressed();

  std::chrono::time_point<high_resolution_clock> t1, t2;
  std::chrono::duration<double, std::milli> ms_fp;

  // Suitesparse version
#ifdef GC_HAVE_SUITESPARSE

  // Convert suitesparse format
  if (internals->cMat != nullptr) {
    cholmod_l_free_sparse(&internals->cMat, internals->context);
  }
  internals->cMat = toCholmod(mat, internals->context, SType::SYMMETRIC);

  // Factor
  t1 = high_resolution_clock::now();
  internals->context.setSimplicial(); // must use simplicial for LDLt
  internals->context.setLDL();        // ensure we get an LDLt internals->factorization
  internals->factorization = cholmod_l_analyze(internals->cMat, internals->context);
  bool success = (bool)cholmod_l_factorize(internals->cMat, internals->factorization, internals->context);
  t2 = high_resolution_clock::now();
  ms_fp = t2 - t1;
  factorTime = ms_fp.count() / 1000.;

  if (!success) {
    throw std::runtime_error("failure in cholmod_l_factorize");
  }
  if (internals->context.context.status == CHOLMOD_NOT_POSDEF) {
    throw std::runtime_error("matrix is not positive definite");
  }


  // Eigen version
#else
  t1 = high_resolution_clock::now();
  internals->solver.compute(mat);
  t2 = high_resolution_clock::now();
  ms_fp = t2 - t1;
  factorTime = ms_fp.count() / 1000.;
  if (internals->solver.info() != Eigen::Success) {
    std::cerr << "Solver internals->factorization error: " << internals->solver.info() << std::endl;
    throw std::invalid_argument("Solver internals->factorization failed");
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

  size_t N = this->nRows;

  // Check some sanity
  if ((size_t)rhs.rows() != N) {
    throw std::logic_error("Vector is not the right length");
  }
#ifndef GC_NLINALG_DEBUG
  checkFinite(rhs);
#endif

  std::chrono::time_point<high_resolution_clock> t1, t2, s1, s2;
  std::chrono::duration<double, std::milli> ms_fp;

  // Suitesparse version
#ifdef GC_HAVE_SUITESPARSE

  // Convert input to suitesparse format
  s1 = high_resolution_clock::now();
  cholmod_dense* inVec = toCholmod(rhs, internals->context);
  s2 = high_resolution_clock::now();
  ms_fp = s2 - s1;
  solveConvertTime = ms_fp.count() / 1000.;

  // Solve
  t1 = high_resolution_clock::now();
  cholmod_dense* outVec = cholmod_l_solve(CHOLMOD_A, internals->factorization, inVec, internals->context);
  t2 = high_resolution_clock::now();
  ms_fp = t2 - t1;
  solveTime = ms_fp.count() / 1000.;

  // Convert back
  s1 = high_resolution_clock::now();
  toEigen(outVec, internals->context, x);

  // Free
  cholmod_l_free_dense(&outVec, internals->context);
  cholmod_l_free_dense(&inVec, internals->context);
  s2 = high_resolution_clock::now();
  ms_fp = s2 - s1;
  solveConvertTime += ms_fp.count() / 1000.;

  // Eigen version
#else
  // Solve
  t1 = high_resolution_clock::now();
  x = internals->solver.solve(rhs);
  t2 = high_resolution_clock::now();
  ms_fp = t2 - t1;
  solveTime = ms_fp.count() / 1000.;
  if (internals->solver.info() != Eigen::Success) {
    std::cerr << "Solver error: " << internals->solver.info() << std::endl;
    throw std::invalid_argument("Solve failed");
  }
#endif
}

template <typename T>
Vector<T> solvePositiveDefinite(SparseMatrix<T>& A, const Vector<T>& rhs) {
  PositiveDefiniteSolver<T> s(A);
  return s.solve(rhs);
}


// Explicit instantiations
template class PositiveDefiniteSolver<double>;
template class PositiveDefiniteSolver<float>;
template class PositiveDefiniteSolver<std::complex<double>>;

template Vector<float> solvePositiveDefinite<float>(SparseMatrix<float>& A, const Vector<float>& rhs);
template Vector<double> solvePositiveDefinite<double>(SparseMatrix<double>& A, const Vector<double>& rhs);
template Vector<std::complex<double>>
solvePositiveDefinite<std::complex<double>>(SparseMatrix<std::complex<double>>& A,
                                            const Vector<std::complex<double>>& rhs);


} // namespace geometrycentral
