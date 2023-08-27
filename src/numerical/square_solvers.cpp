#include "geometrycentral/numerical/linear_solvers.h"

#include "geometrycentral/numerical/linear_algebra_utilities.h"

#ifdef GC_HAVE_SUITESPARSE
#include "geometrycentral/numerical/suitesparse_utilities.h"
#include <umfpack.h>
#endif

namespace geometrycentral {

template <typename T>
struct SquareSolverInternals {
#ifdef GC_HAVE_SUITESPARSE
  CholmodContext context;
  cholmod_sparse* cMat = nullptr;
  void* symbolicFactorization = nullptr;
  void* numericFactorization = nullptr;
#else
  Eigen::SparseLU<SparseMatrix<T>> solver;
#endif
};

template <typename T>
SquareSolver<T>::~SquareSolver() {
#ifdef GC_HAVE_SUITESPARSE
  if (internals->cMat != nullptr) {
    cholmod_l_free_sparse(&internals->cMat, internals->context);
    internals->cMat = nullptr;
  }
  if (internals->symbolicFactorization != nullptr) {
    umfpack_dl_free_symbolic(&internals->symbolicFactorization);
  }
  if (internals->numericFactorization != nullptr) {
    umfpack_dl_free_numeric(&internals->numericFactorization);
  }
#endif
}

// Helper functions to interface with umfpack without explicitly specializing all of constructor and solve(). Different
// function calls are needed for real vs. complex case.
// Note that float case is identical to double; umfpack never uses single precision
// Note that we use packed (interleaved) complex formats, rather than passing separate real and imaginary arrays
// TODO: The old SparseMatrix code in grand-central does seem to do this. Was is a bad bug, or am I missing something?
namespace {

#ifdef GC_HAVE_SUITESPARSE
// = Factorization
template <typename T>
void umfFactor(size_t N, cholmod_sparse* mat, void*& symbolicFac, void*& numericFac);

template <>
void umfFactor<double>(size_t N, cholmod_sparse* mat, void*& symbolicFac, void*& numericFac) {
  SuiteSparse_long* cMat_p = (SuiteSparse_long*)mat->p;
  SuiteSparse_long* cMat_i = (SuiteSparse_long*)mat->i;
  double* cMat_x = (double*)mat->x;
  umfpack_dl_symbolic(N, N, cMat_p, cMat_i, cMat_x, &symbolicFac, NULL, NULL);
  umfpack_dl_numeric(cMat_p, cMat_i, cMat_x, symbolicFac, &numericFac, NULL, NULL);
}
template <>
void umfFactor<float>(size_t N, cholmod_sparse* mat, void*& symbolicFac, void*& numericFac) {
  SuiteSparse_long* cMat_p = (SuiteSparse_long*)mat->p;
  SuiteSparse_long* cMat_i = (SuiteSparse_long*)mat->i;
  double* cMat_x = (double*)mat->x;
  umfpack_dl_symbolic(N, N, cMat_p, cMat_i, cMat_x, &symbolicFac, NULL, NULL);
  umfpack_dl_numeric(cMat_p, cMat_i, cMat_x, symbolicFac, &numericFac, NULL, NULL);
}
template <>
void umfFactor<std::complex<double>>(size_t N, cholmod_sparse* mat, void*& symbolicFac, void*& numericFac) {
  SuiteSparse_long* cMat_p = (SuiteSparse_long*)mat->p;
  SuiteSparse_long* cMat_i = (SuiteSparse_long*)mat->i;
  double* cMat_x = (double*)mat->x;
  umfpack_zl_symbolic(N, N, cMat_p, cMat_i, cMat_x, NULL, &symbolicFac, NULL, NULL);
  umfpack_zl_numeric(cMat_p, cMat_i, cMat_x, NULL, symbolicFac, &numericFac, NULL, NULL);
}

// = Solves
template <typename T>
void umfSolve(size_t N, cholmod_sparse* mat, void* numericFac, Vector<T>& x, const Vector<T>& rhs);

template <>
void umfSolve<float>(size_t N, cholmod_sparse* mat, void* numericFac, Vector<float>& x, const Vector<float>& rhs) {
  Vector<double> xD = Vector<double>(N);
  Vector<double> rhsD = rhs.cast<double>(); // explicitly convert to doubles
  SuiteSparse_long* cMat_p = (SuiteSparse_long*)mat->p;
  SuiteSparse_long* cMat_i = (SuiteSparse_long*)mat->i;
  double* cMat_x = (double*)mat->x;
  umfpack_dl_solve(UMFPACK_A, cMat_p, cMat_i, cMat_x, &(xD[0]), &(rhsD[0]), numericFac, NULL, NULL);
  x = xD.cast<float>();
}
template <>
void umfSolve<double>(size_t N, cholmod_sparse* mat, void* numericFac, Vector<double>& x, const Vector<double>& rhs) {
  x = Vector<double>(N);
  SuiteSparse_long* cMat_p = (SuiteSparse_long*)mat->p;
  SuiteSparse_long* cMat_i = (SuiteSparse_long*)mat->i;
  double* cMat_x = (double*)mat->x;
  umfpack_dl_solve(UMFPACK_A, cMat_p, cMat_i, cMat_x, &(x[0]), &(rhs[0]), numericFac, NULL, NULL);
}
template <>
void umfSolve<std::complex<double>>(size_t N, cholmod_sparse* mat, void* numericFac, Vector<std::complex<double>>& x,
                                    const Vector<std::complex<double>>& rhs) {
  // Note: the ordering of std::complex is specified by the standard, so this certainly works
  x = Vector<std::complex<double>>(N);
  SuiteSparse_long* cMat_p = (SuiteSparse_long*)mat->p;
  SuiteSparse_long* cMat_i = (SuiteSparse_long*)mat->i;
  double* cMat_x = (double*)mat->x;
  umfpack_zl_solve(UMFPACK_A, cMat_p, cMat_i, cMat_x, NULL, (double*)&(x[0]), NULL, (double*)&(rhs[0]), NULL,
                   numericFac, NULL, NULL);
}

#endif

} // namespace


template <typename T>
SquareSolver<T>::SquareSolver(SparseMatrix<T>& mat) : LinearSolver<T>(mat), internals(new SquareSolverInternals<T>()) {

  // Check some sanity
  if (this->nRows != this->nCols) {
    throw std::logic_error("Matrix must be square");
  }
#ifndef GC_NLINALG_DEBUG
  checkFinite(mat);
#endif

  mat.makeCompressed();

// Suitesparse variant
#ifdef GC_HAVE_SUITESPARSE
  // Convert suitesparse format
  if (internals->cMat != nullptr) {
    cholmod_l_free_sparse(&internals->cMat, internals->context);
  }
  internals->cMat = toCholmod(mat, internals->context);

  // Factor
  umfFactor<T>(this->nRows, internals->cMat, internals->symbolicFactorization, internals->numericFactorization);


// Eigen variant
#else
  internals->solver.compute(mat);
  if (internals->solver.info() != Eigen::Success) {
    std::cerr << "Solver factorization error: " << internals->solver.info() << std::endl;
    throw std::invalid_argument("Solver factorization failed");
  }
#endif
};

template <typename T>
Vector<T> SquareSolver<T>::solve(const Vector<T>& rhs) {
  Vector<T> out;
  solve(out, rhs);
  return out;
}

template <typename T>
void SquareSolver<T>::solve(Vector<T>& x, const Vector<T>& rhs) {

  size_t N = this->nRows;

  // Check some sanity
#ifndef GC_NLINALG_DEBUG
  if ((size_t)rhs.rows() != N) {
    throw std::logic_error("Vector is not the right length");
  }
  checkFinite(rhs);
#endif

  // Suitesparse version
#ifdef GC_HAVE_SUITESPARSE

  // Templated helper does all the hard work
  umfSolve<T>(N, internals->cMat, internals->numericFactorization, x, rhs);

  // Eigen version
#else
  // Solve
  x = internals->solver.solve(rhs);
  if (internals->solver.info() != Eigen::Success) {
    std::cerr << "Solver error: " << internals->solver.info() << std::endl;
    std::cerr << "Solver says: " << internals->solver.lastErrorMessage() << std::endl;
    throw std::invalid_argument("Solve failed");
  }
#endif
}

template <typename T>
Vector<T> solveSquare(SparseMatrix<T>& A, const Vector<T>& rhs) {
  SquareSolver<T> s(A);
  return s.solve(rhs);
}

// Explicit instantiations
template class SquareSolver<double>;
template class SquareSolver<float>;
template class SquareSolver<std::complex<double>>;

template Vector<float> solveSquare(SparseMatrix<float>& A, const Vector<float>& rhs);
template Vector<double> solveSquare(SparseMatrix<double>& A, const Vector<double>& rhs);
template Vector<std::complex<double>> solveSquare(SparseMatrix<std::complex<double>>& A,
                                                  const Vector<std::complex<double>>& rhs);


} // namespace geometrycentral
