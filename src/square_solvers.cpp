#include "geometrycentral/linear_solvers.h"

#include "geometrycentral/linear_algebra_utilities.h"


using namespace Eigen;

namespace geometrycentral {


template <typename T>
void SquareSolver<T>::prepare() {

  size_t N = this->mat.rows();

// Check some sanity
#ifndef GC_NLINALG_DEBUG
  assert(this->mat.cols() == static_cast<int>(N) && "Matrix is square");
  checkFinite(this->mat);
#endif

  this->mat.makeCompressed();
  solver.compute(this->mat);
  if (solver.info() != Eigen::Success) {
    std::cerr << "Solver factorization error: " << solver.info() << std::endl;
    throw std::invalid_argument("Solver factorization failed");
  }
};

template <typename T>
void SquareSolver<T>::operator()(Vector<T>& x, const Vector<T>& rhs) {

  size_t N = this->mat.rows();

  // Check some sanity
#ifndef GC_NLINALG_DEBUG
  assert(rhs.rows() == static_cast<int>(N) && "Vector is the right length");
  checkFinite(rhs);
#endif

  // Solve
  x = solver.solve(rhs);
  if (solver.info() != Eigen::Success) {
    std::cerr << "Solver error: " << solver.info() << std::endl;
    std::cerr << "Solver says: " << solver.lastErrorMessage() << std::endl;
    throw std::invalid_argument("Solve failed");
  }

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
Vector<T> SquareSolver<T>::solve(const Eigen::SparseMatrix<T>& A, const Vector<T>& rhs) {
  SquareSolver<T> s(A);
  return static_cast<LinearSolver<T>*>(&s)->operator()(rhs); // lol?
}

template <typename T>
void SquareSolver<T>::solve(const Eigen::SparseMatrix<T>& A, Vector<T>& x, const Vector<T>& rhs) {
  SquareSolver<T> s(A);
  s(x, rhs);
}

// Explicit instantiations
template class SquareSolver<double>;
template class SquareSolver<float>;
template class SquareSolver<Complex>;


} // namespace geometrycentral