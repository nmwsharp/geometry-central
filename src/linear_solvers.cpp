#include "geometrycentral/linear_solvers.h"

#include "geometrycentral/linear_algebra_utilities.h"


// This control various safety chceks in linear algebra code and solvers
// #define NLINALG_DEBUG;

using namespace Eigen;

namespace geometrycentral {

template <typename T>
Vector<T> smallestEigenvectorPositiveDefinite(SparseMatrix<T>& energyMatrix, SparseMatrix<T>& massMatrix,
                                              size_t nIterations) {

  size_t N = energyMatrix.rows();
  PositiveDefiniteSolver<T> solver(energyMatrix);

  Vector<T> u = Vector<T>::Random(N);
  Vector<T> x = u;
  for (size_t iIter = 0; iIter < nIterations; iIter++) {

    // Solve
    solver(x, massMatrix * u);

    // Re-normalize
    double scale = std::sqrt(std::abs((x.transpose() * massMatrix * x)[0]));
    x /= scale;

    // Update
    u = x;
  }

  return x;
}
template Vector<double> smallestEigenvectorPositiveDefinite(SparseMatrix<double>& energyMatrix,
                                                            SparseMatrix<double>& massMatrix, size_t nIterations);
template Vector<float> smallestEigenvectorPositiveDefinite(SparseMatrix<float>& energyMatrix,
                                                            SparseMatrix<float>& massMatrix, size_t nIterations);
template Vector<Complex> smallestEigenvectorPositiveDefinite(SparseMatrix<Complex>& energyMatrix,
                                                            SparseMatrix<Complex>& massMatrix, size_t nIterations);


template <typename T>
Vector<T> LinearSolver<T>::operator()(const Vector<T>& rhs) {
  Vector<T> lhs;
  (*this)(lhs, rhs);
  return lhs;
}
template class LinearSolver<double>;


// ========================================================
// ==========          Positive Definite         ==========
// ========================================================

template <typename T>
void PositiveDefiniteSolver<T>::prepare() {

  size_t N = this->mat.rows();

// Check some sanity
#ifndef NLINALG_DEBUG
  assert(this->mat.cols() == static_cast<int>(N) && "Matrix is square");
  checkFinite(this->mat);
  checkHermitian(this->mat);
#endif

  this->mat.makeCompressed();
  solver.compute(this->mat);
  if (solver.info() != Eigen::Success) {
    std::cerr << "Solver factorization error: " << solver.info() << std::endl;
    throw std::invalid_argument("Solver factorization failed");
  }
};

template <typename T>
void PositiveDefiniteSolver<T>::operator()(Vector<T>& x, const Vector<T>& rhs) {

  size_t N = this->mat.rows();

  // Check some sanity
#ifndef NLINALG_DEBUG
  assert(rhs.rows() == static_cast<int>(N) && "Vector is the right length");
  checkFinite(rhs);
#endif

  // Solve
  x = solver.solve(rhs);
  if (solver.info() != Eigen::Success) {
    std::cerr << "Solver error: " << solver.info() << std::endl;
    // std::cerr << "Solver says: " << solver.lastErrorMessage() << std::endl;
    throw std::invalid_argument("Solve failed");
  }

    // Compute residual to spot bad solves
#ifndef NLINALG_DEBUG
  Matrix<T, Dynamic, 1> residual = this->mat * x - rhs;
  double residualNorm = residual.norm();
  double relativeResidualNorm = residualNorm / rhs.norm();
  std::cout << "  -- Residual norm: " << residualNorm << "   relative residual norm: " << relativeResidualNorm
            << std::endl;
#endif
}

template <typename T>
Vector<T> PositiveDefiniteSolver<T>::solve(const Eigen::SparseMatrix<T>& A, const Vector<T>& rhs) {
  PositiveDefiniteSolver<T> s(A);
  return static_cast<LinearSolver<T>*>(&s)->operator()(rhs); // lol?
}

template <typename T>
void PositiveDefiniteSolver<T>::solve(const Eigen::SparseMatrix<T>& A, Vector<T>& x, const Vector<T>& rhs) {
  PositiveDefiniteSolver<T> s(A);
  s(x, rhs);
}

// Explicit instantiations
template class PositiveDefiniteSolver<double>;
template class PositiveDefiniteSolver<float>;
template class PositiveDefiniteSolver<Complex>;

// ========================================================
// ==========               Square               ==========
// ========================================================

template <typename T>
void SquareSolver<T>::prepare() {

  size_t N = this->mat.rows();

// Check some sanity
#ifndef NLINALG_DEBUG
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
#ifndef NLINALG_DEBUG
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
#ifndef NLINALG_DEBUG
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