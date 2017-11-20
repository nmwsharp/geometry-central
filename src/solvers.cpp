#include <solvers.h>

#include <cassert>
#include <stdexcept>

#include <linear_algebra_utilities.h>
#include <timing.h>

#include <Eigen/SparseQR>
#include <Eigen/OrderingMethods>

#ifdef HAVE_SUITESPARSE
#include <Eigen/SPQRSupport>
#endif

using std::cout;
using std::endl;
using std::cerr;

using namespace Eigen;

template <typename T>
Matrix<T, Dynamic, 1> smallestEigenvectorPositiveDefinite(
    SparseMatrix<T, ColMajor> &energyMatrix, SparseMatrix<T> &massMatrix,
    unsigned int nIterations) {
  START_TIMING(solver)
  unsigned int N = energyMatrix.rows();

  // Check some sanity
  // (none of these do anything with NDEBUG on)
  assert(energyMatrix.cols() == static_cast<int>(N) &&
         "Energy matrix is square");
  assert(massMatrix.rows() == static_cast<int>(N) &&
         "Mass matrix is same size as energy matrix");
  assert(massMatrix.cols() == static_cast<int>(N) &&
         "Mass matrix is same size as energy matrix");
  checkFinite(energyMatrix);
  checkFinite(massMatrix);
  checkSymmetric(energyMatrix);
  checkSymmetric(massMatrix);

  std::cout << std::endl
            << "Solving for smallest eigenvector with inverse power method "
               "using Eigen's sparse Cholesky LDLT"
            << std::endl;

  // Random initial vector
  // Eigen uses system rand without seeding it, so this unless you seed it the
  // behavior here will be deterministic.
  // When debugging this is a feature rather than a bug.
  Matrix<T, Dynamic, 1> u = Matrix<T, Dynamic, 1>::Random(N);

  // Compress the matrices and factor the energy matrix for faster solve
  std::cout << "  -- Factoring matrix " << std::endl;
  energyMatrix.makeCompressed();
  massMatrix.makeCompressed();
  SimplicialLDLT<SparseMatrix<T, ColMajor>> solver;
  solver.compute(energyMatrix);
  if (solver.info() != Success) {
    std::cerr << "Solver factorize error: " << solver.info() << std::endl;
    throw std::invalid_argument("Solver factorization failed");
  }

  // Perform iterative solves to converge on the solution
  std::cout << "  -- Performing " << nIterations << " inverse power iterations"
            << std::endl;
  for (unsigned int iIter = 0; iIter < nIterations; iIter++) {
    // Solve
    Matrix<T, Dynamic, 1> x = solver.solve(massMatrix * u);
    if (solver.info() != Success) {
      std::cerr << "Solver error: " << solver.info() << std::endl;
      throw std::invalid_argument("Solve failed");
    }

    // Re-normalize
    double scale = std::sqrt(std::abs((x.transpose() * massMatrix * x)[0]));
    x /= scale;

    // Update
    u = x;
  }
  std::cout << "  -- Solve complete." << endl;

  // TODO compute residuals or give some output to help us spot bad solves

  std::cout << "  -- Total time: " << pretty_time(FINISH_TIMING(solver)) << endl
            << endl;

  return u;
}
template VectorXd smallestEigenvectorPositiveDefinite(
    SparseMatrix<double, ColMajor> &energyMatrix,
    SparseMatrix<double> &massMatrix, unsigned int nIterations);
template VectorXcd smallestEigenvectorPositiveDefinite(
    SparseMatrix<std::complex<double>, ColMajor> &energyMatrix,
    SparseMatrix<std::complex<double>> &massMatrix, unsigned int nIterations);

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> solve(
    Eigen::SparseMatrix<T, Eigen::ColMajor> &A,
    Eigen::Matrix<T, Eigen::Dynamic, 1> &b) {
  START_TIMING(solver)

  A.makeCompressed();  // some solvers require compressed, and everything is
                       // faster
                       // (this is why the funciton can't be const)

  // Check some sanity
  // (none of these do anything with NDEBUG on)
  unsigned int N = A.rows();
  assert(A.cols() == static_cast<int>(N) && "Matrix is square");
  assert(b.cols() == static_cast<int>(1) && "Vector is column vector");
  assert(b.rows() == static_cast<int>(N) && "Vector is the right length");
  checkFinite(A);
  checkFinite(b);

// Call the best solver we have
#ifdef HAVE_SUITESPARSE
  Matrix<T, Dynamic, 1> x = solve_SPQR(A, b);
#else
  Matrix<T, Dynamic, 1> x = solve_EigenQR(A, b);
#endif

  std::cout << "  -- Total time: " << pretty_time(FINISH_TIMING(solver)) << endl
            << endl;

  return x;
}
template VectorXd solve(SparseMatrix<double, Eigen::ColMajor> &A, VectorXd &b);
template VectorXcd solve(SparseMatrix<std::complex<double>, Eigen::ColMajor> &A,
                         VectorXcd &b);

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> solveSquare(
    Eigen::SparseMatrix<T, Eigen::ColMajor> &A,
    Eigen::Matrix<T, Eigen::Dynamic, 1> &b) {
  START_TIMING(solver)
  unsigned int N = A.rows();

  // Check some sanity
  // (none of these do anything with NDEBUG on)
  assert(A.cols() == static_cast<int>(N) && "Matrix is square");
  assert(b.cols() == static_cast<int>(1) && "Vector is column vector");
  assert(b.rows() == static_cast<int>(N) && "Vector is the right length");
  checkFinite(A);
  checkFinite(b);

  std::cout << std::endl
            << "Solving square Ax = b with Eigen's sparse LU" << std::endl;

  // Compress the matrices and factor the energy matrix for faster solve
  std::cout << "  -- Factoring matrix " << std::endl;
  A.makeCompressed();
  Eigen::SparseLU<Eigen::SparseMatrix<T, Eigen::ColMajor>,
                  Eigen::COLAMDOrdering<int>> solver;
  solver.compute(A);
  if (solver.info() != Eigen::Success) {
    std::cerr << "Solver factorization error: " << solver.info() << std::endl;
    std::cerr << "Solver factorization says: " << solver.lastErrorMessage()
              << std::endl;
    throw std::invalid_argument("Solver factorization failed");
  }

  // Solve
  std::cout << "  -- Solving system" << std::endl;
  Matrix<T, Dynamic, 1> x = solver.solve(b);
  if (solver.info() != Eigen::Success) {
    std::cerr << "Solver error: " << solver.info() << std::endl;
    std::cerr << "Solver says: " << solver.lastErrorMessage() << std::endl;
    throw std::invalid_argument("Solve failed");
  }
  std::cout << "  -- Solve complete" << std::endl;

  // Compute residual to spot bad solves
  Matrix<T, Dynamic, 1> residual = A * x - b;
  double residualNorm = residual.norm();
  double relativeResidualNorm = residualNorm / b.norm();
  std::cout << "  -- Residual norm: " << residualNorm
            << "   relative residual norm: " << relativeResidualNorm
            << std::endl;

  std::cout << "  -- Total time: " << pretty_time(FINISH_TIMING(solver)) << endl
            << endl;

  return x;
}
template VectorXd solveSquare(SparseMatrix<double, Eigen::ColMajor> &A,
                              VectorXd &b);
template VectorXcd solveSquare(
    SparseMatrix<std::complex<double>, Eigen::ColMajor> &A, VectorXcd &b);

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> solvePositiveDefinite(
    Eigen::SparseMatrix<T, Eigen::ColMajor> &A,
    Eigen::Matrix<T, Eigen::Dynamic, 1> &b) {
  START_TIMING(solver)
  unsigned int N = A.rows();

  // Check some sanity
  // (none of these do anything with NDEBUG on)
  assert(A.cols() == static_cast<int>(N) && "Matrix is square");
  assert(b.cols() == static_cast<int>(1) && "Vector is column vector");
  assert(b.rows() == static_cast<int>(N) && "Vector is the right length");
  checkFinite(A);
  checkFinite(b);

  std::cout << std::endl
            << "Solving SPD Ax = b with Eigen's sparse Cholesky LDLT"
            << std::endl;

  // Compress the matrices and factor the energy matrix for faster solve
  std::cout << "  -- Factoring matrix " << std::endl;
  A.makeCompressed();
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<T, Eigen::ColMajor>, Eigen::Upper,
                        Eigen::COLAMDOrdering<int>> solver;
  solver.compute(A);
  if (solver.info() != Eigen::Success) {
    std::cerr << "Solver factorization error: " << solver.info() << std::endl;
    throw std::invalid_argument("Solver factorization failed");
  }

  // Solve
  std::cout << "  -- Solving system" << std::endl;
  Matrix<T, Dynamic, 1> x = solver.solve(b);
  if (solver.info() != Eigen::Success) {
    std::cerr << "Solver error: " << solver.info() << std::endl;
    // std::cerr << "Solver says: " << solver.lastErrorMessage() << std::endl;
    throw std::invalid_argument("Solve failed");
  }
  std::cout << "  -- Solve complete" << std::endl;

  // Compute residual to spot bad solves
  Matrix<T, Dynamic, 1> residual = A * x - b;
  double residualNorm = residual.norm();
  double relativeResidualNorm = residualNorm / b.norm();
  std::cout << "  -- Residual norm: " << residualNorm
            << "   relative residual norm: " << relativeResidualNorm
            << std::endl;

  std::cout << "  -- Total time: " << pretty_time(FINISH_TIMING(solver)) << endl
            << endl;

  return x;
}
template VectorXd solvePositiveDefinite(
    SparseMatrix<double, Eigen::ColMajor> &A, VectorXd &b);
template VectorXcd solvePositiveDefinite(
    SparseMatrix<std::complex<double>, Eigen::ColMajor> &A, VectorXcd &b);

// === Eigen solver variants

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> solve_EigenQR(
    Eigen::SparseMatrix<T, Eigen::ColMajor> &A,
    Eigen::Matrix<T, Eigen::Dynamic, 1> &b) {
  std::cout << std::endl
            << "Solving general Ax = b with Eigen's sparse QR" << std::endl;

  // Compress the matrices and factor the energy matrix for faster solve
  std::cout << "  -- Factoring matrix " << std::endl;
  Eigen::SparseQR<Eigen::SparseMatrix<T, Eigen::ColMajor>,
                  Eigen::COLAMDOrdering<int>> solver;
  solver.compute(A);
  if (solver.info() != Eigen::Success) {
    std::cerr << "Solver factorization error: " << solver.info() << std::endl;
    std::cerr << "Solver factorization says: " << solver.lastErrorMessage()
              << std::endl;
    throw std::invalid_argument("Solver factorization failed");
  }

  // Solve
  std::cout << "  -- Solving system" << std::endl;
  Matrix<T, Dynamic, 1> x = solver.solve(b);
  if (solver.info() != Eigen::Success) {
    std::cerr << "Solver error: " << solver.info() << std::endl;
    std::cerr << "Solver says: " << solver.lastErrorMessage() << std::endl;
    throw std::invalid_argument("Solve failed");
  }
  std::cout << "  -- Solve complete" << std::endl;

  // Compute residual to spot bad solves
  Matrix<T, Dynamic, 1> residual = A * x - b;
  double residualNorm = residual.norm();
  double relativeResidualNorm = residualNorm / b.norm();
  std::cout << "  -- Residual norm: " << residualNorm
            << "   relative residual norm: " << relativeResidualNorm
            << std::endl;

  return x;
}

// === Suitesparse solver variants

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> solve_SPQR(
    Eigen::SparseMatrix<T, Eigen::ColMajor> &A,
    Eigen::Matrix<T, Eigen::Dynamic, 1> &b) {
  std::cout << std::endl
            << "Solving general Ax = b with SuiteSparse SPQR" << std::endl;

#ifdef HAVE_SUITESPARSE

  // Compress the matrices and factor the energy matrix for faster solve
  std::cout << "  -- Factoring matrix " << std::endl;
  Eigen::SPQR<Eigen::SparseMatrix<T, Eigen::ColMajor>> solver;
  solver.compute(A);
  if (solver.info() != Eigen::Success) {
    std::cerr << "Solver factorization error: " << solver.info() << std::endl;
    throw std::invalid_argument("Solver factorization failed");
  }

  // Solve
  std::cout << "  -- Solving system" << std::endl;
  Matrix<T, Dynamic, 1> x = solver.solve(b);
  if (solver.info() != Eigen::Success) {
    std::cerr << "Solver error: " << solver.info() << std::endl;
    throw std::invalid_argument("Solve failed");
  }
  std::cout << "  -- Solve complete" << std::endl;

  // Compute residual to spot bad solves
  Matrix<T, Dynamic, 1> residual = A * x - b;
  double residualNorm = residual.norm();
  double relativeResidualNorm = residualNorm / b.norm();
  std::cout << "  -- Residual norm: " << residualNorm
            << "   relative residual norm: " << relativeResidualNorm
            << std::endl;

  return x;

#else

  throw std::invalid_argument(
      "ERROR: Attempted to call suitesparse solver in code which was compiled "
      "without suitesparse support");
  Matrix<T, Dynamic, 1> x;
  return x;

#endif
}
