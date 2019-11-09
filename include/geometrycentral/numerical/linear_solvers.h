#pragma once

#include "geometrycentral/numerical/linear_algebra_utilities.h"

#include "Eigen/Sparse"

#include <iostream>
#include <memory>

// This disables various safety checks in linear algebra code and solvers
// #define GC_NLINALG_DEBUG

// Note: actual solvers implemented with explicit template instantiation in solvers.cpp

namespace geometrycentral {

// === Utility solvers, which use the classes below

// Returns smallest nontrivial eigenvector
template <typename T>
Vector<T> smallestEigenvectorPositiveDefinite(SparseMatrix<T>& energyMatrix, SparseMatrix<T>& massMatrix,
                                              size_t nIterations = 50);

template <typename T>
std::vector<Vector<T>> smallestKEigenvectorsPositiveDefinite(SparseMatrix<T>& energyMatrix, SparseMatrix<T>& massMatrix,
                                                             size_t kEigenvalues, size_t nIterations = 50);
template <typename T>
std::vector<Vector<T>> smallestKEigenvectorsPositiveDefiniteTol(SparseMatrix<T>& energyMatrix,
                                                                SparseMatrix<T>& massMatrix, size_t kEigenvalues,
                                                                double tol = 1e-8);

// Returns smallest (positive-eigenvalued) nontirivial eigenvector
template <typename T>
Vector<T> smallestEigenvectorSquare(SparseMatrix<T>& energyMatrix, SparseMatrix<T>& massMatrix,
                                    size_t nIterations = 50);

// Mass matrix must be positive definite
template <typename T>
Vector<T> largestEigenvector(SparseMatrix<T>& energyMatrix, SparseMatrix<T>& massMatrix, size_t nIterations = 50);

// Measure L2 residual
template <typename T>
double eigenvectorResidual(const SparseMatrix<T>& energyMatrix, const SparseMatrix<T>& massMatrix, const Vector<T>& v);

// Quick and easy solvers which do not retain factorization
template <typename T>
Vector<T> solve(SparseMatrix<T>& matrix, const Vector<T>& rhs);
template <typename T>
Vector<T> solveSquare(SparseMatrix<T>& matrix, const Vector<T>& rhs);
template <typename T>
Vector<T> solvePositiveDefinite(SparseMatrix<T>& matrix, const Vector<T>& rhs);

// Measure L2 residual
template <typename T>
double residual(const SparseMatrix<T>& matrix, const Vector<T>& lhs, const Vector<T>& rhs);

// Base class for all linear solvers
template <typename T>
class LinearSolver {

public:
  LinearSolver(const SparseMatrix<T>& mat) : nRows(mat.rows()), nCols(mat.cols()) {}
  virtual ~LinearSolver() {}

  // Solve for a particular right hand side
  virtual Vector<T> solve(const Vector<T>& rhs) = 0;

  // Solve for a particular right hand side, and return in an existing vector objects
  virtual void solve(Vector<T>& x, const Vector<T>& rhs) = 0;

protected:
  size_t nRows, nCols;
};

// General solver (uses QR)
// Computes least-squares solution for overdetermined systems, minimum norm solution for underdetermined systems
// TODO name is dumb
template <typename T>
struct QRSolverInternals; // hide implementation details
template <typename T>
class Solver final : public LinearSolver<T> {

public:
  Solver(SparseMatrix<T>& mat);
  ~Solver();

  // Solve!
  void solve(Vector<T>& x, const Vector<T>& rhs) override;
  Vector<T> solve(const Vector<T>& rhs) override;

  // Gets the rank of the system
  size_t rank();

protected:
  bool underdetermined;
  std::unique_ptr<QRSolverInternals<T>> internals;
};

template <typename T>
struct PSDSolverInternals; // hide implementation details
template <typename T>
class PositiveDefiniteSolver final : public LinearSolver<T> {

public:
  PositiveDefiniteSolver(SparseMatrix<T>& mat);
  ~PositiveDefiniteSolver();

  // Solve!
  void solve(Vector<T>& x, const Vector<T>& rhs) override;
  Vector<T> solve(const Vector<T>& rhs) override;

protected:
  std::unique_ptr<PSDSolverInternals<T>> internals;
};

template <typename T>
struct SquareSolverInternals; // hide implementation details
template <typename T>
class SquareSolver final : public LinearSolver<T> {

public:
  SquareSolver(SparseMatrix<T>& mat);
  ~SquareSolver();

  // Solve!
  void solve(Vector<T>& x, const Vector<T>& rhs) override;
  Vector<T> solve(const Vector<T>& rhs) override;

protected:
  // Implementation-specific quantities
  std::unique_ptr<SquareSolverInternals<T>> internals;
};

} // namespace geometrycentral
