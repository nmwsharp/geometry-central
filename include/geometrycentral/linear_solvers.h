#pragma once

#include "Eigen/Sparse"

// Note: actual solvers implemented with explicit template instantiation in solvers.cpp

// Utilitiy typedef
namespace {
template <typename T>
using Vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;
}

namespace geometrycentral {

// Utility solvers, which use the classes below
template <typename T>
Vector<T> smallestEigenvectorPositiveDefinite(Eigen::SparseMatrix<T>& energyMatrix, Eigen::SparseMatrix<T>& massMatrix,
                                              size_t nIterations = 50);

// Base class for all linear solvers
template <typename T>
class LinearSolver {

public:
  LinearSolver(const Eigen::SparseMatrix<T>& mat_) {
    mat = Eigen::SparseMatrix<T, Eigen::ColMajor>(mat_);
  }

  // Solve for a particular right hand side
  Vector<T> operator()(const Vector<T>& rhs);

  // Solve for a particular right hand side, and return in an existing vector objects
  virtual void operator()(Vector<T>& x, const Vector<T>& rhs) = 0;

  // Compute the residual of a solve
  double residual(const Vector<T>& lhs, const Vector<T>& rhs);

protected:
  Eigen::SparseMatrix<T, Eigen::ColMajor> mat;

};


template <typename T>
class PositiveDefiniteSolver : public LinearSolver<T> {

public:
  PositiveDefiniteSolver(const Eigen::SparseMatrix<T>& mat_) : LinearSolver<T>(mat_) {
    prepare(); 
  }

  // Solve!
  virtual void operator()(Vector<T>& x, const Vector<T>& rhs) override;

  // Static solve without retaining the factorization
  static Vector<T> solve(const Eigen::SparseMatrix<T>& A, const Vector<T>& rhs);
  static void solve(const Eigen::SparseMatrix<T>& A, Vector<T>& x, const Vector<T>& rhs);

protected:
  void prepare();
  
  // Implementation-specific quantities
  // #ifdef HAVE_SUITESPARSE

  // #else
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<T>> solver;
  // #endif
};

template <typename T>
class SquareSolver : public LinearSolver<T> {

public:
  SquareSolver(const Eigen::SparseMatrix<T>& mat_) : LinearSolver<T>(mat_) {
    prepare(); 
  }

  // Solve!
  virtual void operator()(Vector<T>& x, const Vector<T>& rhs) override;

  // Static solve without retaining the factorization
  static Vector<T> solve(const Eigen::SparseMatrix<T>& A, const Vector<T>& rhs);
  static void solve(const Eigen::SparseMatrix<T>& A, Vector<T>& x, const Vector<T>& rhs);

protected:
  void prepare();
  
  // Implementation-specific quantities
  // #ifdef HAVE_SUITESPARSE

  // #else
  Eigen::SparseLU<Eigen::SparseMatrix<T>> solver;
  // #endif
};

} // namespace geometrycentral
