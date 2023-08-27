#include "geometrycentral/numerical/linear_solvers.h"

#include "geometrycentral/numerical/linear_algebra_utilities.h"


namespace geometrycentral {

// Helper
namespace {

template <typename T>
double norm(const Vector<T>& x, const SparseMatrix<T>& massMatrix) {
  return std::sqrt(std::abs(x.dot(massMatrix * x)));
}
template double norm(const Vector<double>& x, const SparseMatrix<double>& massMatrix);
template double norm(const Vector<float>& x, const SparseMatrix<float>& massMatrix);
template double norm(const Vector<std::complex<double>>& x, const SparseMatrix<std::complex<double>>& massMatrix);

template <typename T>
void normalize(Vector<T>& x, SparseMatrix<T>& massMatrix) {
  double scale = norm(x, massMatrix);
  x /= scale;
}
template void normalize(Vector<double>& x, SparseMatrix<double>& massMatrix);
template void normalize(Vector<float>& x, SparseMatrix<float>& massMatrix);
template void normalize(Vector<std::complex<double>>& x, SparseMatrix<std::complex<double>>& massMatrix);
} // namespace


template <typename T>
Vector<T> smallestEigenvectorPositiveDefinite(SparseMatrix<T>& energyMatrix, SparseMatrix<T>& massMatrix,
                                              size_t nIterations) {

  // TODO could implement a faster variant in the suitesparse case; as-is this does a copy-convert each iteration

  size_t N = energyMatrix.rows();
  PositiveDefiniteSolver<T> solver(energyMatrix);

  Vector<T> u = Vector<T>::Random(N);
  Vector<T> x = u;
  for (size_t iIter = 0; iIter < nIterations; iIter++) {

    // Solve
    solver.solve(x, massMatrix * u);

    // Re-normalize
    normalize(x, massMatrix);

    // Update
    u = x;
  }

  return x;
}

template <typename T>
std::vector<Vector<T>> smallestKEigenvectorsPositiveDefinite(SparseMatrix<T>& energyMatrix, SparseMatrix<T>& massMatrix,
                                                             size_t kEigenvalues, size_t nIterations) {

  std::vector<Vector<T>> res;

  size_t N = energyMatrix.rows();
  PositiveDefiniteSolver<T> solver(energyMatrix);

  auto projectOutPreviousVectors = [&](Vector<T>& x) {
    for (Vector<T>& v : res) {
      T proj = ((v.dot(massMatrix * x)));
      x -= proj * v;
    }
  };

  for (size_t kEig = 0; kEig < kEigenvalues; kEig++) {

    Vector<T> u = Vector<T>::Random(N);
    projectOutPreviousVectors(u);
    Vector<T> x = u;
    for (size_t iIter = 0; iIter < nIterations; iIter++) {
      // Solve
      solver.solve(x, massMatrix * u);

      projectOutPreviousVectors(x);
      normalize(x, massMatrix);

      // Update
      u = x;
    }

    res.push_back(x);
  }

  return res;
}

template <typename T>
std::vector<Vector<T>> smallestKEigenvectorsPositiveDefiniteTol(SparseMatrix<T>& energyMatrix,
                                                                SparseMatrix<T>& massMatrix, size_t kEigenvalues,
                                                                double tol) {

  std::vector<Vector<T>> res;

  size_t N = energyMatrix.rows();
  PositiveDefiniteSolver<T> solver(energyMatrix);

  auto projectOutPreviousVectors = [&](Vector<T>& x) {
    for (Vector<T>& v : res) {
      T proj = ((v.dot(massMatrix * x)));
      x -= proj * v;
    }
  };

  for (size_t kEig = 0; kEig < kEigenvalues; kEig++) {
    Vector<T> u = Vector<T>::Random(N);
    projectOutPreviousVectors(u);
    Vector<T> x = u;
    double residual = eigenvectorResidual(energyMatrix, massMatrix, x);
    while (residual > tol) {
      // Solve
      solver.solve(x, massMatrix * u);

      projectOutPreviousVectors(x);
      normalize(x, massMatrix);

      // Update
      u = x;
      residual = eigenvectorResidual(energyMatrix, massMatrix, x);
    }

    res.push_back(x);
  }

  return res;
}

template <typename T>
Vector<T> smallestEigenvectorSquare(SparseMatrix<T>& energyMatrix, SparseMatrix<T>& massMatrix, size_t nIterations) {

  // TODO could implement a faster variant in the suitesparse case; as-is this does a copy-convert each iteration

  size_t N = energyMatrix.rows();
  SquareSolver<T> solver(energyMatrix);

  Vector<T> u = Vector<T>::Random(N);
  Vector<T> x = u;
  for (size_t iIter = 0; iIter < nIterations; iIter++) {

    // Solve
    solver.solve(x, massMatrix * u);

    // Re-normalize
    normalize(x, massMatrix);

    // Update
    u = x;
  }

  return x;
}

template <typename T>
Vector<T> largestEigenvector(SparseMatrix<T>& energyMatrix, SparseMatrix<T>& massMatrix, size_t nIterations) {

  size_t N = massMatrix.rows();
  PositiveDefiniteSolver<T> solver(massMatrix);

  Vector<T> u = Vector<T>::Random(N);
  Vector<T> x = u;
  for (size_t iIter = 0; iIter < nIterations; iIter++) {
    solver.solve(x, energyMatrix * u);
    normalize(x, massMatrix);
    u = x;
  }

  return x;
}

// Measure L2 residual
template <typename T>
double eigenvectorResidual(const SparseMatrix<T>& energyMatrix, const SparseMatrix<T>& massMatrix, const Vector<T>& v) {
  T candidateEigenvalue = v.dot(energyMatrix * v);
  Vector<T> err = energyMatrix * v - candidateEigenvalue * massMatrix * v;
  return norm(err, massMatrix);
}

// Explicit instantiations
template Vector<double> smallestEigenvectorPositiveDefinite(SparseMatrix<double>& energyMatrix,
                                                            SparseMatrix<double>& massMatrix, size_t nIterations);
template Vector<float> smallestEigenvectorPositiveDefinite(SparseMatrix<float>& energyMatrix,
                                                           SparseMatrix<float>& massMatrix, size_t nIterations);
template Vector<std::complex<double>>
smallestEigenvectorPositiveDefinite(SparseMatrix<std::complex<double>>& energyMatrix,
                                    SparseMatrix<std::complex<double>>& massMatrix, size_t nIterations);

template std::vector<Vector<float>> smallestKEigenvectorsPositiveDefinite(SparseMatrix<float>& energyMatrix,
                                                                          SparseMatrix<float>& massMatrix,
                                                                          size_t kEigenvalues, size_t nIterations);
template std::vector<Vector<double>> smallestKEigenvectorsPositiveDefinite(SparseMatrix<double>& energyMatrix,
                                                                           SparseMatrix<double>& massMatrix,
                                                                           size_t kEigenvalues, size_t nIterations);
template std::vector<Vector<std::complex<double>>>
smallestKEigenvectorsPositiveDefinite(SparseMatrix<std::complex<double>>& energyMatrix,
                                      SparseMatrix<std::complex<double>>& massMatrix, size_t kEigenvalues,
                                      size_t nIterations);


template std::vector<Vector<float>> smallestKEigenvectorsPositiveDefiniteTol(SparseMatrix<float>& energyMatrix,
                                                                             SparseMatrix<float>& massMatrix,
                                                                             size_t kEigenvalues, double tol);
template std::vector<Vector<double>> smallestKEigenvectorsPositiveDefiniteTol(SparseMatrix<double>& energyMatrix,
                                                                              SparseMatrix<double>& massMatrix,
                                                                              size_t kEigenvalues, double tol);
template std::vector<Vector<std::complex<double>>>
smallestKEigenvectorsPositiveDefiniteTol(SparseMatrix<std::complex<double>>& energyMatrix,
                                         SparseMatrix<std::complex<double>>& massMatrix, size_t kEigenvalues,
                                         double tol);

template Vector<double> smallestEigenvectorSquare(SparseMatrix<double>& energyMatrix, SparseMatrix<double>& massMatrix,
                                                  size_t nIterations);
template Vector<float> smallestEigenvectorSquare(SparseMatrix<float>& energyMatrix, SparseMatrix<float>& massMatrix,
                                                 size_t nIterations);
template Vector<std::complex<double>> smallestEigenvectorSquare(SparseMatrix<std::complex<double>>& energyMatrix,
                                                                SparseMatrix<std::complex<double>>& massMatrix,
                                                                size_t nIterations);

template Vector<double> largestEigenvector(SparseMatrix<double>& energyMatrix, SparseMatrix<double>& massMatrix,
                                           size_t nIterations);
template Vector<float> largestEigenvector(SparseMatrix<float>& energyMatrix, SparseMatrix<float>& massMatrix,
                                          size_t nIterations);
template Vector<std::complex<double>> largestEigenvector(SparseMatrix<std::complex<double>>& energyMatrix,
                                                         SparseMatrix<std::complex<double>>& massMatrix,
                                                         size_t nIterations);


} // namespace geometrycentral
