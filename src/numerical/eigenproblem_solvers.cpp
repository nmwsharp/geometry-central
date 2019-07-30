#include "geometrycentral/numerical/linear_solvers.h"

#include "geometrycentral/numerical/linear_algebra_utilities.h"


using namespace Eigen;


namespace geometrycentral {

// Helper
namespace {

template <typename T>
double norm(Vector<T>& x, SparseMatrix<T>& massMatrix) {
  return std::sqrt(std::abs((x.transpose() * massMatrix * x)[0]));
}
template double norm(Vector<double>& x, SparseMatrix<double>& massMatrix);
template double norm(Vector<float>& x, SparseMatrix<float>& massMatrix);
template double norm(Vector<std::complex<double>>& x, SparseMatrix<std::complex<double>>& massMatrix);

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
      T proj = (((x.transpose() * massMatrix * v)[0]));
      x -= v * proj;
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
