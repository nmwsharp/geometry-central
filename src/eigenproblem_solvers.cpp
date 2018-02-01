#include "geometrycentral/linear_solvers.h"

#include "geometrycentral/linear_algebra_utilities.h"



using namespace Eigen;

namespace geometrycentral {

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
    double scale = std::sqrt(std::abs((x.transpose() * massMatrix * x)[0]));
    x /= scale;

    // Update
    u = x;
  }

  return x;
}

// Explicit instantiations
template Vector<double> smallestEigenvectorPositiveDefinite(SparseMatrix<double>& energyMatrix,
                                                            SparseMatrix<double>& massMatrix, size_t nIterations);
template Vector<float> smallestEigenvectorPositiveDefinite(SparseMatrix<float>& energyMatrix,
                                                            SparseMatrix<float>& massMatrix, size_t nIterations);
template Vector<Complex> smallestEigenvectorPositiveDefinite(SparseMatrix<Complex>& energyMatrix,
                                                            SparseMatrix<Complex>& massMatrix, size_t nIterations);




} // namespace geometrycentral
