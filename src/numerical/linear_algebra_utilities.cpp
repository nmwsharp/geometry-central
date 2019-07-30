#include "geometrycentral/numerical/linear_algebra_utilities.h"

namespace geometrycentral {

SparseMatrix<double> complexToReal(const SparseMatrix<std::complex<double>>& m) {

  size_t nRow = m.rows();
  size_t nCol = m.cols();

  SparseMatrix<double> realM(2 * nRow, 2 * nCol);
  std::vector<Eigen::Triplet<double>> triplets;

  for (int k = 0; k < m.outerSize(); ++k) {
    for (typename SparseMatrix<std::complex<double>>::InnerIterator it(m, k); it; ++it) {


      std::complex<double> val = it.value();
      size_t iRow = it.row();
      size_t iCol = it.col();

      triplets.emplace_back(2 * iRow + 0, 2 * iCol + 0, val.real());
      triplets.emplace_back(2 * iRow + 0, 2 * iCol + 1, -val.imag());
      triplets.emplace_back(2 * iRow + 1, 2 * iCol + 0, val.imag());
      triplets.emplace_back(2 * iRow + 1, 2 * iCol + 1, val.real());
    }
  }

  realM.setFromTriplets(triplets.begin(), triplets.end());
  realM.makeCompressed();

  return realM;
}

Vector<double> complexToReal(const Vector<std::complex<double>>& vec) {

  size_t N = vec.rows();

  Vector<double> realVec(2 * N);

  for (size_t i = 0; i < (size_t)N; i++) {
    realVec(2 * i) = vec(i).real();
    realVec(2 * i + 1) = vec(i).imag();
  }

  return realVec;
}

} // namespace geometrycentral
