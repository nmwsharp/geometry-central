#include <linear_algebra_utilities.h>

#include <utilities.h>

void checkSymmetric(const Eigen::SparseMatrix<double>& m) {
#ifndef NDEBUG

  // Compute a scale factor for the matrix to use for closeness tests
  double sum = 0;
  long long nEntries = 0;
  for (int k = 0; k < m.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(m, k); it; ++it) {
      sum += std::abs(it.value());
      nEntries++;
    }
  }
  double scale = sum / nEntries;
  double eps = scale * 1e-8;

  // Test each symmtric pair in the matrix (actually tests each twice)
  for (int k = 0; k < m.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(m, k); it; ++it) {
      double thisVal = it.value();
      double otherVal = m.coeff(it.col(), it.row());

      if (!approxEqualsAbsolute(thisVal, otherVal, eps)) {
        std::cerr << std::endl
                  << "Uh oh. Non-symmtric matrix entry at [" << it.row() << ","
                  << it.col() << "]." << std::endl
                  << "    [" << it.row() << "," << it.col() << "] = " << thisVal
                  << std::endl
                  << "    [" << it.col() << "," << it.row()
                  << "] = " << otherVal << std::endl
                  << std::endl
                  << std::endl;

        throw std::invalid_argument("Matrix has non-symmtric entries");
      }
    }
  }
#endif
}

void checkSymmetric(const Eigen::SparseMatrix<std::complex<double>>& m) {
#ifndef NDEBUG

  // Compute a scale factor for the matrix to use for closeness tests
  double sum = 0;
  long long nEntries = 0;
  for (int k = 0; k < m.outerSize(); ++k) {
    for (Eigen::SparseMatrix<std::complex<double>>::InnerIterator it(m, k); it;
         ++it) {
      sum += std::abs(it.value());
      nEntries++;
    }
  }
  double scale = sum / nEntries;
  double eps = scale * 1e-8;

  // Test each symmtric pair in the matrix (actually tests each twice)
  for (int k = 0; k < m.outerSize(); ++k) {
    for (Eigen::SparseMatrix<std::complex<double>>::InnerIterator it(m, k); it;
         ++it) {
      std::complex<double> thisVal = it.value();
      std::complex<double> otherVal = m.coeff(it.col(), it.row());

      if (!approxEqualsAbsolute(thisVal, std::conj(otherVal), eps)) {
        std::cerr << std::endl
                  << "Uh oh. Non-symmtric matrix entry at [" << it.row() << ","
                  << it.col() << "]." << std::endl
                  << "    [" << it.row() << "," << it.col() << "] = " << thisVal
                  << std::endl
                  << "    [" << it.col() << "," << it.row()
                  << "] = " << otherVal << std::endl
                  << std::endl
                  << std::endl;

        throw std::invalid_argument("Matrix has non-symmtric entries");
      }
    }
  }
#endif
}
