#include <iostream>

template <typename T>
inline void checkFinite(const Eigen::SparseMatrix<T>& m) {
    for (int k = 0; k < m.outerSize(); ++k) {
        for (typename Eigen::SparseMatrix<T>::InnerIterator it(m, k); it; ++it) {
            if (!std::isfinite(it.value())) {
                std::cerr << std::endl
                          << "Uh oh. Non-finite matrix entry [" << it.row() << "," << it.col() << "] = " << it.value()
                          << std::endl
                          << std::endl;

                throw std::invalid_argument("Matrix has non-finite entries");
            }
        }
    }
}



// General form of checkFinite(Matrix m)
template <typename T, int R, int C>
inline void checkFinite(const Eigen::Matrix<T, R, C>& m) {
    for (unsigned int i = 0; i < m.rows(); i++) {
        for (unsigned int j = 0; j < m.cols(); j++) {
            if (!std::isfinite(m(i, j))) {
                std::cerr << std::endl
                          << "Uh oh. Non-finite vector entry [" << i << "," << j << "] = " << m(i, j) << std::endl
                          << std::endl;
                throw std::invalid_argument("Matrix has non-finite entries");
            }
        }
    }
}


// Specialization of checkFinite(Matrix m) for row vectors
template <typename T, int C>
inline void checkFinite(const Eigen::Matrix<T, 1, C>& m) {
    for (unsigned int j = 0; j < m.cols(); j++) {
        if (!std::isfinite(m(1, j))) {
            std::cerr << std::endl
                      << "Uh oh. Non-finite row vector entry [" << j << "] = " << m(j) << std::endl
                      << std::endl;
            throw std::invalid_argument("Matrix has non-finite entries");
        }
    }
}


// Specialization of checkFinite(Matrix m) for column vectors
template <typename T, int R>
inline void checkFinite(const Eigen::Matrix<T, R, 1>& m) {
    for (unsigned int i = 0; i < m.rows(); i++) {
        if (!std::isfinite(m(i))) {
            std::cerr << std::endl
                      << "Uh oh. Non-finite column vector entry [" << i << "] = " << m(i, 1) << std::endl
                      << std::endl;
            throw std::invalid_argument("Matrix has non-finite entries");
        }
    }
}

template <typename T>
void checkHermitian(const Eigen::SparseMatrix<T>& m) {

  // Compute a scale factor for the matrix to use for closeness tests
  double sum = 0;
  size_t nEntries = 0;
  for (int k = 0; k < m.outerSize(); ++k) {
    for (typename Eigen::SparseMatrix<T>::InnerIterator it(m, k); it; ++it) {
      sum += std::abs(it.value());
      nEntries++;
    }
  }
  double scale = sum / nEntries;
  double eps = scale * 1e-8;

  // Test each symmtric pair in the matrix (actually tests each twice)
  for (int k = 0; k < m.outerSize(); ++k) {
    for (typename Eigen::SparseMatrix<T>::InnerIterator it(m, k); it; ++it) {

      T thisVal = it.value();
      T otherVal = m.coeff(it.col(), it.row());

      if (!approxEqualsAbsolute(thisVal, otherVal, eps)) {
        std::cerr << std::endl
                  << "Uh oh. Non-symmtric matrix entry at [" << it.row() << "," << it.col() << "]." << std::endl
                  << "    [" << it.row() << "," << it.col() << "] = " << thisVal << std::endl
                  << "    [" << it.col() << "," << it.row() << "] = " << otherVal << std::endl
                  << std::endl
                  << std::endl;

        throw std::logic_error("Matrix has non-symmtric entries");
      }
    }
  }
}

template<>
void checkHermitian(const Eigen::SparseMatrix<std::complex<double>>& m) {

  // Compute a scale factor for the matrix to use for closeness tests
  double sum = 0;
  long long nEntries = 0;
  for (int k = 0; k < m.outerSize(); ++k) {
    for (Eigen::SparseMatrix<std::complex<double>>::InnerIterator it(m, k); it; ++it) {
      sum += std::abs(it.value());
      nEntries++;
    }
  }
  double scale = sum / nEntries;
  double eps = scale * 1e-8;

  // Test each symmtric pair in the matrix (actually tests each twice)
  for (int k = 0; k < m.outerSize(); ++k) {
    for (Eigen::SparseMatrix<std::complex<double>>::InnerIterator it(m, k); it; ++it) {

      std::complex<double> thisVal = it.value();
      std::complex<double> otherVal = m.coeff(it.col(), it.row());

      if (!approxEqualsAbsolute(thisVal, std::conj(otherVal), eps)) {
        std::cerr << std::endl
                  << "Uh oh. Non-symmtric matrix entry at [" << it.row() << "," << it.col() << "]." << std::endl
                  << "    [" << it.row() << "," << it.col() << "] = " << thisVal << std::endl
                  << "    [" << it.col() << "," << it.row() << "] = " << otherVal << std::endl
                  << std::endl
                  << std::endl;

        throw std::logic_error("Matrix has non-symmtric entries");
      }
    }
  }
}
