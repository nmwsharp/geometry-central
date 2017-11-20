#include <iostream>

template <typename T>
inline void checkFinite(const Eigen::SparseMatrix<T>& m) {
#ifndef NDEBUG
  for (int k = 0; k < m.outerSize(); ++k) {
    for (typename Eigen::SparseMatrix<T>::InnerIterator it(m, k); it; ++it) {
      if (!std::isfinite(it.value())) {
        std::cerr << std::endl
                  << "Uh oh. Non-finite matrix entry [" << it.row() << ","
                  << it.col() << "] = " << it.value() << std::endl
                  << std::endl;

        throw std::invalid_argument("Matrix has non-finite entries");
      }
    }
  }
#endif
}

// General form of checkFinite(Matrix m)
template <typename T, int R, int C>
inline void checkFinite(const Eigen::Matrix<T, R, C>& m) {
#ifndef NDEBUG
  for (unsigned int i = 0; i < m.rows(); i++) {
    for (unsigned int j = 0; j < m.cols(); j++) {
      if (!std::isfinite(m(i, j))) {
        std::cerr << std::endl
                  << "Uh oh. Non-finite vector entry [" << i << "," << j
                  << "] = " << m(i, j) << std::endl
                  << std::endl;
        throw std::invalid_argument("Matrix has non-finite entries");
      }
    }
  }
#endif
}

// Specialization of checkFinite(Matrix m) for row vectors
template <typename T, int C>
inline void checkFinite(const Eigen::Matrix<T, 1, C>& m) {
#ifndef NDEBUG
  for (unsigned int j = 0; j < m.cols(); j++) {
    if (!std::isfinite(m(1, j))) {
      std::cerr << std::endl
                << "Uh oh. Non-finite row vector entry [" << j << "] = " << m(j)
                << std::endl
                << std::endl;
      throw std::invalid_argument("Matrix has non-finite entries");
    }
  }
#endif
}

// Specialization of checkFinite(Matrix m) for column vectors
template <typename T, int R>
inline void checkFinite(const Eigen::Matrix<T, R, 1>& m) {
#ifndef NDEBUG
  for (unsigned int i = 0; i < m.rows(); i++) {
    if (!std::isfinite(m(i))) {
      std::cerr << std::endl
                << "Uh oh. Non-finite column vector entry [" << i
                << "] = " << m(i, 1) << std::endl
                << std::endl;
      throw std::invalid_argument("Matrix has non-finite entries");
    }
  }
#endif
}
