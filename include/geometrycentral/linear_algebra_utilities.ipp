#include <iostream>

template <typename T>
Eigen::SparseMatrix<T> identityMatrix(size_t N) {
  Eigen::SparseMatrix<T> eye(N, N);
  eye.setIdentity();
  return eye;
}


template <typename T>
void shiftDiagonal(Eigen::SparseMatrix<T>& m, T shiftAmount) {

  // Check square
  size_t N = m.rows();
  if ((size_t)m.cols() != N) {
    throw std::logic_error("Can only shift diagonal of square matrix");
  }

  m += shiftAmount * identityMatrix<T>(N);
}


inline Eigen::SparseMatrix<double> complexToReal(const Eigen::SparseMatrix<Complex>& m) {

  size_t nRow = m.rows();
  size_t nCol = m.cols();

  Eigen::SparseMatrix<double> realM(2*nRow, 2*nCol);
  std::vector<Eigen::Triplet<double>> triplets;

  for (int k = 0; k < m.outerSize(); ++k) {
    for (typename Eigen::SparseMatrix<Complex>::InnerIterator it(m, k); it; ++it) {


      Complex val = it.value();
      size_t iRow = it.row();
      size_t iCol = it.col();

      triplets.emplace_back(2*iRow + 0, 2*iCol + 0, val.real());
      triplets.emplace_back(2*iRow + 0, 2*iCol + 1, -val.imag());
      triplets.emplace_back(2*iRow + 1, 2*iCol + 0, val.imag());
      triplets.emplace_back(2*iRow + 1, 2*iCol + 1, val.real());

    }
  }
  
  realM.setFromTriplets(triplets.begin(), triplets.end());
  realM.makeCompressed();

  return realM;
}


template <typename T>
inline void checkFinite(const Eigen::SparseMatrix<T>& m) {
  for (int k = 0; k < m.outerSize(); ++k) {
    for (typename Eigen::SparseMatrix<T>::InnerIterator it(m, k); it; ++it) {
      // std::cout << "checking "
      //<< " [" << it.row() << "," << it.col() << "] = " << it.value() << std::endl;
      if (!std::isfinite(it.value())) {
        std::cerr << std::endl
                  << "Uh oh. Non-finite matrix entry [" << it.row() << "," << it.col() << "] = " << it.value()
                  << std::endl
                  << std::endl;

        throw std::logic_error("Matrix has non-finite entries");
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
        throw std::logic_error("Matrix has non-finite entries");
      }
    }
  }
}


// Specialization of checkFinite(Matrix m) for row vectors
template <typename T, int C>
inline void checkFinite(const Eigen::Matrix<T, 1, C>& m) {
  for (unsigned int j = 0; j < m.cols(); j++) {
    if (!std::isfinite(m(1, j))) {
      std::cerr << std::endl << "Uh oh. Non-finite row vector entry [" << j << "] = " << m(j) << std::endl << std::endl;
      throw std::logic_error("Matrix has non-finite entries");
    }
  }
}


// Specialization of checkFinite(Matrix m) for column vectors
template <typename T, int R>
inline void checkFinite(const Eigen::Matrix<T, R, 1>& m) {
  for (unsigned int i = 0; i < m.rows(); i++) {
    if (!std::isfinite(m(i))) {
      std::cerr << std::endl
                << "Uh oh. Non-finite column vector entry [" << i << "] = " << m(i) << std::endl
                << std::endl;
      throw std::logic_error("Matrix has non-finite entries");
    }
  }
}

template <typename T>
inline void checkHermitian(const Eigen::SparseMatrix<T>& m) {

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

template <>
inline void checkHermitian(const Eigen::SparseMatrix<std::complex<double>>& m) {

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


template <typename T>
BlockDecompositionResult<T> blockDecomposeSquare(const Eigen::SparseMatrix<T>& m, const Vector<bool>& Aset,
                                                 bool buildBuildBside) {


  if (m.rows() != m.cols()) throw std::logic_error("blockDecomposeSquare must be called on square matrix");

  // Count sizes
  size_t initSize = m.rows();
  size_t Asize = 0;
  size_t Bsize = 0;
  for (size_t i = 0; i < initSize; i++) {
    if (Aset[i]) {
      Asize++;
    } else {
      Bsize++;
    }
  }

  // Create the result object
  BlockDecompositionResult<T> r;
  r.isA = Aset;
  r.newInds = Vector<size_t>(initSize);
  r.origIndsA = Vector<size_t>(Asize);
  r.origIndsB = Vector<size_t>(Bsize);
  r.AA = Eigen::SparseMatrix<T>(Asize, Asize);
  r.AB = Eigen::SparseMatrix<T>(Asize, Bsize);
  if (buildBuildBside) {
    r.BA = Eigen::SparseMatrix<T>(Bsize, Asize);
    r.BB = Eigen::SparseMatrix<T>(Bsize, Bsize);
  } else {
    r.BA = Eigen::SparseMatrix<T>(0, 0);
    r.BB = Eigen::SparseMatrix<T>(0, 0);
  }

  // Index
  size_t Aind = 0;
  size_t Bind = 0;
  for (size_t i = 0; i < initSize; i++) {
    if (Aset[i]) {
      r.origIndsA[Aind] = i;
      r.newInds[i] = Aind;
      Aind++;
    } else {
      r.origIndsB[Bind] = i;
      r.newInds[i] = Bind;
      Bind++;
    }
  }

  // Split
  std::vector<Eigen::Triplet<T>> AAtrip;
  std::vector<Eigen::Triplet<T>> ABtrip;
  std::vector<Eigen::Triplet<T>> BAtrip;
  std::vector<Eigen::Triplet<T>> BBtrip;
  for (size_t k = 0; k < (size_t)m.outerSize(); k++) {
    for (typename Eigen::SparseMatrix<T>::InnerIterator it(m, k); it; ++it) {

      size_t rowInd = it.row();
      size_t colInd = it.col();
      T val = it.value();

      bool rowA = Aset[rowInd];
      bool colA = Aset[colInd];

      if (rowA && colA) {
        AAtrip.emplace_back(r.newInds[rowInd], r.newInds[colInd], val);
      }
      if (rowA && !colA) {
        ABtrip.emplace_back(r.newInds[rowInd], r.newInds[colInd], val);
      }
      if (buildBuildBside) {
        if (!rowA && colA) {
          BAtrip.emplace_back(r.newInds[rowInd], r.newInds[colInd], val);
        }
        if (!rowA && !colA) {
          BBtrip.emplace_back(r.newInds[rowInd], r.newInds[colInd], val);
        }
      }
    }
  }

  // Build new matrices
  r.AA.setFromTriplets(AAtrip.begin(), AAtrip.end());
  r.AB.setFromTriplets(ABtrip.begin(), ABtrip.end());
  if (buildBuildBside) {
    r.BA.setFromTriplets(BAtrip.begin(), BAtrip.end());
    r.BB.setFromTriplets(BBtrip.begin(), BBtrip.end());
  }

  return r;
}


template <typename T>
void decomposeVector(BlockDecompositionResult<T>& decomp, const Vector<T>& vec, Vector<T>& vecAOut,
                     Vector<T>& vecBOut) {

  vecAOut = Vector<T>(decomp.origIndsA.rows());
  vecBOut = Vector<T>(decomp.origIndsB.rows());

  for (size_t i = 0; i < (size_t)vecAOut.rows(); i++) {
    vecAOut[i] = vec[decomp.origIndsA[i]];
  }
  for (size_t i = 0; i < (size_t)vecBOut.rows(); i++) {
    vecBOut[i] = vec[decomp.origIndsB[i]];
  }
}


template <typename T>
Vector<T> reassembleVector(BlockDecompositionResult<T>& decomp, const Vector<T>& vecA, const Vector<T>& vecB) {

  Vector<T> vecOut(decomp.newInds.rows());

  for (size_t i = 0; i < (size_t)vecA.rows(); i++) {
    vecOut[decomp.origIndsA[i]] = vecA[i];
  }
  for (size_t i = 0; i < (size_t)vecB.rows(); i++) {
    vecOut[decomp.origIndsB[i]] = vecB[i];
  }

  return vecOut;
}
