template <typename T>
SparseMatrix<T> identityMatrix(size_t N) {
  SparseMatrix<T> eye(N, N);
  eye.setIdentity();
  return eye;
}


template <typename T>
void shiftDiagonal(SparseMatrix<T>& m, T shiftAmount) {

  // Check square
  size_t N = m.rows();
  if ((size_t)m.cols() != N) {
    throw std::logic_error("Can only shift diagonal of square matrix");
  }

  m += shiftAmount * identityMatrix<T>(N);
}


template <typename T>
SparseMatrix<T> verticalStack(const std::vector<SparseMatrix<T>, Eigen::aligned_allocator<SparseMatrix<T>>>& mats) {
  if (mats.size() == 0) throw std::logic_error("must have at least one matrix to stack");

  std::vector<Eigen::Triplet<T>> triplets;
  long nCols = mats[0].cols();

  size_t nRowsTot = 0;

  for (const SparseMatrix<T>& mat : mats) {

    if (mat.cols() != nCols) throw std::logic_error("all matrices must have same column size");

    // Copy entries
    for (int k = 0; k < mat.outerSize(); ++k) {
      for (typename SparseMatrix<T>::InnerIterator it(mat, k); it; ++it) {

        T thisVal = it.value();
        int i = it.row();
        int j = it.col();

        triplets.emplace_back(i + nRowsTot, j, thisVal);
      }
    }

    nRowsTot += mat.rows();
  }

  // Build the matrix
  SparseMatrix<T> result(nRowsTot, nCols);
  result.setFromTriplets(triplets.begin(), triplets.end());

  return result;
}

template <typename T>
SparseMatrix<T> horizontalStack(const std::vector<SparseMatrix<T>, Eigen::aligned_allocator<SparseMatrix<T>>>& mats) {
  if (mats.size() == 0) throw std::logic_error("must have at least one matrix to stack");

  std::vector<Eigen::Triplet<T>> triplets;
  long nRows = mats[0].rows();

  size_t nColsTot = 0;

  for (const SparseMatrix<T>& mat : mats) {

    if (mat.rows() != nRows) throw std::logic_error("all matrices must have same row size");

    // Copy entries
    for (int k = 0; k < mat.outerSize(); ++k) {
      for (typename SparseMatrix<T>::InnerIterator it(mat, k); it; ++it) {

        T thisVal = it.value();
        int i = it.row();
        int j = it.col();

        triplets.emplace_back(i, j + nColsTot, thisVal);
      }
    }

    nColsTot += mat.cols();
  }

  // Build the matrix
  SparseMatrix<T> result(nRows, nColsTot);
  result.setFromTriplets(triplets.begin(), triplets.end());

  return result;
}


template <typename T>
std::vector<std::vector<T>> unpackMatrixToStdVector(const DenseMatrix<T>& mat) {
  size_t N = static_cast<size_t>(mat.rows());
  size_t M = static_cast<size_t>(mat.cols());

  // Copy to vector representation
  std::vector<std::vector<size_t>> vectors(N);
  for (size_t i = 0; i < N; i++) {
    vectors[i].resize(M);
    for (size_t j = 0; j < M; j++) {
      vectors[i][j] = mat(i, j);
    }
  }

  return vectors;
}

template <typename T>
inline void checkFinite(const SparseMatrix<T>& m) {
  for (int k = 0; k < m.outerSize(); ++k) {
    for (typename SparseMatrix<T>::InnerIterator it(m, k); it; ++it) {
      if (!isfinite(it.value())) {
        std::ostringstream msg;
        msg << "checkFinite() failure: Non-finite matrix entry [" << it.row() << "," << it.col()
            << "] = " << it.value();
        // std::cerr << msg.str() << std::endl;
        throw std::logic_error(msg.str());
      }
    }
  }
}


// General form of checkFinite(Matrix m)
template <typename T, int R, int C>
inline void checkFinite(const Eigen::Matrix<T, R, C>& m) {
  for (unsigned int i = 0; i < m.rows(); i++) {
    for (unsigned int j = 0; j < m.cols(); j++) {
      if (!isfinite(m(i, j))) {
        std::ostringstream msg;
        msg << "checkFinite() failure. Non-finite vector entry [" << i << "," << j << "] = " << m(i, j);
        // std::cerr << msg.str() << std::endl;
        throw std::logic_error(msg.str());
      }
    }
  }
}


// Specialization of checkFinite(Matrix m) for row vectors
template <typename T, int C>
inline void checkFinite(const Eigen::Matrix<T, 1, C>& m) {
  for (unsigned int j = 0; j < m.cols(); j++) {
    if (!std::isfinite(m(1, j))) {
      std::ostringstream msg;
      msg << "checkFinite() failure. Non-finite row vector entry [" << j << "] = " << m(j);
      // std::cerr << msg.str() << std::endl;
      throw std::logic_error(msg.str());
    }
  }
}


// Specialization of checkFinite(Matrix m) for column vectors
template <typename T, int R>
inline void checkFinite(const Eigen::Matrix<T, R, 1>& m) {

  for (unsigned int i = 0; i < m.rows(); i++) {
    if (!isfinite(m(i))) {
      std::ostringstream msg;
      msg << "checkFinite() failure. Non-finite row vector entry [" << i << "] = " << m(i);
      // std::cerr << msg.str() << std::endl;
      throw std::logic_error(msg.str());
    }
  }
}

template <typename T>
inline void checkSymmetric(const SparseMatrix<T>& m, double absoluteEPS) {

  double eps = absoluteEPS;
  // Compute a scale factor for the matrix to use for closeness tests
  if (eps == -1.) {
    double sum = 0;
    size_t nEntries = 0;
    for (int k = 0; k < m.outerSize(); ++k) {
      for (typename SparseMatrix<T>::InnerIterator it(m, k); it; ++it) {
        sum += std::abs(it.value());
        nEntries++;
      }
    }
    double scale = sum / nEntries;
    eps = scale * 1e-8;
  }

  // Test each symmtric pair in the matrix (actually tests each twice)
  for (int k = 0; k < m.outerSize(); ++k) {
    for (typename SparseMatrix<T>::InnerIterator it(m, k); it; ++it) {

      T thisVal = it.value();
      T otherVal = m.coeff(it.col(), it.row());

      if (std::abs(thisVal - otherVal) > eps) {
        std::ostringstream msg;
        msg << "checkSymmetric() error. Non-symmtric matrix entry at [" << it.row() << "," << it.col() << "]."
            << "    [" << it.row() << "," << it.col() << "] = " << thisVal << "    [" << it.col() << "," << it.row()
            << "] = " << otherVal;
        // std::cerr << msg.str() << std::endl;
        throw std::logic_error(msg.str());
      }
    }
  }
}


template <typename T>
inline void checkHermitian(const SparseMatrix<T>& m, double absoluteEPS) {

  double eps = absoluteEPS;
  // Compute a scale factor for the matrix to use for closeness tests
  if (eps == -1.) {
    double sum = 0;
    size_t nEntries = 0;
    for (int k = 0; k < m.outerSize(); ++k) {
      for (typename SparseMatrix<T>::InnerIterator it(m, k); it; ++it) {
        sum += std::abs(it.value());
        nEntries++;
      }
    }
    double scale = sum / nEntries;
    eps = scale * 1e-8;
  }

  // Test each symmtric pair in the matrix (actually tests each twice)
  for (int k = 0; k < m.outerSize(); ++k) {
    for (typename SparseMatrix<T>::InnerIterator it(m, k); it; ++it) {

      T thisVal = it.value();
      T otherVal = m.coeff(it.col(), it.row());

      if (std::abs(thisVal - conj(otherVal)) > eps) {
        std::ostringstream msg;
        msg << "checkHermitian() error. Non-Hermitian matrix entry at [" << it.row() << "," << it.col() << "]."
            << "    [" << it.row() << "," << it.col() << "] = " << thisVal << "    [" << it.col() << "," << it.row()
            << "] = " << otherVal;
        // std::cerr << msg.str() << std::endl;
        throw std::logic_error(msg.str());
      }
    }
  }
}


template <typename T>
BlockDecompositionResult<T> blockDecomposeSquare(const SparseMatrix<T>& m, const Vector<bool>& Aset,
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
  r.AA = SparseMatrix<T>(Asize, Asize);
  r.AB = SparseMatrix<T>(Asize, Bsize);
  if (buildBuildBside) {
    r.BA = SparseMatrix<T>(Bsize, Asize);
    r.BB = SparseMatrix<T>(Bsize, Bsize);
  } else {
    r.BA = SparseMatrix<T>(0, 0);
    r.BB = SparseMatrix<T>(0, 0);
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
    for (typename SparseMatrix<T>::InnerIterator it(m, k); it; ++it) {

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

template <typename T>
void saveSparseMatrix(std::string filename, SparseMatrix<T>& matrix) {
  // WARNING: this follows matlab convention and thus is 1-indexed

  std::ofstream outFile(filename);
  if (!outFile) {
    throw std::runtime_error("failed to open output file " + filename);
  }
  saveSparseMatrix(outFile, matrix);
  outFile.close();
}

template <typename T>
void saveSparseMatrix(std::ostream& out, SparseMatrix<T>& matrix) {
  // Write a comment on the first line giving the dimensions
  out << "# sparse " << matrix.rows() << " " << matrix.cols() << std::endl;

  out << std::setprecision(16);

  for (int k = 0; k < matrix.outerSize(); ++k) {
    for (typename SparseMatrix<T>::InnerIterator it(matrix, k); it; ++it) {
      T val = it.value();
      size_t iRow = it.row();
      size_t iCol = it.col();

      out << (iRow + 1) << " " << (iCol + 1) << " " << val << std::endl;
    }
  }
}

template <typename T>
void saveDenseMatrix(std::string filename, DenseMatrix<T>& matrix) {
  std::ofstream outFile(filename);
  if (!outFile) {
    throw std::runtime_error("failed to open output file " + filename);
  }

  saveDenseMatrix(outFile, matrix);
  outFile.close();
}

template <typename T>
void saveDenseMatrix(std::ostream& out, DenseMatrix<T>& matrix) {
  // Write a comment on the first line giving the dimensions
  out << "# dense " << matrix.rows() << " " << matrix.cols() << std::endl;

  out << std::setprecision(16);

  for (size_t iRow = 0; iRow < (size_t)matrix.rows(); iRow++) {
    for (size_t iCol = 0; iCol < (size_t)matrix.cols(); iCol++) {
      T val = matrix(iRow, iCol);
      out << val;
      if (iCol + 1 != (size_t)matrix.cols()) {
        out << " ";
      }
    }
    out << std::endl;
  }
}


template <typename T>
SparseMatrix<T> loadSparseMatrix(std::string filename) {
  std::ifstream inFile(filename);
  if (!inFile) {
    throw std::runtime_error("failed to open input file " + filename);
  }

  SparseMatrix<T> M;
  try {
    M = loadSparseMatrix<T>(inFile);
  } catch (const std::runtime_error& e) {
    throw std::runtime_error("failed to parse sparse matrix " + filename);
  }

  inFile.close();
  return M;
}

template <typename T>
SparseMatrix<T> loadSparseMatrix(std::istream& in) {
  std::string line;
  std::vector<Eigen::Triplet<T>> triplets;

  // Read fields from comment, if present
  int cols = -1;
  int rows = -1;

  // Also keep track of max triplet entries in case comment is absent
  int maxCol = -1;
  int maxRow = -1;

  while (getline(in, line)) {
    std::stringstream ss(line);
    std::string token;
    ss >> token;
    if (token == "#") {
      // Parse comment of the form (# sparse rows cols)
      ss >> token >> rows >> cols;
    } else {
      // Read triplet
      int iRow = std::stoi(token);
      int iCol;
      T val;
      ss >> iCol >> val;
      // Shift from 1-indexed to 0-indexed
      iRow -= 1;
      iCol -= 1;
      triplets.emplace_back(iRow, iCol, val);
      maxCol = std::max(maxCol, iCol);
      maxRow = std::max(maxRow, iRow);
    }
  }

  if (cols > 0 && (maxCol >= cols || maxRow >= rows)) {
    std::cout << "Error: stated cols (" << cols << ") and rows (" << rows << ") are insufficient to fit all entries."
              << std::endl;
    std::cout << "       triplets are present with iRow=" << maxRow << " or iCol=" << maxCol << std::endl;
    throw std::runtime_error("failed to parse matrix");
  }
  if (cols < 0) {
    cols = maxCol + 1;
    rows = maxRow + 1;
  }

  SparseMatrix<T> M(rows, cols);
  M.setFromTriplets(triplets.begin(), triplets.end());
  return M;
}

template <typename T>
DenseMatrix<T> loadDenseMatrix(std::string filename) {
  std::ifstream inFile(filename);
  if (!inFile) {
    throw std::runtime_error("failed to open input file " + filename);
  }

  DenseMatrix<T> M;
  try {
    M = loadDenseMatrix<T>(inFile);
  } catch (const std::runtime_error& e) {
    throw std::runtime_error("failed to parse dense matrix " + filename);
  }
  inFile.close();
  return M;
}

template <typename T>
DenseMatrix<T> loadDenseMatrix(std::istream& in) {
  std::string line;
  // Parse first line separately - it should be a comment of the form # dense rows cols
  getline(in, line);
  std::stringstream ss(line);
  std::string token, matrixType;
  int rows, cols;
  ss >> token >> matrixType >> rows >> cols;
  if (token != "#" || matrixType != "dense" || rows < 0 || cols < 0) {
    std::cout << "Error: first line of file must be of the form # dense cols rows." << std::endl;
    std::cout << "       instead, it was " << line << std::endl;
    throw std::runtime_error("failed to parse matrix");
  }

  DenseMatrix<T> M(rows, cols);

  int iRow = 0;
  // Parse matrix entries
  while (getline(in, line)) {
    std::stringstream ss(line);
    // Read entries
    T val;
    for (int iCol = 0; iCol < cols; iCol++) {
      ss >> val;
      M(iRow, iCol) = val;
    }
    iRow++;
  }

  return M;
}
