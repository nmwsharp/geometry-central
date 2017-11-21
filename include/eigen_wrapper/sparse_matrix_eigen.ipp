#include <algorithm>
#include <cassert>
#include <iostream>
#include <cmath>

namespace geometrycentral {

// initialize an mxn matrix
template <typename T>
SparseMatrix<T>::SparseMatrix(size_t m_, size_t n_)
    : factors(*this), m(m_), n(n_) {}

// copy constructor
template <typename T>
SparseMatrix<T>::SparseMatrix(const SparseMatrix<T>& B)
    : factors(*this) {
  *this = B;
}

// destructor
template <typename T>
SparseMatrix<T>::~SparseMatrix(void) {}

// copies B
template <typename T>
const SparseMatrix<T>& SparseMatrix<T>::operator=(const SparseMatrix<T>& B) {
  // if the new matrix has the same sparsity as
  // the old one, keep the current symbolic
  // factorizations and just rebuild the numeric
  // ones as necessary
  if (sameSparsity(*this, B)) {
    factors.clearNumeric();
  } else {
    // otherwise, get rid of the current factors
    factors.clear();
  }

  m = B.m;
  n = B.n;
  data = B.data;

  return *this;
}

template <typename T>
Eigen::SparseMatrix<T> SparseMatrix<T>::toEigen(void) {
  // Convert to Eigen format by converting our std::map from indices to values
  // into a flat list of triplets, then calling
  // Eigen::SparseMatrix<T>::setFromTriplets()
  // to build the final Eigen matrix.

  // TODO Since setFromTriplets() takes a begin() and end() iterator,
  // we could instead just iterate over the map with an alternative
  // iterator that dereferences to an Eigen::Triplet<T>.  This way
  // we would save both time and memory.

  typedef Eigen::Triplet<T> Triplet;
  std::vector<Triplet> entries(data.size());

  int k = 0;
  for (auto e = begin(); e != end(); e++) {
    entries[k] = Triplet(e->first.row, e->first.col, e->second);
    k++;
  }

  Eigen::SparseMatrix<T> A(m, n);
  A.setFromTriplets(entries.begin(), entries.end());

  return A;
}

// returns product of this matrix with sparse B
template <typename T>
SparseMatrix<T> SparseMatrix<T>::operator*(const SparseMatrix<T>& B) const {
  const SparseMatrix<T>& A(*this);

  // make sure matrix dimensions agree
  assert(A.nColumns() == B.nRows());

  // collect nonzeros in each row
  vector<vector<size_t> > Bcol(B.nRows());
  vector<vector<T> > Bval(B.nRows());
  for (const_iterator e = B.begin(); e != B.end(); e++) {
    size_t row = e->first.row;
    size_t col = e->first.col;
    T val = e->second;

    Bcol[row].push_back(col);
    Bval[row].push_back(val);
  }

  // multiply C = A*B
  SparseMatrix<T> C(A.nRows(), B.nColumns());
  for (const_iterator e = begin(); e != end(); e++) {
    size_t i = e->first.row;
    size_t j = e->first.col;

    for (size_t n = 0; n < Bcol[j].size(); n++) {
      size_t k = Bcol[j][n];

      C(i, k) += e->second * Bval[j][n];
    }
  }

  return C;
}

// returns product of this matrix with dense B
template <typename T>
DenseMatrix<T> SparseMatrix<T>::operator*(const DenseMatrix<T>& B) const {
  const SparseMatrix<T>& A(*this);

  // make sure matrix dimensions agree
  assert(A.nColumns() == B.nRows());

  // multiply C = A*B
  DenseMatrix<T> C(A.nRows(), B.nColumns());
  for (const_iterator e = begin(); e != end(); e++) {
    size_t i = e->first.row;
    size_t j = e->first.col;

    for (size_t k = 0; k < B.nColumns(); k++) {
      C(i, k) += e->second * B(j, k);
    }
  }

  return C;
}

template <typename T>
void SparseMatrix<T>::operator*=(const T& c) {
  for (iterator e = begin(); e != end(); e++) {
    e->second *= c;
  }

  factors.clearNumeric();
}

template <typename T>
void SparseMatrix<T>::operator/=(const T& c) {
  for (iterator e = begin(); e != end(); e++) {
    e->second /= c;
  }

  factors.clearNumeric();
}

// adds B to this matrix
template <typename T>
void SparseMatrix<T>::operator+=(const SparseMatrix<T>& B) {
  SparseMatrix<T>& A(*this);

  // make sure matrix dimensions agree
  assert(A.nRows() == B.nRows());
  assert(A.nColumns() == B.nColumns());

  for (const_iterator e = B.begin(); e != B.end(); e++) {
    size_t i = e->first.row;
    size_t j = e->first.col;
    const T& Bij(e->second);

    A(i, j) += Bij;
  }

  if (sameSparsity(*this, B)) {
    factors.clearNumeric();
  } else {
    factors.clear();
  }
}

// subtracts B from this matrix
template <typename T>
void SparseMatrix<T>::operator-=(const SparseMatrix<T>& B) {
  SparseMatrix<T>& A(*this);

  // make sure matrix dimensions agree
  assert(A.nRows() == B.nRows());
  assert(A.nColumns() == B.nColumns());

  for (const_iterator e = B.begin(); e != B.end(); e++) {
    size_t i = e->first.row;
    size_t j = e->first.col;
    const T& Bij(e->second);

    A(i, j) -= Bij;
  }

  if (sameSparsity(*this, B)) {
    factors.clearNumeric();
  } else {
    factors.clear();
  }
}

// returns sum of this matrix with B
template <typename T>
SparseMatrix<T> SparseMatrix<T>::operator+(const SparseMatrix<T>& B) const {
  SparseMatrix<T> C(nRows(), nColumns());

  C += *this;
  C += B;

  return C;
}

// returns sum of this matrix with B
template <typename T>
SparseMatrix<T> SparseMatrix<T>::operator-(const SparseMatrix<T>& B) const {
  SparseMatrix<T> C(nRows(), nColumns());

  C += *this;
  C -= B;

  return C;
}

template <typename T>
SparseMatrix<T> operator*(const T& c, const SparseMatrix<T>& A) {
  SparseMatrix<T> cA = A;

  for (typename SparseMatrix<T>::iterator e = cA.begin(); e != cA.end(); e++) {
    e->second = c * e->second;
  }

  return cA;
}

template <typename T>
SparseMatrix<T> operator*(const SparseMatrix<T>& A, const T& c) {
  SparseMatrix<T> Ac = A;

  Ac *= c;

  return Ac;
}

template <typename T>
SparseMatrix<T> operator/(const SparseMatrix<T>& A, T c) {
  SparseMatrix<T> Ac = A;

  Ac /= c;

  return Ac;
}

template <typename T>
void SparseMatrix<T>::resize(size_t m_, size_t n_) {
  m = m_;
  n = n_;

  data.clear();

  factors.clear();
}

template <typename T>
SparseMatrix<T> SparseMatrix<T>::transpose(void) const {
  SparseMatrix<T> AT(n, m);

  for (const_iterator e = begin(); e != end(); e++) {
    size_t i = e->first.row;
    size_t j = e->first.col;
    T Aij = e->second;

    AT(j, i) = conjugate(Aij);
  }

  return AT;
}

// returns the number of rows
template <typename T>
size_t SparseMatrix<T>::nRows(void) const {
  return m;
}

// returns the number of columns
template <typename T>
size_t SparseMatrix<T>::nColumns(void) const {
  return n;
}

// returns the size of the largest dimension
template <typename T>
size_t SparseMatrix<T>::length(void) const {
  return max(m, n);
}

// sets all nonzero elements to val
template <typename T>
void SparseMatrix<T>::zero(const T& val) {
  for (iterator i = begin(); i != end(); i++) {
    i->second = val;
  }

  factors.clearNumeric();
}

template <typename T>
size_t SparseMatrix<T>::nNonZeros(void) const {
  return data.size();
}

template <typename T>
void SparseMatrix<T>::invertDiagonal(void) {
  SparseMatrix<T>& A(*this);

  for (size_t i = 0; i < max(m, n); i++) {
    A(i, i) = A(i, i).inv();
  }

  factors.clearNumeric();
}

template <typename T>
SparseMatrix<T> SparseMatrix<T>::identity(size_t N) {
  SparseMatrix<T> I(N, N);

  for (size_t i = 0; i < N; i++) {
    I(i, i) = 1.;
  }

  return I;
}

template <typename T>
DenseMatrix<T> SparseMatrix<T>::toDense(void) const
// converts to a dense matrix
{
  const size_t maxSize = 1048576;
  if (m * n > maxSize) {
    cerr << "Error: refusing to convert sparse to dense (too big!)"
         << "\n";
    exit(1);
  }

  const SparseMatrix<T>& A(*this);
  DenseMatrix<T> B(m, n);

  for (size_t i = 0; i < m; i++)
    for (size_t j = 0; j < n; j++) {
      B(i, j) = A(i, j);
    }

  return B;
}

template <typename T>
T& SparseMatrix<T>::operator()(size_t row, size_t col) {
  assert(row >= 0);
  assert(row < m);
  assert(col >= 0);
  assert(col < n);

  EntryIndex index{row, col};
  const_iterator entry = data.find(index);

  if (entry == end()) {
    // data[ index ] = T( 0. );
    data.insert(pair<EntryIndex, T>(index, T(0.)));

    // creating a nonzero invalidates both the
    // current symbolic and numeric factorization(s)
    factors.clear();
  } else {
    // non-const access potentially invalidates
    // the current numeric factorization(s)
    factors.clearNumeric();
  }

  return data[index];
}

template <typename T>
T SparseMatrix<T>::operator()(size_t row, size_t col) const {
  assert(row >= 0);
  assert(row < m);
  assert(col >= 0);
  assert(col < n);

  EntryIndex index{row, col};
  const_iterator entry = data.find(index);

  if (entry == end()) {
    return T(0.);
  }

  return entry->second;
}

template <typename T>
typename SparseMatrix<T>::iterator SparseMatrix<T>::begin(void) {
  return data.begin();
}

template <typename T>
typename SparseMatrix<T>::const_iterator SparseMatrix<T>::begin(void) const {
  return data.begin();
}

template <typename T>
typename SparseMatrix<T>::iterator SparseMatrix<T>::end(void) {
  return data.end();
}

template <typename T>
typename SparseMatrix<T>::const_iterator SparseMatrix<T>::end(void) const {
  return data.end();
}

template <typename T>
void SparseMatrix<T>::shift(double c)
// adds c times the identity matrix to this matrix
{
  assert(m == n);
  SparseMatrix<T>& A(*this);

  for (size_t i = 0; i < m; i++) {
    A(i, i) += c;
  }
}

template <typename T>
SparseMatrix<T> SparseMatrix<T>::sub(size_t r0, size_t r1, size_t c0,
                                     size_t c1) const {
  // check bounds
  assert(0 <= r0 && r0 <= r1 && r1 < m);
  assert(0 <= c0 && c0 <= c1 && c1 < n);

  size_t M = r1 - r0 + 1;
  size_t N = c1 - c0 + 1;
  SparseMatrix<T> B(M, N);

  for (const_iterator e = begin(); e != end(); e++) {
    size_t i = e->first.row;
    size_t j = e->first.col;

    if (r0 <= i && i <= r1 && c0 <= j && j <= c1) {
      B(i - r0, j - c0) = e->second;
    }
  }

  return B;
}

// replaces entries in the inclusive range (r0,r1) x (c0,c1) with the matrix B
template <typename T>
void SparseMatrix<T>::setSubMatrix(size_t r0, size_t r1, size_t c0, size_t c1,
                                   const SparseMatrix<T>& B) {
  // check bounds
  assert(0 <= r0 && r0 <= r1 && r1 < m);
  assert(0 <= c0 && c0 <= c1 && c1 < n);
  assert(r1 - r0 + 1 == B.m);
  assert(c1 - c0 + 1 == B.n);

  // erase any nonzero entries within this block
  iterator e = begin(), eNext;
  while (e != end()) {
    eNext = e;
    eNext++;

    size_t i = e->first.row;
    size_t j = e->first.col;
    if (r0 <= i && i <= r1 && c0 <= j && j <= c1) {
      data.erase(e);
    }

    e = eNext;
  }

  // copy new entries
  SparseMatrix<T>& A(*this);
  for (const_iterator e = B.begin(); e != B.end(); e++) {
    size_t i = e->first.row;
    size_t j = e->first.col;

    A(i + r0, j + c0) = e->second;
  }

  factors.clear();
}

// adds the 3x3 matrix B to the 3x3 block starting at (r0,c0)
template <typename T>
void SparseMatrix<T>::increment(size_t r0, size_t c0, const Matrix3x3& B) {
  SparseMatrix<T>& A(*this);

  A(r0 + 0, c0 + 0) += B(0, 0);
  A(r0 + 0, c0 + 1) += B(0, 1);
  A(r0 + 0, c0 + 2) += B(0, 2);
  A(r0 + 1, c0 + 0) += B(1, 0);
  A(r0 + 1, c0 + 1) += B(1, 1);
  A(r0 + 1, c0 + 2) += B(1, 2);
  A(r0 + 2, c0 + 0) += B(2, 0);
  A(r0 + 2, c0 + 1) += B(2, 1);
  A(r0 + 2, c0 + 2) += B(2, 2);
}

// returns true if the matrix contains any INFs or NaNs
template <typename T>
bool SparseMatrix<T>::hasBadEntries(void) const {
  for (const_iterator e = begin(); e != end(); e++) {
    if (!std::isfinite(e->second)) {
      return true;
    }
  }
  return false;
}

template <typename T>
double SparseMatrix<T>::norm(void) const {
  double M = 0.;

  for (const_iterator e = begin(); e != end(); e++) {
    M = max(M, std::norm(e->second));
  }

  return sqrt(M);
}

// returns the max residual of the linear problem A x = b relative to the
// largest entry of the solution
template <typename T>
double residual(const SparseMatrix<T>& A, const DenseMatrix<T>& x,
                const DenseMatrix<T>& b) {
  return (A * x - b).norm(lInfinity) / b.norm(lInfinity);
}

// returns the max residual of the eigenvalue problem A x = lambda x relative to
// the largest entry of the solution
template <typename T>
double residual(const SparseMatrix<T>& A, const DenseMatrix<T>& x) {
  T lambda = rayleighQuotient(A, x);
  return (A * x - lambda * x).norm() / x.norm();
}

// returns the max residual of the generalized eigenvalue problem A x = lambda B
// x relative to the largest entry of the solution
template <typename T>
double residual(const SparseMatrix<T>& A, const SparseMatrix<T>& B,
                const DenseMatrix<T>& x) {
  T lambda = rayleighQuotient(A, B, x);
  return (A * x - lambda * B * x).norm() / x.norm();
}

// returns the max residual of the generalized eigenvalue problem A x = lambda
// (B - EE^T) x relative to the largest entry of the solution
template <class T>
double residual(const SparseMatrix<T>& A, const SparseMatrix<T>& B,
                const DenseMatrix<T>& E, const DenseMatrix<T>& x) {
  T lambda = rayleighQuotient(A, B, E, x);
  return (A * x - lambda * (B * x - E * (E.transpose() * x))).norm() / x.norm();
}

// returns <x,Ax>/<x,x>
template <typename T>
T rayleighQuotient(const SparseMatrix<T>& A, const DenseMatrix<T>& x) {
  DenseMatrix<T> xT = x.transpose();
  return (xT * (A * x))(0) / (xT * x)(0);
}

// returns <Ax,x>/<Bx,x>
template <typename T>
T rayleighQuotient(const SparseMatrix<T>& A, const SparseMatrix<T>& B,
                   const DenseMatrix<T>& x) {
  DenseMatrix<T> xT = x.transpose();
  return (xT * (A * x))(0) / (xT * (B * x))(0);
}

// returns <Ax,x>/<(B-EE^T)x,x>
template <class T>
T rayleighQuotient(const SparseMatrix<T>& A, const SparseMatrix<T>& B,
                   const DenseMatrix<T>& E, const DenseMatrix<T>& x) {
  return (x.transpose() * (A * x))(0) /
         (x.transpose() * (B * x - E * (E.transpose() * x)))(0);
}

template <class T, class U>
bool sameSparsity(const SparseMatrix<T>& A, const SparseMatrix<U>& B) {
  if (A.nRows() != B.nRows()) return false;
  if (A.nColumns() != B.nColumns()) return false;

  typename SparseMatrix<T>::const_iterator eA, endA;
  typename SparseMatrix<U>::const_iterator eB, endB;

  eA = A.begin();
  endA = A.end();
  eB = B.begin();
  endB = B.end();

  while (eA != endA && eB != endB) {
    if (eA->first.row != eB->first.row || eA->first.col != eB->first.col) {
      return false;
    }

    eA++;
    eB++;
  }

  // TODO We could create a hash of the
  // sparsity pattern every time the symbolic
  // factorization is performed.  However, we still
  // need to iterate over all nonzeros to construct
  // the hash, so it may not be a big win unless we
  // are doing a huge number of matrix-level operations.
  //
  // Hash-on-compress makes sense, since that's what
  // we're caching (and we have to iterate over all
  // entries at that point anyway).

  return true;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const SparseMatrix<T>& o) {
  os.precision(3);

  for (typename SparseMatrix<T>::const_iterator e = o.begin(); e != o.end();
       e++) {
    size_t row = e->first.row;
    size_t col = e->first.col;

    os << "( " << row << ", " << col << " ): " << e->second << "\n";
  }

  return os;
}
}
