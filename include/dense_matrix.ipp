#include <cassert>
#include <cmath>

#include <algorithm>
#include <iostream>
using namespace std;

namespace geometrycentral {
#ifdef HAVE_SUITESPARSE
extern LinearContext context;
#endif

template <typename T>
DenseMatrix<T>::DenseMatrix(size_t m_, size_t n_)
    // initialize an mxn matrix
    : m(m_),
      n(n_) {
#ifdef HAVE_SUITESPARSE
  cData = NULL;
#endif
  data.resize(m * n);
  zero();
}

template <typename T>
DenseMatrix<T>::DenseMatrix(const DenseMatrix<T>& A)
// copy constructor
{
#ifdef HAVE_SUITESPARSE
  cData = NULL;
#endif
  *this = A;
}

template <typename T>
DenseMatrix<T>::~DenseMatrix(void)
// destructor
{
#ifdef HAVE_SUITESPARSE
  if (cData) {
    cholmod_l_free_dense(&cData, context);
  }
#endif
}

template <typename T>
void DenseMatrix<T>::resize(size_t m_, size_t n_) {
#ifdef HAVE_SUITESPARSE
  if (cData) {
    cholmod_l_free_dense(&cData, context);
    cData = NULL;
  }
#endif

  m = m_;
  n = n_;
  data.resize(m * n);
  zero();
}

template <typename T>
DenseMatrix<T> DenseMatrix<T>::transpose(void) const {
  const DenseMatrix<T>& A(*this);
  DenseMatrix<T> AT(n, m);

  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < m; j++) {
      AT(i, j) = A(j, i);
    }

  return AT;
}

template <typename T>
SparseMatrix<T> DenseMatrix<T>::toSparse(void)
// converts to a sparse matrix
{
  SparseMatrix<T> B(n, m);

#ifdef HAVE_SUITESPARSE
  B = cholmod_l_dense_to_sparse(this->to_cholmod(), true, context);

#else
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < m; j++) {
      B(i, j) = (*this)(i, j);
    }

#endif
  return B;
}

template <typename T>
DenseMatrix<T> DenseMatrix<T>::operator*(const DenseMatrix<T>& B) const
// returns product of this matrix with B
{
  const DenseMatrix<T>& A(*this);

#ifndef NDEBUG
  // make sure matrix dimensions agree
  assert(A.nColumns() == B.nRows());
#endif

  DenseMatrix<T> AB(A.nRows(), B.nColumns());

  for (size_t i = 0; i < A.nRows(); i++)
    for (size_t j = 0; j < B.nColumns(); j++)
      for (size_t k = 0; k < A.nColumns(); k++) {
        AB(i, j) += A(i, k) * B(k, j);
      }

  return AB;
}

template <typename T>
void DenseMatrix<T>::operator*=(const T& c) {
  DenseMatrix<T>& A(*this);

  for (size_t i = 0; i < m; i++)
    for (size_t j = 0; j < n; j++) {
      A(i, j) *= c;
    }
}

template <typename T>
void DenseMatrix<T>::operator/=(const T& c) {
  DenseMatrix<T>& A(*this);

  for (size_t i = 0; i < m; i++)
    for (size_t j = 0; j < n; j++) {
      A(i, j) /= c;
    }
}

template <typename T>
DenseMatrix<T> DenseMatrix<T>::operator+(const DenseMatrix<T>& B) const
// returns sum of this matrix with B
{
  const DenseMatrix& A(*this);

#ifndef NDEBUG
  // make sure matrix dimensions agree
  assert(A.nRows() == B.nRows());
  assert(A.nColumns() == B.nColumns());
#endif

  DenseMatrix<T> C(nRows(), nColumns());

  for (size_t i = 0; i < nRows(); i++)
    for (size_t j = 0; j < nColumns(); j++) {
      C(i, j) = A(i, j) + B(i, j);
    }

  return C;
}

template <typename T>
void DenseMatrix<T>::operator+=(const DenseMatrix<T>& B) {
  DenseMatrix<T>& A(*this);

#ifndef NDEBUG
  // make sure matrix dimensions agree
  assert(A.nRows() == B.nRows());
  assert(A.nColumns() == B.nColumns());
#endif

  for (size_t i = 0; i < nRows(); i++)
    for (size_t j = 0; j < nColumns(); j++) {
      A(i, j) += B(i, j);
    }
}

template <typename T>
DenseMatrix<T> DenseMatrix<T>::operator-(const DenseMatrix<T>& B) const
// returns difference of this matrix with B
{
  const DenseMatrix<T>& A(*this);

#ifndef NDEBUG
  // make sure matrix dimensions agree
  assert(A.nRows() == B.nRows());
  assert(A.nColumns() == B.nColumns());
#endif

  DenseMatrix C(nRows(), nColumns());

  for (size_t i = 0; i < nRows(); i++)
    for (size_t j = 0; j < nColumns(); j++) {
      C(i, j) = A(i, j) - B(i, j);
    }

  return C;
}

template <typename T>
void DenseMatrix<T>::operator-=(const DenseMatrix<T>& B) {
  DenseMatrix<T>& A(*this);

#ifndef NDEBUG
  // make sure matrix dimensions agree
  assert(A.nRows() == B.nRows());
  assert(A.nColumns() == B.nColumns());
#endif

  for (size_t i = 0; i < nRows(); i++)
    for (size_t j = 0; j < nColumns(); j++) {
      A(i, j) -= B(i, j);
    }
}

template <typename T>
DenseMatrix<T> operator*(const T& c, const DenseMatrix<T>& A) {
  DenseMatrix<T> cA = A;

  cA *= c;

  return cA;
}

template <typename T>
DenseMatrix<T> operator*(const DenseMatrix<T>& A, const T& c) {
  return c * A;
}

template <typename T>
DenseMatrix<T> operator/(const DenseMatrix<T>& A, const T& c) {
  DenseMatrix<T> Ac = A;

  Ac /= c;

  return Ac;
}

template <typename T>
const DenseMatrix<T>& DenseMatrix<T>::operator=(const DenseMatrix<T>& B)
// copies B
{
#ifdef HAVE_SUITESPARSE
  if (cData) {
    cholmod_l_free_dense(&cData, context);
    cData = NULL;
  }
#endif

  m = B.m;
  n = B.n;
  data = B.data;

  return *this;
}

template <typename T>
const DenseMatrix<T>& DenseMatrix<T>::operator=(
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat) {
#ifdef HAVE_SUITESPARSE
  if (cData) {
    cholmod_l_free_dense(&cData, context);
    cData = NULL;
  }
#endif

  resize(mat.rows(), mat.cols());

  for (size_t i = 0; i < m; i++) {
    for (size_t j = 0; j < n; j++) {
      (*this)(i, j) = mat(i, j);
    }
  }

  return *this;
}

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> DenseMatrix<T>::toEigen(
    void) const {
  const DenseMatrix<T>& A(*this);
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> B(m, n);

  for (size_t i = 0; i < m; i++)
    for (size_t j = 0; j < n; j++) {
      B(i, j) = A(i, j);
    }

  return B;
}

template <typename T>
SparseMatrix<T> DenseMatrix<T>::asDiagonal() const {
  size_t N = nRows();
  SparseMatrix<T> result(N, N);

  for (size_t i = 0; i < N; i++) {
    result(i, i) = (*this)(i);
  }

  return result;
}

template <typename T>
DenseMatrix<T> DenseMatrix<T>::sub(size_t r0, size_t r1) {
#ifndef NDEBUG
  // check bounds
  assert(0 <= r0 && r0 <= r1 && r1 < m);
#endif

  size_t M = r1 - r0 + 1;
  DenseMatrix<T> b(M);

  for (size_t i = r0; i <= r1; i++) {
    b(i - r0) = data[i];
  }

  return b;
}

template <typename T>
void DenseMatrix<T>::setSubMatrix(size_t r0, size_t r1, size_t c0, size_t c1,
                                  const DenseMatrix<T>& B)
// replaces entries in the inclusive range (r0,r1) x (c0,c1) with the matrix B
{
#ifndef NDEBUG
  // check bounds
  assert(0 <= r0 && r0 <= r1 && r1 < m);
  assert(0 <= c0 && c0 <= c1 && c1 < n);
  assert(r1 - r0 + 1 == B.m);
  assert(c1 - c0 + 1 == B.n);
#endif

  // copy new entries
  DenseMatrix<T>& A(*this);
  for (size_t r = r0; r <= r1; r++)
    for (size_t c = c0; c <= c1; c++) {
      A(r, c) = B(r - r0, c - c0);
    }
}

template <typename T>
void DenseMatrix<T>::increment(size_t r0, const Vector3& u)
// adds the 3-vector u to the 3x1 block starting at r0
{
#ifndef NDEBUG
  assert(0 <= r0 && r0 <= m - 3);
#endif

  data[r0 + 0] = u.x;
  data[r0 + 1] = u.y;
  data[r0 + 2] = u.z;
}

template <typename T>
size_t DenseMatrix<T>::nRows(void) const
// returns the number of rows
{
  return m;
}

template <typename T>
size_t DenseMatrix<T>::nColumns(void) const
// returns the number of columns
{
  return n;
}

template <typename T>
size_t DenseMatrix<T>::length(void) const
// returns the size of the largest dimension
{
  return max(m, n);
}

template <typename T>
void DenseMatrix<T>::zero(const T& val)
// sets all elements to val
{
  for (size_t i = 0; i < m * n; i++) {
    data[i] = val;
  }
}

template <typename T>
double DenseMatrix<T>::norm(NormType type) const {
  double r = 0.;

  if (type == lInfinity) {
    for (size_t i = 0; i < m * n; i++) {
      r = max(r, sqrt(std::norm(data[i])));
    }
  } else if (type == lOne) {
    for (size_t i = 0; i < m * n; i++) {
      r += sqrt(std::norm(data[i]));
    }
  } else if (type == lTwo) {
    for (size_t i = 0; i < m * n; i++) {
      r += std::norm(data[i]);
    }
    r = sqrt(r);
  }

  return r;
}

template <typename T>
void DenseMatrix<T>::normalize(void)
// divides by l2 norm
{
  *this /= norm(lTwo);
}

template <typename T>
T& DenseMatrix<T>::operator()(size_t row, size_t col) {
  assert(row < m && col < n);
  return data[row + m * col];
}

template <typename T>
T DenseMatrix<T>::operator()(size_t row, size_t col) const {
  assert(row < m && col < n);
  return data[row + m * col];
}

template <typename T>
T& DenseMatrix<T>::operator()(size_t index) {
  return data[index];
}

template <typename T>
T DenseMatrix<T>::operator()(size_t index) const {
  return data[index];
}

template <typename T>
T DenseMatrix<T>::sum(void) const
// returns the sum of all entries
{
  T total(0.0);

  for (size_t i = 0; i < m * n; i++) {
    total += data[i];
  }

  return total;
}

template <typename T>
void DenseMatrix<T>::removeMean(void) {
  T mean = 0.;
  size_t N = m * n;

  for (size_t i = 0; i < N; i++) {
    mean += data[i];
  }

  mean /= (double)N;

  for (size_t i = 0; i < N; i++) {
    data[i] -= mean;
  }
}

template <typename T>
double dot(const DenseMatrix<T>& x, const DenseMatrix<T>& y)
// returns Euclidean inner product of x and y
{
  return (x.transpose() * y)(0);
}

template <typename T>
DenseMatrix<T> DenseMatrix<T>::operator-(void) const
// returns additive inverse of this matrix
{
  const DenseMatrix<T>& A(*this);
  DenseMatrix<T> B(m, n);

  for (size_t i = 0; i < m; i++)
    for (size_t j = 0; j < n; j++) {
      B(i, j) = -A(i, j);
    }

  return B;
}

template <typename T>
T inner(const DenseMatrix<T>& x, const DenseMatrix<T>& y)
// standard inner product
{
  T sum = 0.;

#ifndef NDEBUG
  assert(x.nRows() == y.nRows() && x.nColumns() == y.nColumns());
#endif

  for (size_t i = 0; i < x.nRows() * x.nColumns(); i++) {
    sum += conj(x(i)) * y(i);
  }

  return sum;
}

template <typename T>
T inner(const DenseMatrix<T>& x, const DenseMatrix<T>& B,
        const DenseMatrix<T>& y)
// inner product with respect a diagonal inner
// product B represented as a dense vector
{
  T sum = 0.;

#ifndef NDEBUG
  assert(x.nRows() == y.nRows() && x.nRows() == B.nRows() &&
         x.nColumns() == 1 && B.nColumns() == 1 && y.nColumns() == 1);
#endif

  for (size_t i = 0; i < x.nRows() * x.nColumns(); i++) {
    sum += conj(x(i)) * B(i) * y(i);
  }

  return sum;
}
}
