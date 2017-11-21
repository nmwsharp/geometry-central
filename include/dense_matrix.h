#pragma once

// DenseMatrix represents an m by n matrix where every entry -- including
// zero-valued entries -- is stored explicitly.
//
// --THIS CLASS SHOULD *NOT* BE USED FOR LARGE, SPARSE MATRICES!--
//             --PERFORMANCE **WILL** BE HORRIBLE---
//
// For matrices with many zeros, see geometrycentral::SparseMatrix.
// geometrycentral::DenseMatrix is
// most commonly used to represent dense vectors in sparse linear systems
// (i.e., the right hand side and the solution vector).
//
// A real or complex matrix is allocated via
//
//    DenseMatrix<entryType> A( nRows, nColumns );
//
// Matrix elements are then accessed using parenthesis, e.g.,
//
//    A(i,j) = 1;
//    A(i,j) += 2;
//    a = A(i,j);
//
// etc.  geometrycentral::SparseMatrix is interoperable with the Eigen numerical
// linear
// algebra library.  In particular, the method
// geometrycentral::SparseMatrix::toEigen
// returns an Eigen::SparseMatrix which can be used by routines in Eigen.  For
// basic operations, however, you should not need to perform this conversion
// explicitly.
//

#include <vector3.h>
#include <Eigen/Core>
#include <vector>

#ifdef HAVE_SUITESPARSE
#include <suitesparse_wrapper/linear_context.h>
#endif

// For real entry types we need to define a conjugation
// operation (which does nothing) in order to implement
// generic methods that call conj()
inline double conjugate(double x) { return x; }

namespace geometrycentral {

template <typename T>
class SparseMatrix;

enum NormType { lInfinity, lOne, lTwo };

template <class T>
class DenseMatrix {
 public:
  DenseMatrix(size_t m = 0, size_t n = 1);
  // initialize an mxn matrix (specifying just m yields a column vector)

  DenseMatrix(const DenseMatrix<T>& A);
  // copy constructor

  const DenseMatrix<T>& operator=(const DenseMatrix<T>& B);
  // copies B

  ~DenseMatrix(void);
  // destructor

  void resize(size_t m = 0, size_t n = 1);
  // resize to an mxn matrix (all previous values are set to zero)

  SparseMatrix<T> toSparse(void);
  // converts to a sparse matrix

  size_t nRows(void) const;
  // returns the number of rows

  size_t nColumns(void) const;
  // returns the number of columns

  size_t length(void) const;
  // returns the size of the largest dimension

  void zero(const T& val = 0.);
  // sets all elements to val

  double norm(NormType type = lInfinity) const;
  // returns the maximum magnitude of any entry

  T& operator()(size_t row, size_t col);
  T operator()(size_t row, size_t col) const;
  // access the specified element of the matrix (uses 0-based indexing)

  T& operator()(size_t index);
  T operator()(size_t index) const;
  // access the specified element of a vector (uses 0-based indexing)

  DenseMatrix<T> transpose(void) const;
  // returns the transpose of this matrix

  DenseMatrix<T> operator*(const DenseMatrix<T>& B) const;
  // returns product of this matrix with B

  void operator*=(const T& c);
  // multiplies this matrix by the scalar c

  void operator/=(const T& c);
  // divides this matrix by the scalar c

  DenseMatrix operator+(const DenseMatrix& B) const;
  // returns sum of this matrix with B

  void operator+=(const DenseMatrix& B);
  // adds B to this matrix

  DenseMatrix operator-(const DenseMatrix& B) const;
  // returns difference of this matrix with B

  void operator-=(const DenseMatrix& B);
  // subtracts B from this matrix

  DenseMatrix<T> operator-(void) const;
  // returns additive inverse of this matrix

  SparseMatrix<T> asDiagonal() const;
  // returns the N x N matrix with this COLUMN VECTOR on the diagonal

  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> toEigen(void) const;
  // returns copy of the matrix in Eigen format, performing
  // conversions as necessary (e.g., for quaternionic matrices)

  const DenseMatrix<T>& operator=(
      const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat);
// resizes and fills this matrix from the specified Eigen matrix

#ifdef HAVE_SUITESPARSE
  cholmod_dense* to_cholmod(void);
  // returns pointer to copy of matrix in CHOLMOD format

  const DenseMatrix<T>& operator=(cholmod_dense* B);
// copies a cholmod_dense* into a DenseMatrix;
// takes responsibility for deallocating B
#endif

  DenseMatrix<T> sub(size_t r0, size_t r1);
  // for a dense vector, returns [r0,r1]

  void setSubMatrix(size_t r0, size_t r1, size_t c0, size_t c1,
                    const DenseMatrix<T>& B);
  // replaces entries in the inclusive range (r0,r1) x (c0,c1) with the matrix B

  void increment(size_t r0, const Vector3& u);
  // adds the 3-vector u to the 3x1 block starting at r0

  void normalize(void);
  // divides real part by Frobenius norm of real part

  T sum(void) const;
  // returns the sum of all entries

  void removeMean(void);
  // removes the mean of the real part from the real part

  void randomize(void);
  // replaces entries with uniformly distributed random real numbers in the
  // interval [-1,1]

 protected:
  size_t m, n;
  std::vector<T> data;

#ifdef HAVE_SUITESPARSE
  cholmod_dense* cData;
#endif
};

template <class T>
DenseMatrix<T> operator*(const DenseMatrix<T>& A, const T& c);
// right scalar multiplication

template <class T>
DenseMatrix<T> operator*(const T& c, const DenseMatrix<T>& A);
// left scalar multiplication

template <class T>
DenseMatrix<T> operator/(const DenseMatrix<T>& A, const T& c);
// scalar division

template <class T>
double dot(const DenseMatrix<T>& x, const DenseMatrix<T>& y);
// returns Euclidean inner product of x and y

template <class T>
std::ostream& operator<<(std::ostream& os, const DenseMatrix<T>& o);
// prints entries

template <class T>
T inner(const DenseMatrix<T>& x, const DenseMatrix<T>& y);
// standard inner product

template <class T>
T inner(const DenseMatrix<T>& x, const DenseMatrix<T>& B,
        const DenseMatrix<T>& y);
// inner product with respect to a diagonal inner
// product B represented as a dense vector

// Syntatic sugar for vector type
// Note that the burden is on the user not to store a general matrix and call it
// a vector.
template <typename T>
using DenseVector = DenseMatrix<T>;

}  // namespace geometrycentral

#include "dense_matrix.ipp"
