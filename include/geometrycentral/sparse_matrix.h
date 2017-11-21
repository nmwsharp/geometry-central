#pragma once

// geometrycentral::SparseMatrix represents an m by n matrix where only nonzero
// entries are
// stored explicitly.  This class is most commonly used to represent the linear
// term in sparse linear systems (i.e., the matrix part).
//
// A matrix is allocated via
//
//     SparseMatrix<entryType> A(nRows, nColumns);
//
// Matrix elements are then accessed using parenthesis, e.g.,
//
//     A(i,j)  = 1;
//     A(i,j) += 2;
//     a = A(i,j);
//
// etc.  geometrycentral::SparseMatrix is interoperable with the Eigen numerical
// linear
// algebra library.  In particular, the method
// geometrycentral::SparseMatrix::toEigen
// returns an Eigen::SparseMatrix which can be used by routines in Eigen.  For
// basic operations, however, you should not need to perform this conversion
// explicitly.
//
// Internally SparseMatrix stores nonzero entries in a heap data structure; the
// amortized cost of insertion is therefore no worse than the sorting cost of
// putting the matrix in compressed-column order.
//

#include <map>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include "geometrycentral/dense_matrix.h"
#include "geometrycentral/quaternion.h"
#include "geometrycentral/timing.h"

#ifdef HAVE_SUITESPARSE
#include "suitesparse_wrapper/linear_context.h"
#include "suitesparse_wrapper/sparse_factorization.h"
#else
#include "eigen_wrapper/sparse_factorization.h"
#endif

namespace geometrycentral {
// XXX dummy class; still need to port from libDDG
template <typename T>
class SparseMatrix;
typedef Eigen::Matrix<double, 3, 3> Matrix3x3;

class EntryIndex {
 public:
  size_t row;
  size_t col;

  // compare by column-major order
  bool operator<(const EntryIndex& e) const {
    if (col < e.col) return true;
    if (col > e.col) return false;
    if (row < e.row) return true;
    return false;
  }
};

template <typename T>
class SparseMatrix {
 public:
  SparseMatrix(size_t m = 0, size_t n = 0);
  // initialize an mxn matrix

  SparseMatrix(const SparseMatrix<T>& B);
  // copy constructor

  ~SparseMatrix(void);
  // destructor

  const SparseMatrix<T>& operator=(const SparseMatrix<T>& B);
  // copies B

  void resize(size_t m, size_t n);
  // clears and resizes to mxn matrix

  SparseMatrix<T> transpose(void) const;
  // returns the transpose of this matrix

  SparseMatrix<double> toReal();
  // converts this matrix to its real representation
  // (currently only defined for quaternionic matrices)

  SparseMatrix<Quaternion> toQuaternionic();
  // converts this matrix to its quaternionic representation
  // (currently only defined for real matrices)

  Eigen::SparseMatrix<T> toEigen(void);
  // returns copy of the matrix in Eigen format, performing
  // conversions as necessary (e.g., for quaternionic matrices)

  SparseMatrix<T> operator*(const SparseMatrix<T>& B) const;
  // returns product of this matrix with sparse B

  DenseMatrix<T> operator*(const DenseMatrix<T>& B) const;
  // returns product of this matrix with dense B

  DenseMatrix<T> multiply(DenseMatrix<T>& x, bool transpose = false);
  // returns product of this matrix with dense x; if
  // transpose is true, uses transpose of this matrix

  SparseMatrix<T> multiply(SparseMatrix<T>& B, bool transposeA = false,
                           bool transposeB = false);
  // returns product of this matrix with sparse B

  void operator*=(const T& c);
  // multiplies this matrix by the scalar c

  void operator/=(const T& c);
  // divides this matrix by the scalar c

  void operator+=(const SparseMatrix<T>& B);
  // adds B to this matrix

  void operator-=(const SparseMatrix<T>& B);
  // subtracts B from this matrix

  SparseMatrix<T> operator+(const SparseMatrix<T>& B) const;
  // returns sum of this matrix with B

  SparseMatrix<T> operator-(const SparseMatrix<T>& B) const;
  // returns difference of this matrix with B

  size_t nRows(void) const;
  // returns the number of rows

  size_t nColumns(void) const;
  // returns the number of columns

  size_t length(void) const;
  // returns the size of the largest dimension

  size_t nNonZeros(void) const;
  // returns the number of nonzeros

  void zero(const T& val);
  // sets all nonzero elements to val

  void invertDiagonal(void);
  // inverts diagonal elements

  static SparseMatrix<T> identity(size_t N);
  // returns the N x N identity matrix

  DenseMatrix<T> toDense(void) const;
  // converts to a dense matrix

  T& operator()(size_t row, size_t col);
  T operator()(size_t row, size_t col) const;
  // access the specified element (uses 0-based indexing)

  typedef std::map<EntryIndex, T> EntryMap;

  typedef typename EntryMap::iterator iterator;
  typedef typename EntryMap::const_iterator const_iterator;
  // convenience types for storing and accessing entries

  iterator begin(void);
  iterator end(void);
  const_iterator begin(void) const;
  const_iterator end(void) const;
  // return iterators to first and last nonzero entries

  void shift(double c);
  // adds c times the identity matrix to this matrix

  SparseMatrix<T> sub(size_t r0, size_t r1, size_t c0, size_t c1) const;
  // returns the sub-block within the specified range [r0,r1]x[c0,c1]

  void setSubMatrix(size_t r0, size_t r1, size_t c0, size_t c1,
                    const SparseMatrix<T>& B);
  // replaces entries in the inclusive range (r0,r1) x (c0,c1) with the matrix B

  void increment(size_t r0, size_t c0, const Matrix3x3& B);
  // adds the 3x3 matrix B to the 3x3 block starting at (r0,c0)

  bool hasBadEntries(void) const;
  // returns true if the matrix contains any INFs or NaNs

  double norm(void) const;
  // returns the maximum magnitude of any entry
  // TODO add option to get Frobenius or other matrix norms

  SparseFactorization<T> factors;
  // Factorizations for various direct solvers (built only
  // if needed). Symbolic part is invalidated by creation
  // of new nonzeros, numeric part is invalidated by non-
  // const access.  Matrix-level operations invalidate both
  // parts. Note that the code for this class comes from
  // one of several implementations based on which matrix
  // backend is being used.

 protected:
  size_t m, n;
  EntryMap data;

#ifdef HAVE_SUITESPARSE

 public:
  const SparseMatrix<T>& operator=(cholmod_sparse* B);
  // copies a cholmod_sparse* into a SparseMatrix;
  // takes responsibility for deallocating B
  cholmod_sparse* to_cholmod(void);
  // returns pointer to copy of matrix in compressed-column CHOLMOD format

 protected:
  cholmod_sparse* cData;

  void clearCompressed(void);
  void buildCompressed(void);

  void allocateSparse(void);
  void setEntry(const_iterator e, int i, double* pr);
#endif
};

template <typename T>
SparseMatrix<T> operator*(const SparseMatrix<T>& A, const T& c);
// right scalar multiplication

template <typename T>
SparseMatrix<T> operator*(const T& c, const SparseMatrix<T>& A);
// left scalar multiplication

template <typename T>
SparseMatrix<T> operator/(const SparseMatrix<T>& A, const T& c);
// scalar division

template <typename T, typename U>
bool sameSparsity(const SparseMatrix<T>& A, const SparseMatrix<U>& B);
// returns true if and only if A have the same size and pattern of nonzeros

template <typename T>
std::ostream& operator<<(std::ostream& os, const SparseMatrix<T>& o);
// prints entries

template <typename T>
void solve(SparseMatrix<T>& A, DenseMatrix<T>& x, DenseMatrix<T> b);
// solves the sparse linear system Ax = b using sparse QR factorization

template <typename T>
void solveSquare(SparseMatrix<T>& A, DenseMatrix<T>& x, DenseMatrix<T> b);
// solves the sparse linear system Ax = b using sparse LU factorization

// THIS IS __SYMMETRIC__ positive definite
template <typename T>
void solvePositiveDefinite(SparseMatrix<T>& A, DenseMatrix<T>& x,
                           DenseMatrix<T> b);
// solves the positive definite sparse linear system Ax = b using sparse
// Cholesky factorization

template <typename T>
void solveIndefinite(SparseMatrix<T>& A, DenseMatrix<T>& x, DenseMatrix<T> b);
// solves the indefinite sparse linear system Ax = b using sparse LDL^T
// factorization

template <typename T>
void solveModifiedIndefinite(SparseMatrix<T>& A, DenseMatrix<T>& x,
                             DenseMatrix<T> b);
// solves the system LD'Lx = b, where LDL = A and D' is a modification
// of the diagonal such that the resulting system is positive-definite

template <typename T>
void solveLeastSquares(SparseMatrix<T>& A, DenseMatrix<T>& x, DenseMatrix<T> b);
// solves the sparse linear system Ax = b, yielding the least-squares
// solution for over- or under-determined systems

template <typename T>
void solveConjugateGradient(SparseMatrix<T>& A, SparseMatrix<T>& M,
                            DenseMatrix<T>& x, DenseMatrix<T>& b,
                            double relativeTolerance = 1e-7,
                            int maxIterations = 100);
// solves the sparse positive semidefinite linear system Ax = b using the
// conjugate gradient method with preconditioner M; x is used as an
// initial guess

template <typename T>
void smallestEig(SparseMatrix<T>& A, DenseMatrix<T>& x,
                 bool ignoreConstantVector = true);
// solves A x = lambda x for the smallest nonzero eigenvalue lambda
// A must be square; x is used as an initial guess

template <typename T>
void smallestEig(SparseMatrix<T>& A, SparseMatrix<T>& B, DenseMatrix<T>& x);
// solves A x = lambda B x for the smallest nonzero generalized eigenvalue
// lambda
// A and B must be symmetric; x is used as an initial guess

template <typename T>
void smallestEigPositiveDefinite(SparseMatrix<T>& A, DenseMatrix<T>& x,
                                 bool ignoreConstantVector = true);
// solves A x = lambda x for the smallest nonzero eigenvalue lambda
// A must be positive (semi-)definite; x is used as an initial guess

template <typename T>
void smallestEigPositiveDefinite(SparseMatrix<T>& A, SparseMatrix<T>& B,
                                 DenseMatrix<T>& x,
                                 bool ignoreConstantVector = true);
// solves A x = lambda B x for the smallest nonzero eigenvalue lambda
// A must be positive (semi-)definite, B must be symmetric; x is used as an
// initial guess

template <typename T>
void smallestEigPositiveDefinite(SparseMatrix<T>& A, SparseMatrix<T>& B,
                                 DenseMatrix<T>& E, DenseMatrix<T>& x);
// solves A x = lambda (B - EE^T) x for the smallest nonzero eigenvalue lambda
// A must be positive (semi-)definite, B must be symmetric; EE^T is a low-rank
// matrix, and
// x is used as an initial guess

template <typename T>
void smallestEigPCG(SparseMatrix<T>& A,
                    SparseMatrix<T>& M,  // preconditioner
                    DenseMatrix<T>& x, bool ignoreConstantVector = true);
// solves A x = lambda x for the smallest nonzero eigenvalue lambda
// A must be positive (semi-)definite; x is used as an initial guess;
// power iterations are computed via preconditioned conjugate gradient

template <typename T>
void solveLeastSquares(const SparseMatrix<T>& A, const DenseMatrix<T>& b,
                       DenseMatrix<T>& x);
// solves the positive definite sparse linear system A^TAx = A^Tb using sparse
// Cholesky factorization

template <typename T>
void smallestEigPositiveDefinite(SparseMatrix<T>& A, SparseMatrix<T>& B,
                                 DenseMatrix<T>& e, DenseMatrix<T>& x);
// solves A x = lambda (B - EE^T) x for the smallest nonzero eigenvalue lambda
// A must be positive (semi-)definite, B must be symmetric; EE^T is a low-rank
// matrix, and
// x is used as an initial guess

template <typename T>
double residual(const SparseMatrix<T>& A, const DenseMatrix<T>& x,
                const DenseMatrix<T>& b);
// returns the max residual of the linear problem A x = b relative to the
// largest entry of the solution

template <typename T>
double residual(const SparseMatrix<T>& A, const DenseMatrix<T>& x);
// returns the max residual of the eigenvalue problem A x = lambda x relative to
// the largest entry of the solution

template <typename T>
double residual(const SparseMatrix<T>& A, const SparseMatrix<T>& B,
                const DenseMatrix<T>& x);
// returns the max residual of the generalized eigenvalue problem A x = lambda B
// x relative to the largest entry of the solution

template <typename T>
double residual(const SparseMatrix<T>& A, const SparseMatrix<T>& B,
                const DenseMatrix<T>& E, const DenseMatrix<T>& x);
// returns the max residual of the generalized eigenvalue problem A x = lambda
// (B - EE^T) x relative to the largest entry of the solution

template <typename T>
T rayleighQuotient(const SparseMatrix<T>& A, const DenseMatrix<T>& x);
// returns <Ax,x>/<x,x>

template <typename T>
T rayleighQuotient(const SparseMatrix<T>& A, const SparseMatrix<T>& B,
                   const DenseMatrix<T>& x);
// returns <Ax,x>/<Bx,x>

template <typename T>
T rayleighQuotient(const SparseMatrix<T>& A, const SparseMatrix<T>& B,
                   const DenseMatrix<T>& E, const DenseMatrix<T>& x);
// returns <Ax,x>/<(B-EE^T)x,x>

}  // namespace geometrycentral

#ifdef HAVE_SUITESPARSE
#include "suitesparse_wrapper/sparse_matrix_suitesparse.ipp"
#include "suitesparse_wrapper/sparse_solvers_suitesparse.ipp"
#else
#include "eigen_wrapper/sparse_matrix_eigen.ipp"
#include "eigen_wrapper/sparse_solvers_eigen.ipp"
#endif
