#pragma once

namespace geometrycentral {

const int DEFAULT_MAX_EIG_ITER = 32;
// maximum number of iterations used to solve eigenvalue problems

// solves the sparse linear system Ax = b using sparse QR factorization
template <typename T>
void solve(SparseMatrix<T>& A, DenseMatrix<T>& x, DenseMatrix<T> b) {
  START_TIMING(solve)

  A.factors.qr.solve(x, b);

  cout << "[qr] time: " << pretty_time(FINISH_TIMING(solve)) << "\n";
  cout << "[qr] max residual: " << residual(A, x, b) << "\n";
  cout << "[qr] size: " << A.nRows() << " x " << A.nColumns() << "\n";
  cout << "[qr] rank: " << A.factors.qr.rank() << "\n";
}

// solves the sparse linear system Ax = b using sparse LU factorization
template <typename T>
void solveSquare(SparseMatrix<T>& A, DenseMatrix<T>& x, DenseMatrix<T> b) {
  START_TIMING(solve)
  A.factors.lu.solve(x, b);

  cout << "[lu] solve: " << pretty_time(FINISH_TIMING(solve)) << "\n";
  cout << "[lu] max residual: " << residual(A, x, b) << "\n";
}

// solves the positive definite sparse linear system Ax = b using sparse
// Cholesky factorization
template <typename T>
void solvePositiveDefinite(SparseMatrix<T>& A, DenseMatrix<T>& x,
                           DenseMatrix<T> b) {
  START_TIMING(solve)
  A.factors.ll.solve(x, b);

  cout << "[ll] time: " << pretty_time(FINISH_TIMING(solve)) << "\n";
  cout << "[ll] max residual: " << residual(A, x, b) << "\n";
}

// solves the indefinite sparse linear system Ax = b using sparse LDL^T
// factorization
template <typename T>
void solveIndefinite(SparseMatrix<T>& A, DenseMatrix<T>& x, DenseMatrix<T> b) {
  START_TIMING(solve)
  A.factors.ldl.solve(x, b);

  cout << "[ldl] time: " << pretty_time(FINISH_TIMING(solve)) << "\n";
  cout << "[ldl] max residual: " << residual(A, x, b) << "\n";
}

// solves the system LD'Lx = b, where LDL = A and D' is a modification
// of the diagonal such that the resulting system is positive-definite
template <typename T>
void solveModifiedIndefinite(SparseMatrix<T>& A, DenseMatrix<T>& x,
                             DenseMatrix<T> b) {
  throw geometrycentral::FunctionalityException("Not implmented with this matrix backend");
}

// solves the sparse linear system Ax = b, yielding the least-squares
// solution for over- or under-determined systems
template <typename T>
void solveLeastSquares(SparseMatrix<T>& A, DenseMatrix<T>& x,
                       DenseMatrix<T> b) {
  START_TIMING(solve)
  A.factors.qr.solve(x, b);

  cout << "[qrls] time: " << pretty_time(FINISH_TIMING(solve)) << "\n";
  cout << "[qrls] max residual: " << residual(A, x, b) << "\n";
  cout << "[qrls] size: " << A.nRows() << " x " << A.nColumns() << "\n";
  cout << "[qrls] rank: " << A.factors.qr.rank() << "\n";
}

template <typename T>
void solveConjugateGradient(SparseMatrix<T>& A, SparseMatrix<T>& M,
                            DenseMatrix<T>& x, DenseMatrix<T>& b,
                            double relativeTolerance, int maxIterations) {
  double relativeResidual = -1.;
  DenseMatrix<T> r = b - A * x;

  DenseMatrix<T> z;
  solvePositiveDefinite(M, z, r);

  DenseMatrix<T> p = z;

  int iter;
  for (iter = 1; iter < maxIterations; iter++) {
    T zr0 = inner(z, r);
    T alpha = inner(p, A * p).inv() * zr0;
    DenseMatrix<T> alphap = alpha * p;
    x += alphap;
    r -= A * alphap;

    relativeResidual = r.norm(lInfinity) / b.norm(lInfinity);
    if (relativeResidual < relativeTolerance) {
      break;
    }

    solvePositiveDefinite(M, z, r);

    T beta = zr0.inv() * inner(z, r);
    DenseMatrix<T> betap = beta * p;
    p = z + betap;
  }

  cout << "[pcg] achieved a residual of " << relativeResidual << " after "
       << iter << " iterations" << endl;
}

// solves A x = lambda x for the smallest nonzero eigenvalue lambda
// A must be square; x is used as an initial guess
template <typename T>
void smallestEig(SparseMatrix<T>& A, DenseMatrix<T>& x,
                 bool ignoreConstantVector) {
  START_TIMING(solve)

  for (int iter = 0; iter < DEFAULT_MAX_EIG_ITER; iter++) {
    solveSquare(A, x, x);
    if (ignoreConstantVector) {
      x.removeMean();
    }
    x.normalize();
  }

  cout << "[eig] time: " << pretty_time(FINISH_TIMING(solve)) << "\n";
  cout << "[eig] max residual: " << residual(A, x) << "\n";
}

// solves A x = lambda x for the smallest nonzero eigenvalue lambda
// A must be positive (semi-)definite; x is used as an initial guess
template <typename T>
void smallestEigPositiveDefinite(SparseMatrix<T>& A, DenseMatrix<T>& x,
                                 bool ignoreConstantVector) {
  START_TIMING(solve)
  for (int iter = 0; iter < DEFAULT_MAX_EIG_ITER; iter++) {
    solvePositiveDefinite(A, x, x);
    if (ignoreConstantVector) {
      x.removeMean();
    }
    x.normalize();
  }

  cout << "[eig] time: " << pretty_time(FINISH_TIMING(solve)) << "\n";
  cout << "[eig] max residual: " << residual(A, x) << "\n";
}

template <typename T>
void smallestEigPCG(SparseMatrix<T>& A,
                    SparseMatrix<T>& M,  // preconditioner
                    DenseMatrix<T>& x, bool ignoreConstantVector) {
  START_TIMING(solve)

  for (int iter = 0; iter < DEFAULT_MAX_EIG_ITER; iter++) {
    DenseMatrix<T> b = x;
    solveConjugateGradient(A, M, x, b, 1e-2, 10000);
    if (ignoreConstantVector) {
      x.removeMean();
    }
    x.normalize();
  }

  cout << "[eig] time: " << pretty_time(FINISH_TIMING(solve)) << "\n";
  cout << "[eig] max residual: " << residual(A, x) << "\n";
}

// solves A x = lambda (B - EE^T) x for the smallest nonzero eigenvalue lambda
// A must be positive (semi-)definite, B must be symmetric; EE^T is a low-rank
// matrix, and
// x is used as an initial guess
template <typename T>
void smallestEigPositiveDefinite(SparseMatrix<T>& A, SparseMatrix<T>& B,
                                 DenseMatrix<T>& e, DenseMatrix<T>& x) {
  int iter;
  DenseMatrix<T> eT = e.transpose();

  for (iter = 0; iter < DEFAULT_MAX_EIG_ITER; iter++) {
    x = B * x - e * (eT * x);
    solvePositiveDefinite(A, x, x);
    x.normalize();
  }

  cout << "[eig] max residual: " << residual(A, B, e, x) << "\n";
}
}