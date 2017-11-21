#pragma once

namespace geometrycentral {

const int DEFAULT_MAX_EIG_ITER = 32;
// maximum number of iterations used to solve eigenvalue problems

// Extra helper method declaration that SuiteSparse solvers make use of
template <class T>
void backsolvePositiveDefinite(SparseFactor<T>& L, DenseMatrix<T>& x,
                               DenseMatrix<T>& b);

// solves the sparse linear system Ax = b using sparse QR factorization
template <typename T>
void solve(SparseMatrix<T>& A, DenseMatrix<T>& x, DenseMatrix<T> b) {
  SuiteSparseQR_factorization<typename QRFactorization<T>::AUTO_ENTRYTYPE>*
      factorization = A.factors.qr.get();
  DenseMatrix<T> y;

  START_TIMING(solve)

  // solve y = R'\(E'*b)
  y = SuiteSparseQR_solve(SPQR_RTX_EQUALS_ETB, factorization, b.to_cholmod(),
                          context);

  // compute w = Q*y
  x = SuiteSparseQR_qmult(SPQR_QX, factorization, y.to_cholmod(), context);

  cout << "[qr] time: " << pretty_time(FINISH_TIMING(solve)) << "\n";
  cout << "[qr] max residual: " << residual(A, x, b) << "\n";
  cout << "[qr] size: " << A.nRows() << " x " << A.nColumns() << "\n";
  cout << "[qr] rank: " << (*context).SPQR_istat[4] << "\n";
}

// solves the sparse linear system Ax = b using sparse LU factorization
template <typename T>
void solveSquare(SparseMatrix<T>& A, DenseMatrix<T>& x, DenseMatrix<T> b) {
  cholmod_sparse* Ac = A.to_cholmod();
  SuiteSparse_long* Ap = (SuiteSparse_long*)Ac->p;
  SuiteSparse_long* Ai = (SuiteSparse_long*)Ac->i;
  double* Ax = (double*)Ac->x;
  void* factorization = A.factors.lu.get();

  START_TIMING(solve)
  umfpack_dl_solve(UMFPACK_A, Ap, Ai, Ax, (double*)&x(0), (double*)&b(0),
                   factorization, NULL, NULL);

  cout << "[lu] solve: " << pretty_time(FINISH_TIMING(solve)) << "\n";
  cout << "[lu] max residual: " << residual(A, x, b) << "\n";
}

// solves the positive definite sparse linear system Ax = b using sparse
// Cholesky factorization
template <typename T>
void solvePositiveDefinite(SparseMatrix<T>& A, DenseMatrix<T>& x,
                           DenseMatrix<T> b) {
  cholmod_factor* factorization = A.factors.ll.get();

  START_TIMING(solve)
  x = cholmod_l_solve(CHOLMOD_A, factorization, b.to_cholmod(), context);

  cout << "[ll] time: " << pretty_time(FINISH_TIMING(solve)) << "\n";
  cout << "[ll] max residual: " << residual(A, x, b) << "\n";
}

// solves the indefinite sparse linear system Ax = b using sparse LDL^T
// factorization
template <typename T>
void solveIndefinite(SparseMatrix<T>& A, DenseMatrix<T>& x, DenseMatrix<T> b) {
  cholmod_factor* factorization = A.factors.ldl.get();

  START_TIMING(solve)
  x = cholmod_l_solve(CHOLMOD_A, factorization, b.to_cholmod(), context);

  cout << "[ldl] time: " << pretty_time(FINISH_TIMING(solve)) << "\n";
  cout << "[ldl] max residual: " << residual(A, x, b) << "\n";
}

// solves the system LD'Lx = b, where LDL = A and D' is a modification
// of the diagonal such that the resulting system is positive-definite
template <typename T>
void solveModifiedIndefinite(SparseMatrix<T>& A, DenseMatrix<T>& x,
                             DenseMatrix<T> b) {
  cholmod_factor* factorization = A.factors.ldl.get();
  A.factors.ldl.makePositive();

  START_TIMING(solve)
  x = cholmod_l_solve(CHOLMOD_A, factorization, b.to_cholmod(), context);

  cout << "[ldl] time: " << pretty_time(FINISH_TIMING(solve)) << "\n";
  cout << "[ldl] max residual: " << residual(A, x, b) << "\n";
}

// solves the sparse linear system Ax = b, yielding the least-squares
// solution for over- or under-determined systems
template <typename T>
void solveLeastSquares(SparseMatrix<T>& A, DenseMatrix<T>& x,
                       DenseMatrix<T> b) {
  START_TIMING(solve)
  auto aa = A.to_cholmod();
  auto bb = b.to_cholmod();
  // x = SuiteSparseQR<typename QRFactorization<T>::AUTO_ENTRYTYPE>(
  // A.to_cholmod(), b.to_cholmod(), context );
  x = SuiteSparseQR<typename QRFactorization<T>::AUTO_ENTRYTYPE>(aa, bb,
                                                                 context);

  cout << "[qrls] time: " << pretty_time(FINISH_TIMING(solve)) << "\n";
  cout << "[qrls] max residual: " << residual(A, x, b) << "\n";
  cout << "[qrls] size: " << A.nRows() << " x " << A.nColumns() << "\n";
  cout << "[qrls] rank: " << (*context).SPQR_istat[4] << "\n";
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

// backsolves the prefactored positive definite sparse linear system LL'x = b
template <class T>
void backsolvePositiveDefinite(SparseFactor<T>& L, DenseMatrix<T>& x,
                               DenseMatrix<T>& b) {
  x = cholmod_l_solve(CHOLMOD_A, L.to_cholmod(), b.to_cholmod(), context);
}

// template<typename T>
// void solveLeastSquares(SparseMatrix<T>& A,
//                        DenseMatrix<T>& b,
//                        DenseMatrix<T>& x)
// {
//     SparseMatrix<T> At = A.transpose();
//     SparseMatrix<T> AtA = At * A;
//     DenseMatrix<T> Atb = At * b;
//     solvePositiveDefinite(AtA, Atb, x);
// }

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

// ## Sparse Factor methods
template <class T>
SparseFactor<T>::SparseFactor(void) : L(NULL) {}

template <class T>
SparseFactor<T>::~SparseFactor(void) {
  if (L) {
    cholmod_l_free_factor(&L, context);
  }
}

template <class T>
void SparseFactor<T>::printPermutation() {
  size_t n = L->n;
  size_t* P = (size_t*)L->Perm;
  size_t nnz = L->nzmax;

  cout << "PERMUTATION" << endl;
  cout << "===========" << endl;
  for (size_t i = 0; i < n; i++) {
    cout << i << ": " << P[i] << endl;
  }
  cout << "NONZEROS: " << nnz << endl;
  cout << endl;
}

template <class T>
void SparseFactor<T>::printOrderingMethod() {
  cout << "selected method: ";
  int method = ((cholmod_common*)context)->selected;
  switch (method) {
    case 0:
      cout << "[0] user-provided" << endl;
      break;
    case 1:
      cout << "[1] AMD" << endl;
      break;
    case 2:
      cout << "[2] METIS" << endl;
      break;
    case 3:
      cout << "[3] NESDIS (default)" << endl;
      break;
    case 4:
      cout << "[4] natural" << endl;
      break;
    case 5:
      cout << "[5] NESDIS, nd_small=20000" << endl;
      break;
    case 6:
      cout << "[6] NESDIS, nd_small=4, no constrained minimum degree" << endl;
      break;
    case 7:
      cout << "[7] NESDIS, no dense node removal" << endl;
      break;
    case 8:
      cout << "[8] AMD for A, COLAMD for A*A'" << endl;
      break;
    default:
      cout << "(unknown)" << endl;
      break;
  }
}

template <class T>
void SparseFactor<T>::build(SparseMatrix<T>& A, size_t nI) {
  if (L) {
    cholmod_l_free_factor(&L, context);
    L = NULL;
  }

  cholmod_sparse* Ac = A.to_cholmod();
  Ac->stype = 1;

  // TODO modify permutation
  if (0)  // nI > 0 )
  {
    ((cholmod_common*)context)->nmethods = 1;
    ((cholmod_common*)context)->method[0].ordering = 0;  // user-provided
    SuiteSparse_long n = Ac->nrow;
    vector<SuiteSparse_long> P(n);
    for (SuiteSparse_long i = 0; i < n; i++) {
      P[i] = i;  // XXX need better ordering!
    }
    L = cholmod_l_analyze_p(Ac, &P[0], NULL, 0, context);
  } else {
    // Seems that for cotan-Laplace, AMD is typically (always?) the best choice
    ((cholmod_common*)context)->nmethods = 1;
    ((cholmod_common*)context)->method[0].ordering = CHOLMOD_AMD;
    L = cholmod_l_analyze(Ac, context);
  }

  cholmod_l_factorize(Ac, L, context);
}

template <class T>
void SparseFactor<T>::print() const {
  cout << "=================================================" << endl;
  cout << "SPARSE FACTOR" << endl;
  cout << "=================================================" << endl;
  cout << "for both simplicial and supernodal factorizations" << endl;
  cout << "          n: " << L->n << endl;
  cout << "      minor: " << L->minor << endl;
  cout << "symbolic ordering and analysis" << endl;
  cout << "          Perm: " << L->Perm << endl;
  cout << "      ColCount: " << L->ColCount << endl;
  cout << "         IPerm: " << L->IPerm << endl;
  cout << "simplicial factorization" << endl;
  cout << "      nzmax: " << L->nzmax << endl;
  cout << "         *p: " << L->p << endl;
  cout << "         *i: " << L->i << endl;
  cout << "         *x: " << L->x << endl;
  cout << "         *z: " << L->z << endl;
  cout << "        *nz: " << L->nz << endl;
  cout << "      *next: " << L->next << endl;
  cout << "      *prev: " << L->prev << endl;
  cout << "supernodal factorization" << endl;
  cout << "        nsuper: " << L->nsuper << endl;
  cout << "         ssize: " << L->ssize << endl;
  cout << "         xsize: " << L->xsize << endl;
  cout << "      maxcsize: " << L->maxcsize << endl;
  cout << "      maxesize: " << L->maxesize << endl;
  cout << "        *super: " << L->super << endl;
  cout << "           *pi: " << L->pi << endl;
  cout << "           *px: " << L->px << endl;
  cout << "            *s: " << L->s << endl;
  cout << "factorization type" << endl;
  cout << "          ordering: " << L->ordering << endl;
  cout << "             is_ll: " << L->is_ll << endl;
  cout << "          is_super: " << L->is_super << endl;
  cout << "      is_monotonic: " << L->is_monotonic << endl;
  cout << "             itype: " << L->itype << endl;
  cout << "             xtype: " << L->xtype << endl;
  cout << "             dtype: " << L->dtype << endl;
  cout << endl;

  if (L->nzmax != 0) {
    cout << "--------------------------------------" << endl;
    cout << " Simplicial" << endl;
    cout << "--------------------------------------" << endl;
    cout << "p: " << endl;
    for (int i = 0; i < L->n + 1; i++) {
      cout << i << ": " << ((SuiteSparse_long*)(L->p))[i] << endl;
    }
    cout << "--------------------------------------" << endl;
    cout << "i: " << endl;
    for (int i = 0; i < L->nzmax; i++) {
      cout << i << ": " << ((SuiteSparse_long*)(L->i))[i] << endl;
    }
    cout << "--------------------------------------" << endl;
    cout << "x: " << endl;
    for (int i = 0; i < L->nzmax; i++) {
      cout << i << ": " << ((double*)(L->x))[i] << endl;
    }
    cout << "--------------------------------------" << endl;
    cout << "nz: " << endl;
    for (int i = 0; i < L->n; i++) {
      cout << i << ": " << ((SuiteSparse_long*)(L->nz))[i] << endl;
    }
    cout << "--------------------------------------" << endl;
    cout << "next: " << endl;
    for (int i = 0; i < L->n + 2; i++) {
      cout << i << ": " << ((SuiteSparse_long*)(L->next))[i] << endl;
    }
    cout << "--------------------------------------" << endl;
    cout << "prev: " << endl;
    for (int i = 0; i < L->n + 2; i++) {
      cout << i << ": " << ((SuiteSparse_long*)(L->prev))[i] << endl;
    }
    cout << endl;
  } else {
    cout << "--------------------------------------" << endl;
    cout << " Supernodal" << endl;
    cout << "--------------------------------------" << endl;
    cout << "super: " << endl;
    for (int i = 0; i < L->nsuper + 1; i++) {
      cout << i << ": " << ((SuiteSparse_long*)(L->super))[i] << endl;
    }
    cout << "--------------------------------------" << endl;
    cout << "pi: " << endl;
    for (int i = 0; i < L->nsuper + 1; i++) {
      cout << i << ": " << ((SuiteSparse_long*)(L->pi))[i] << endl;
    }
    cout << "--------------------------------------" << endl;
    cout << "px: " << endl;
    for (int i = 0; i < L->nsuper + 1; i++) {
      cout << i << ": " << ((SuiteSparse_long*)(L->px))[i] << endl;
    }
    cout << "--------------------------------------" << endl;
    cout << "s: " << endl;
    for (int i = 0; i < L->ssize; i++) {
      cout << i << ": " << ((SuiteSparse_long*)(L->s))[i] << endl;
    }
    cout << "--------------------------------------" << endl;
    cout << "x: " << endl;
    for (int i = 0; i < L->xsize; i++) {
      cout << i << ": " << ((double*)(L->x))[i] << endl;
    }
    cout << endl;

    cout << "L = {";
    SuiteSparse_long* Lsuper = (SuiteSparse_long*)L->super;
    SuiteSparse_long* Lpi = (SuiteSparse_long*)L->pi;
    SuiteSparse_long* Ls = (SuiteSparse_long*)L->s;
    double* Lx = (double*)L->x;
    for (size_t node = 0; node < L->nsuper; node++)  // iterate over supernodes
    {
      for (int col = Lsuper[node]; col < Lsuper[node + 1];
           col++)  // iterate over columns in this supernode
      {
        for (int k = Lpi[node]; k < Lpi[node + 1];
             k++)  // iterate over nonzero row indices for this supernode
        {
          size_t row = Ls[k];  // get the current row index

          // emit the current entry
          cout << "{" << (1 + row) << "," << (1 + col) << "}->" << (*Lx) << ",";
          Lx++;
        }
      }
    }
    cout << "\b};" << endl;
  }

  cout << "=================================================" << endl;
  cout << endl;
}

template <class T>
bool SparseFactor<T>::valid(void) const {
  if (L == NULL) {
    return false;
  }
  return true;
}

template <class T>
cholmod_factor* SparseFactor<T>::to_cholmod(void) {
  return L;
}
}