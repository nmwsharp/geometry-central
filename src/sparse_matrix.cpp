#include "sparse_matrix.h"

namespace GC {
#ifdef HAVE_SUITESPARSE
// This routine has been used during debugging to make sure I'm
// not doing anything stupid during memory allocation.
void safealloc(void*& p, size_t size) {
  if (p)  // If this pointer was allocated already, deallocate and set to NULL.
  {
    free(p);
    p = NULL;
  }
  // XXX TOTAL KLUDGE!!! we allocate 2x as much memory
  // as we actually need to avoid some memory-related error...
  p = malloc(size * 2);
  // A good way to debug would be to see *which* array
  // needs this "double allocation," then track down a more
  // specific addressing error.
  // XXX TOTAL KLUDGE!!!
}

template <>
void SparseFactor<double>::subfactor_simplicial(size_t m,
                                                SparseFactor<double>& factor) {
  cholmod_factor*& M = factor.L;

  // Allocate the factor
  M = cholmod_l_allocate_factor(m, context);

  if (m > L->n) {
    cerr << "Error: cannot extract a subfactor larger than the original factor!"
         << endl;
    exit(1);
  }

  if (L->minor != L->n) {
    cerr << "Warning: attempting to extract a subfactor of a failed "
            "factorization."
         << endl;
  }

  M->n = m;
  M->minor = m;

  // Copy permutation used for original factorization---here we make the
  // assumption that the original permutation preserves the upper-left m x
  // m block of the matrix.
  safealloc(M->Perm, m * sizeof(SuiteSparse_long));
  SuiteSparse_long* LPerm = (SuiteSparse_long*)(L->Perm);
  SuiteSparse_long* MPerm = (SuiteSparse_long*)(M->Perm);
  for (size_t c = 0; c < m; c++) {
    MPerm[c] = LPerm[c];
  }

  size_t nzmax = 0;

  vector<SuiteSparse_long> i;
  vector<double> x;
  vector<SuiteSparse_long> p(m + 1), nz(m);
  SuiteSparse_long* Lp = (SuiteSparse_long*)(L->p);
  SuiteSparse_long* Lnz = (SuiteSparse_long*)(L->nz);
  SuiteSparse_long* Li = (SuiteSparse_long*)(L->i);
  double* Lx = (double*)(L->x);
  for (size_t c = 0; c < m; c++)  // iterate over first m columns
  {
    p[c] = nzmax;  // start index for this column
    nz[c] = 0;     // count the number of nonzeros in this column

    // iterate over nonzeros in current column of L
    for (int e = Lp[c]; e < Lp[c] + Lnz[c]; e++) {
      // if the row index is in the sub-block
      if (Li[e] < (int)m) {
        i.push_back(Li[e]);
        x.push_back(Lx[e]);
        nz[c]++;
        nzmax++;
      }
    }
  }
  p[m] = nzmax;  // end index for the final column

  M->nzmax = nzmax;

  // Copy temporary arrays into factor data structure
  safealloc(M->p, (m + 1) * sizeof(SuiteSparse_long));
  SuiteSparse_long* Mp = (SuiteSparse_long*)(M->p);
  for (size_t c = 0; c <= m; c++) {
    Mp[c] = p[c];
  }

  safealloc(M->i, nzmax * sizeof(SuiteSparse_long));
  safealloc(M->x, nzmax * sizeof(double));
  SuiteSparse_long* Mi = (SuiteSparse_long*)(M->i);
  double* Mx = (double*)(M->x);
  for (size_t e = 0; e < nzmax; e++) {
    Mi[e] = i[e];
    Mx[e] = x[e];
  }

  safealloc(M->nz, m * sizeof(SuiteSparse_long));
  SuiteSparse_long* Mnz = (SuiteSparse_long*)(M->nz);
  for (size_t c = 0; c < m; c++) {
    Mnz[c] = nz[c];
  }

  safealloc(M->next, (m + 2) * sizeof(SuiteSparse_long));
  safealloc(M->prev, (m + 2) * sizeof(SuiteSparse_long));

  SuiteSparse_long* next = (SuiteSparse_long*)(M->next);
  for (size_t c = 0; c < m; c++) {
    next[c] = c + 1;
  }
  next[m] = -1;
  next[m + 1] = 0;

  SuiteSparse_long* prev = (SuiteSparse_long*)(M->prev);
  prev[0] = m + 1;
  for (size_t c = 1; c < m + 1; c++) {
    prev[c] = c - 1;
  }
  prev[m + 1] = -1;

  // Set the column count to the number of nonzeros in each column
  safealloc(M->ColCount, m * sizeof(SuiteSparse_long));
  SuiteSparse_long* ColCount = (SuiteSparse_long*)M->ColCount;
  for (size_t c = 0; c < m; c++) {
    ColCount[c] = nz[c];
  }

  // Zero out parameters for supernodal factorization
  M->nsuper = 0;
  M->ssize = 0;
  M->xsize = 0;
  M->maxcsize = 0;
  M->maxesize = 0;
  M->super = NULL;
  M->pi = NULL;
  M->px = NULL;
  M->s = NULL;

  // Specify factorization type
  M->ordering = 1;
  M->is_ll = 0;
  M->is_super = 0;
  M->is_monotonic = 1;
  M->itype = 2;
  M->xtype = 1;
  M->dtype = 0;

  // factor.print();
}

template <>
void SparseFactor<double>::subfactor_supernodal(size_t m,
                                                SparseFactor<double>& factor) {
  cholmod_factor*& M = factor.L;

  // Allocate the factor
  M = cholmod_l_allocate_factor(m, context);

  if (m > L->n) {
    cerr << "Error: cannot extract a subfactor larger than the original factor!"
         << endl;
    exit(1);
  }

  if (L->minor != L->n) {
    cerr << "Warning: attempting to extract a subfactor of a failed "
            "factorization."
         << endl;
  }

  M->n = m;
  M->minor = m;

  // Copy permutation used for original factorization---here we make the
  // assumption that the original permutation preserves the upper-left m x
  // m block of the matrix.
  safealloc(M->Perm, m * sizeof(SuiteSparse_long));
  SuiteSparse_long* LPerm = (SuiteSparse_long*)(L->Perm);
  SuiteSparse_long* MPerm = (SuiteSparse_long*)(M->Perm);
  for (size_t c = 0; c < m; c++) {
    MPerm[c] = LPerm[c];
  }

  // Zero out parameters for simplicial factorization
  M->nzmax = 0;
  M->p = NULL;
  M->i = NULL;
  M->x = NULL;
  M->z = NULL;
  M->nz = NULL;
  M->next = NULL;
  M->prev = NULL;

  // Zero out parameters for supernodal factorization
  M->nsuper = 0;
  M->ssize = 0;
  M->xsize = 0;
  M->maxcsize = 0;
  M->maxesize = 0;
  M->super = NULL;
  M->pi = NULL;
  M->px = NULL;
  M->s = NULL;
  M->x = NULL;

  SuiteSparse_long* Lsuper = (SuiteSparse_long*)L->super;
  SuiteSparse_long* Lpi = (SuiteSparse_long*)L->pi;
  SuiteSparse_long* Lpx = (SuiteSparse_long*)L->px;
  SuiteSparse_long* Ls = (SuiteSparse_long*)L->s;
  double* Lx = (double*)L->x;

  // Get the list of supernodes whose first column is contained in the extracted
  // block
  vector<SuiteSparse_long> super;
  for (size_t i = 0; i < L->nsuper; i++) {
    if (Lsuper[i] < (int)m) {
      super.push_back(Lsuper[i]);
    }
  }
  M->nsuper = super.size();  // number of supernodes in the subfactor
  super.push_back(m);        // pad the end with "one past the last column"
  safealloc(M->super, (M->nsuper + 1) * sizeof(SuiteSparse_long));
  SuiteSparse_long* Msuper = (SuiteSparse_long*)M->super;

  // Now that we know how many supernodes are in the subfactor,
  // copy our dynamic array super into the static array M->super.
  for (size_t i = 0; i < M->nsuper + 1; i++) {
    Msuper[i] = super[i];
  }

  vector<SuiteSparse_long> s;
  safealloc(M->pi, (M->nsuper + 1) * sizeof(SuiteSparse_long));
  SuiteSparse_long* Mpi = (SuiteSparse_long*)M->pi;
  Mpi[0] = 0;
  for (size_t node = 0; node < M->nsuper;
       node++)  // iterate over supernodes in the subfactor
  {
    for (int k = Lpi[node]; k < Lpi[node + 1];
         k++)  // iterate over row indices in the original supernode
    {
      int row = Ls[k];
      if (row < (int)m) {
        s.push_back(row);
      }
    }
    Mpi[node + 1] = s.size();
  }
  // Set the number of row indices to the size of "s".
  M->ssize = s.size();

  // Now that we know the number of row indices, copy our dynamic
  // array s into the static array M->s.  Note that the index
  // pointers for each supernode were already set above (Mpi).
  safealloc(M->s, M->ssize * sizeof(SuiteSparse_long));
  SuiteSparse_long* Ms = (SuiteSparse_long*)M->s;
  for (size_t i = 0; i < M->ssize; i++) {
    Ms[i] = s[i];
  }

  // Set x, px, and xsize; allocate the first two
  safealloc(M->px, (M->nsuper + 1) * sizeof(SuiteSparse_long));
  SuiteSparse_long* Mpx = (SuiteSparse_long*)M->px;
  vector<double> x;
  for (size_t node = 0; node < M->nsuper; node++) {
    size_t nLRows = Lpi[node + 1] - Lpi[node];
    size_t nMRows = Mpi[node + 1] - Mpi[node];
    size_t nSkip = nLRows - nMRows;
    double* Lxi = Lx + Lpx[node];
    for (int col = Msuper[node]; col < Msuper[node + 1]; col++) {
      for (size_t k = 0; k < nMRows; k++) {
        x.push_back(*Lxi);
        Lxi++;
      }
      Lxi += nSkip;
    }
    Mpx[node + 1] = x.size();
  }
  M->xsize = x.size();
  safealloc(M->x, (M->xsize) * sizeof(double));
  double* Mx = (double*)M->x;
  for (size_t i = 0; i < M->xsize; i++) {
    Mx[i] = x[i];
  }

  // Construct a map from columns to their (M-)supernodes
  vector<SuiteSparse_long> SuperMap(m);
  size_t currentNode = 0;
  for (size_t column = 0; column < m; column++) {
    if ((int)column == Msuper[currentNode + 1]) {
      currentNode++;
    }
    SuperMap[column] = currentNode;
  }

  // Compute value of maxcsize and maxesize
  M->maxcsize = 1;
  M->maxesize = 1;
  for (size_t d = 0; d < M->nsuper; d++)  // iterate over supernodes
  {
    size_t nscol = Msuper[d + 1] - Msuper[d];
    size_t P = Mpi[d] + nscol;
    size_t plast = P;
    size_t pend = Mpi[d + 1];
    size_t esize = pend - P;
    M->maxesize = max(M->maxesize, esize);
    size_t slast = (P == pend) ? (-1) : (SuperMap[Ls[P]]);
    for (; P <= pend; P++) {
      size_t S = (P == pend) ? (-1) : (SuperMap[Ls[P]]);
      if (S != slast) {
        /* row i is the start of a new supernode */
        size_t ndrow1 = P - plast;
        size_t ndrow2 = pend - plast;
        size_t csize = ndrow2 * ndrow1;
        M->maxcsize = max(M->maxcsize, csize);
        plast = P;
        slast = S;
      }
    }
  }

  // Specify factorization type
  M->ordering = 1;
  M->is_ll = 1;
  M->is_super = 1;
  M->is_monotonic = 1;
  M->itype = 2;
  M->xtype = 1;
  M->dtype = 0;

  // factor.print();
}

template <>
void SparseFactor<double>::subfactor(size_t m, SparseFactor<double>& factor) {
  if (L->nzmax != 0)  // is it simplicial?
  {
    subfactor_simplicial(m, factor);
  } else  // otherwise, it's supernodal
  {
    subfactor_supernodal(m, factor);
  }
}
#endif

template <>
SparseMatrix<double> SparseMatrix<Quaternion>::toReal()
// converts this matrix to its real representation
{
  SparseMatrix<double> A(m * 4, n * 4);

  for (SparseMatrix<Quaternion>::const_iterator e = begin(); e != end(); e++) {
    int i = e->first.row;
    int j = e->first.col;
    const Quaternion& q(e->second);

    A(i * 4 + 0, j * 4 + 0) = q[0];
    A(i * 4 + 0, j * 4 + 1) = -q[1];
    A(i * 4 + 0, j * 4 + 2) = -q[2];
    A(i * 4 + 0, j * 4 + 3) = -q[3];
    A(i * 4 + 1, j * 4 + 0) = q[1];
    A(i * 4 + 1, j * 4 + 1) = q[0];
    A(i * 4 + 1, j * 4 + 2) = -q[3];
    A(i * 4 + 1, j * 4 + 3) = q[2];
    A(i * 4 + 2, j * 4 + 0) = q[2];
    A(i * 4 + 2, j * 4 + 1) = q[3];
    A(i * 4 + 2, j * 4 + 2) = q[0];
    A(i * 4 + 2, j * 4 + 3) = -q[1];
    A(i * 4 + 3, j * 4 + 0) = q[3];
    A(i * 4 + 3, j * 4 + 1) = -q[2];
    A(i * 4 + 3, j * 4 + 2) = q[1];
    A(i * 4 + 3, j * 4 + 3) = q[0];
  }

  return A;
}

template <>
SparseMatrix<Quaternion> SparseMatrix<double>::toQuaternionic()
// converts this matrix to its quaternionic representation
{
  SparseMatrix<Quaternion> A(m, n);

  for (SparseMatrix<double>::const_iterator e = begin(); e != end(); e++) {
    int i = e->first.row;
    int j = e->first.col;
    double val = e->second;

    A(i, j) = val;
  }

  return A;
}

#ifdef HAVE_SUITESPARSE
template <>
const SparseMatrix<double>& SparseMatrix<double>::operator=(cholmod_sparse* B) {
  assert(B);
  assert(B->xtype == CHOLMOD_REAL);

  if (cData) {
    cholmod_l_free_sparse(&cData, context);
  }
  cData = B;

  m = cData->nrow;
  n = cData->ncol;
  resize(m, n);

  double* pr = (double*)cData->x;
  SuiteSparse_long* ir = (SuiteSparse_long*)cData->i;
  SuiteSparse_long* jc = (SuiteSparse_long*)cData->p;

  // iterate over columns
  for (size_t col = 0; col < n; col++) {
    // iterate over nonzero rows
    for (int k = jc[col]; k < jc[col + 1]; k++) {
      int row = ir[k];

      (*this)(row, col) = pr[k];
    }
  }

  return *this;
}

template <>
const SparseMatrix<Complex>& SparseMatrix<Complex>::operator=(
    cholmod_sparse* B) {
  assert(B);
  assert(B->xtype == CHOLMOD_COMPLEX);

  if (cData) {
    cholmod_l_free_sparse(&cData, context);
  }
  cData = B;

  m = cData->nrow;
  n = cData->ncol;
  resize(m, n);

  double* pr = (double*)cData->x;
  SuiteSparse_long* ir = (SuiteSparse_long*)cData->i;
  SuiteSparse_long* jc = (SuiteSparse_long*)cData->p;

  // iterate over columns
  for (size_t col = 0; col < n; col++) {
    // iterate over nonzero rows
    for (int k = jc[col]; k < jc[col + 1]; k++) {
      int row = ir[k];

      (*this)(row, col) = Complex(pr[k * 2 + 0], pr[k * 2 + 1]);
    }
  }

  return *this;
}

template <>
void SparseMatrix<double>::allocateSparse(void) {
  int nzmax = data.size();
  int sorted = true;
  int packed = true;
  int stype = 0;
  cData = cholmod_l_allocate_sparse(m, n, nzmax, sorted, packed, stype,
                                    CHOLMOD_REAL, context);
}

template <>
void SparseMatrix<Complex>::allocateSparse(void) {
  int nzmax = data.size();
  int sorted = true;
  int packed = true;
  int stype = 0;
  cData = cholmod_l_allocate_sparse(m, n, nzmax, sorted, packed, stype,
                                    CHOLMOD_COMPLEX, context);
}

template <>
void SparseMatrix<double>::setEntry(const_iterator e, int i, double* pr) {
  pr[i] = e->second;
}

template <>
void SparseMatrix<Complex>::setEntry(const_iterator e, int i, double* pr) {
  pr[i * 2 + 0] = e->second.real();
  pr[i * 2 + 1] = e->second.imag();
}
#endif

// TODO re-implement via Eigen
//    template <>
//    void solveSquare( SparseMatrix<Complex>& A,
//                       DenseMatrix<Complex>& x,
//                       DenseMatrix<Complex>  b )
//    // solves the sparse linear system Ax = b using sparse LU factorization
//    {
//       cholmod_sparse* Ac = A.to_cholmod();
//       UF_long* Ap = (UF_long*) Ac->p;
//       UF_long* Ai = (UF_long*) Ac->i;
//       double*  Ax =  (double*) Ac->x;
//       void* factorization  = A.factors.lu.get();
//
// #ifndef NDEBUG
//       int t0 = clock();
// #endif
//       umfpack_zl_solve( UMFPACK_A, Ap, Ai, Ax, NULL, (double*) &x(0), NULL,
//       (double*) &b(0), NULL, factorization, NULL, NULL );
// #ifndef NDEBUG
//       int t1 = clock();
//
//       cout << "[lu] solve: " << seconds( t0, t1 ) << "s" << "\n";
//       cout << "[lu] max residual: " << residual( A, x, b ) << "\n";
// #endif
//    }
//
//    template <>
//    void solve( SparseMatrix<double>& A,
//                 DenseMatrix<double>& x,
//                 DenseMatrix<double>  b )
//    // solves the sparse linear system Ax = b using sparse QR factorization
//    {
//       int t0 = clock();
//       x = SuiteSparseQR<double>( A.to_cholmod(), b.to_cholmod(), context );
//       int t1 = clock();
//
//       cout << "[qr] time: " << seconds( t0, t1 ) << "s" << "\n";
//       cout << "[qr] max residual: " << residual( A, x, b ) << "\n";
//       cout << "[qr] size: " << A.nRows() << " x " << A.nColumns() << "\n";
//       cout << "[qr] rank: " << (*context).SPQR_istat[4] << "\n";
//    }
}
