#ifdef GC_HAVE_SUITESPARSE
#include "geometrycentral/numerical/suitesparse_utilities.h"

#include "geometrycentral/utilities/utilities.h"

#include <complex>
#include <iostream>

using std::cout;
using std::endl;

namespace geometrycentral {

// === Context
CholmodContext::CholmodContext(void) { cholmod_l_start(&context); }

CholmodContext::~CholmodContext(void) { cholmod_l_finish(&context); }

void CholmodContext::setSimplicial(void) { context.supernodal = CHOLMOD_SIMPLICIAL; }

void CholmodContext::setSupernodal(void) { context.supernodal = CHOLMOD_SUPERNODAL; }

void CholmodContext::setLL(void) { context.final_ll = true; }

void CholmodContext::setLDL(void) { context.final_ll = false; }


CholmodContext::operator cholmod_common*(void) { return &context; }

// === Conversion functions

// Helper to manage stypes
// note that if the matrix if symmetric, we still actually store the whole thing, but we need to make the stype
// symmetric flag because some Cholmod routines have different behavior depending on the flag. In the future, we could
// omit lower triangular in this case.
namespace {
int flagForStype(SType s) {
  switch (s) {
  case SType::UNSYMMETRIC:
    return 0; // unsymmetric, use whole
    break;
  case SType::SYMMETRIC:
    return 1; // symmetric, use upper triangular only
    break;
  }
  return -1;
}
} // namespace

// double-valued sparse matrices
template <>
cholmod_sparse* toCholmod(SparseMatrix<double>& A, CholmodContext& context, SType stype) {

  A.makeCompressed();

  // Allocate spase
  size_t Nentries = A.nonZeros();
  size_t Ncols = A.cols();
  size_t Nrows = A.rows();

  cholmod_sparse* cMat =
      cholmod_l_allocate_sparse(Nrows, Ncols, Nentries, true, true, flagForStype(stype), CHOLMOD_REAL, context);

  // Pull out useful pointers
  double* values = (double*)cMat->x;
  SuiteSparse_long* rowIndices = (SuiteSparse_long*)cMat->i;
  SuiteSparse_long* colStart = (SuiteSparse_long*)cMat->p;

  // Copy
  for (size_t iEntry = 0; iEntry < Nentries; iEntry++) {
    values[iEntry] = A.valuePtr()[iEntry];
    rowIndices[iEntry] = A.innerIndexPtr()[iEntry];
  }
  for (size_t iCol = 0; iCol < Ncols; iCol++) {
    colStart[iCol] = A.outerIndexPtr()[iCol];
  }
  colStart[Ncols] = Nentries;

  return cMat;
}

// float-valued sparse matrices
// CHOLMOD only uses CHOLMOD_REAL precision, so you always get one of those back regardless of float/double input)
template <>
cholmod_sparse* toCholmod(SparseMatrix<float>& A, CholmodContext& context, SType stype) {

  A.makeCompressed();

  // Allocate spase
  size_t Nentries = A.nonZeros();
  size_t Ncols = A.cols();
  size_t Nrows = A.rows();

  cholmod_sparse* cMat =
      cholmod_l_allocate_sparse(Nrows, Ncols, Nentries, true, true, flagForStype(stype), CHOLMOD_REAL, context);

  // Pull out useful pointers
  double* values = (double*)cMat->x;
  SuiteSparse_long* rowIndices = (SuiteSparse_long*)cMat->i;
  SuiteSparse_long* colStart = (SuiteSparse_long*)cMat->p;

  // Copy
  for (size_t iEntry = 0; iEntry < Nentries; iEntry++) {
    values[iEntry] = static_cast<double>(A.valuePtr()[iEntry]);
    rowIndices[iEntry] = A.innerIndexPtr()[iEntry];
  }
  for (size_t iCol = 0; iCol < Ncols; iCol++) {
    colStart[iCol] = A.outerIndexPtr()[iCol];
  }
  colStart[Ncols] = Nentries;

  return cMat;
}

// Complex-valued sparse matrices
template <>
cholmod_sparse* toCholmod(SparseMatrix<std::complex<double>>& A, CholmodContext& context, SType stype) {

  A.makeCompressed();

  // Allocate spase
  size_t Nentries = A.nonZeros();
  size_t Ncols = A.cols();
  size_t Nrows = A.rows();

  cholmod_sparse* cMat =
      cholmod_l_allocate_sparse(Nrows, Ncols, Nentries, true, true, flagForStype(stype), CHOLMOD_COMPLEX, context);

  // Pull out useful pointers
  std::complex<double>* values = (std::complex<double>*)cMat->x;
  SuiteSparse_long* rowIndices = (SuiteSparse_long*)cMat->i;
  SuiteSparse_long* colStart = (SuiteSparse_long*)cMat->p;

  // Copy
  for (size_t iEntry = 0; iEntry < Nentries; iEntry++) {
    values[iEntry] = A.valuePtr()[iEntry];
    rowIndices[iEntry] = A.innerIndexPtr()[iEntry];
  }
  for (size_t iCol = 0; iCol < Ncols; iCol++) {
    colStart[iCol] = A.outerIndexPtr()[iCol];
  }
  colStart[Ncols] = Nentries;

  return cMat;
}

// Double-valued vector
template <>
cholmod_dense* toCholmod(const Eigen::Matrix<double, Eigen::Dynamic, 1>& v, CholmodContext& context) {

  size_t N = v.rows();

  cholmod_dense* cVec = cholmod_l_allocate_dense(N, 1, N, CHOLMOD_REAL, context);
  double* cVecD = (double*)cVec->x;
  for (size_t i = 0; i < N; i++) {
    cVecD[i] = v(i);
  }

  return cVec;
}

// Float-valued vector
template <>
cholmod_dense* toCholmod(const Eigen::Matrix<float, Eigen::Dynamic, 1>& v, CholmodContext& context) {

  size_t N = v.rows();

  cholmod_dense* cVec = cholmod_l_allocate_dense(N, 1, N, CHOLMOD_REAL, context);
  double* cVecD = (double*)cVec->x;
  for (size_t i = 0; i < N; i++) {
    cVecD[i] = v(i);
  }

  return cVec;
}

// Complex-valued vector
template <>
cholmod_dense* toCholmod(const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>& v, CholmodContext& context) {

  size_t N = v.rows();

  cholmod_dense* cVec = cholmod_l_allocate_dense(N, 1, N, CHOLMOD_COMPLEX, context);
  std::complex<double>* cVecC = (std::complex<double>*)cVec->x;
  for (size_t i = 0; i < N; i++) {
    cVecC[i] = v(i);
  }

  return cVec;
}
// Convert a vector
template <typename T>
void toEigen(cholmod_dense* cVec, CholmodContext& context, Eigen::Matrix<T, Eigen::Dynamic, 1>& xOut) {

  if (cVec->ncol != 1) {
    throw std::logic_error("Input is not a vector");
  }
  size_t N = cVec->nrow;

  // Ensure output is large enough
  xOut = Eigen::Matrix<T, Eigen::Dynamic, 1>(N);

  // Type wizardry. This type is 'double' if T == 'float', and T otherwise
  // Needed because cholmod always uses double precision
  typedef typename std::conditional<std::is_same<T, float>::value, double, T>::type SCALAR_TYPE;

  SCALAR_TYPE* cVecS = (SCALAR_TYPE*)cVec->x;
  for (size_t i = 0; i < N; i++) {
    xOut(i) = cVecS[i];
  }
}
template void toEigen(cholmod_dense* cVec, CholmodContext& context, Eigen::Matrix<double, Eigen::Dynamic, 1>& xOut);
template void toEigen(cholmod_dense* cVec, CholmodContext& context, Eigen::Matrix<float, Eigen::Dynamic, 1>& xOut);
template void toEigen(cholmod_dense* cVec, CholmodContext& context,
                      Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>& xOut);

} // namespace geometrycentral
#endif
