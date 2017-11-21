namespace geometrycentral {
template <class T>
SparseFactorization<T>::SparseFactorization(SparseMatrix<T>& A)
    : ll(A), ldl(A), lu(A), qr(A) {}

template <class T>
void SparseFactorization<T>::clear(void) {
  ll.clear();
  ldl.clear();
  lu.clear();
  qr.clear();
}

template <class T>
void SparseFactorization<T>::clearNumeric(void) {
  ll.clearNumeric();
  ldl.clearNumeric();
  lu.clearNumeric();
  qr.clearNumeric();
}
}