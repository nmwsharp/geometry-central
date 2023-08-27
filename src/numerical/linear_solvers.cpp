#include "geometrycentral/numerical/linear_solvers.h"

#include "geometrycentral/numerical/linear_algebra_utilities.h"
#include "geometrycentral/utilities/vector2.h"


namespace geometrycentral {

template class LinearSolver<double>;
template class LinearSolver<float>;
template class LinearSolver<std::complex<double>>;


template <typename T>
double residual(const SparseMatrix<T>& matrix, const Vector<T>& lhs, const Vector<T>& rhs) {
  Vector<T> residVec = matrix * lhs - rhs;
  double resid = std::abs((residVec.conjugate().transpose() * residVec)(0));
  return std::sqrt(resid);
}


template double residual(const SparseMatrix<float>& matrix, const Vector<float>& lhs, const Vector<float>& rhs);
template double residual(const SparseMatrix<double>& matrix, const Vector<double>& lhs, const Vector<double>& rhs);
template double residual(const SparseMatrix<std::complex<double>>& matrix, const Vector<std::complex<double>>& lhs,
                         const Vector<std::complex<double>>& rhs);

} // namespace geometrycentral
