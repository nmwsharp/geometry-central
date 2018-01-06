#include "geometrycentral/linear_solvers.h"

#include "geometrycentral/linear_algebra_utilities.h"


using namespace Eigen;

namespace geometrycentral {

template <typename T>
Vector<T> LinearSolver<T>::operator()(const Vector<T>& rhs) {
  Vector<T> lhs;
  (*this)(lhs, rhs);
  return lhs;
}
template class LinearSolver<double>;
template class LinearSolver<float>;
template class LinearSolver<Complex>;


} // namespace geometrycentral