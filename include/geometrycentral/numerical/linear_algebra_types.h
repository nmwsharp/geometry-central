#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace geometrycentral {

// Convenience typedefs

// Nicer name for dynamic matrix
template <typename T>
using Vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;

// Nicer name for sparse matrix
template <typename T>
using SparseMatrix = Eigen::SparseMatrix<T>;

// Nicer name for dense matrix
template <typename T>
using DenseMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

} // namespace geometrycentral
