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


/// Type alias for aligned std::vectors
// template <typename T>
// using AlignedVector_T = std::vector<T, Eigen::aligned_allocator<T>>;

/**
 * @brief Typename for an Eigen map to an raw buffer
 *
 * @tparam T            Typename of the contained data
 * @tparam k            Number of columns
 * @tparam Options      Storage order \b Eigen::RowMajor or \b Eigen::ColMajor
 * @tparam Alignment    Alignment of the underlying data
 */
template <typename T, std::size_t k, int Options = Eigen::ColMajor, int Alignment = Eigen::AlignedMax>
using EigenVectorMap_T = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, k, Options>, Alignment>;

/**
 * @brief Typename for an Eigen map to a const raw buffer
 *
 * @tparam T            Typename of the contained data
 * @tparam k            Number of columns
 * @tparam Options      Storage order \b Eigen::RowMajor or \b Eigen::ColMajor
 * @tparam Alignment    Alignment of the underlying data
 */
template <typename T, std::size_t k, int Options = Eigen::ColMajor, int Alignment = Eigen::AlignedMax>
using ConstEigenVectorMap_T = Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, k, Options>, Alignment>;

} // namespace geometrycentral
