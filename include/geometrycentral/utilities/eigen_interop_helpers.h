/**
 * @file eigen_interop_helpers.h
 * @brief Utility functions primarily for interfacing with Eigen and raw buffers
 *
 */
#pragma once

#include <Eigen/Core>

#include "geometrycentral/numerical/linear_algebra_types.h"
#include "geometrycentral/surface/surface_mesh.h"

namespace geometrycentral {

using namespace surface;

template <typename E, typename O>
using EigenTraits = Eigen::internal::traits<typename MeshData<E, O>::DATA_T>;

/// Type alias for aligned vectors
template <typename T>
using AlignedVector_T = std::vector<T, Eigen::aligned_allocator<T>>;

/**
 * @brief Generate an Eigen Map to an raw buffer.

 *
 * @tparam T          Typename of the output data
 * @tparam k          Number of columns in the output
 * @tparam Options    Storage order \b Eigen::RowMajor or \b Eigen::ColMajor
 * @tparam Alignment  Data alignment
 * @param vec         The data
 * @return EigenVectorMap_T<T, k, Options, Alignment>    Eigen map object to the buffer
 */
template <typename T, std::size_t k, int Options = Eigen::RowMajor, int Alignment>
EigenVectorMap_T<T, k, Options, Alignment> EigenMap(T* vec, std::size_t size) {
  return EigenVectorMap_T<T, k, Options, Alignment>(vec, size, k);
}

/**
 * @brief Generate an Eigen Map to an raw buffer.

 *
 * @tparam T          Typename of the output data
 * @tparam k          Number of columns in the output
 * @tparam Options    Storage order \b Eigen::RowMajor or \b Eigen::ColMajor
 * @tparam Alignment  Data alignment
 * @param vec         The data
 * @return EigenVectorMap_T<T, k, Options, Alignment>    Eigen map object to the buffer
 */
template <typename T, std::size_t k, int Options = Eigen::RowMajor, int Alignment>
ConstEigenVectorMap_T<T, k, Options, Alignment> EigenMap(const T* vec, std::size_t size) {
  return ConstEigenVectorMap_T<T, k, Options, Alignment>(vec, size, k);
}

/**
 * @brief Generate an Eigen Map to a raw buffer.
 *
 * This function will force cast the data from O to type T if O is a POD type.
 *
 * @tparam T        Typename of the output data
 * @tparam k        Number of columns in the output
 * @tparam Options  Storage order \b Eigen::RowMajor or \b Eigen::ColMajor
 * @tparam E        Typename of the mesh element
 * @tparam O        Typename of the input data

 * @param vec
 * @return EigenVectorMap_T<T, k, Options, Alignment>
 */
template <typename T, std::size_t k, int Options = Eigen::RowMajor, typename E, typename O>
auto EigenMap(MeshData<E, O>& vec) -> EigenVectorMap_T<T, k, Options, EigenTraits<E, O>::Alignment> {
  static_assert(std::is_standard_layout<T>::value && std::is_trivially_copyable<O>::value, "O must be a POD type.");
  static_assert(sizeof(O) == k * sizeof(T), "sizeof(O) must be a k multiple of sizeof(T)");
  return EigenMap<T, k, Options, EigenTraits<E, O>::Alignment>(reinterpret_cast<T*>(vec.raw().data()), vec.size());
}


/**
 * @brief Generate an Eigen Map to an aligned raw buffer.
 *
 * This function will force cast the data from O to type T if O is a POD type.
 *
 * @tparam T        Typename of the output data
 * @tparam k        Number of columns in the output
 * @tparam Options  Storage order \b Eigen::RowMajor or \b Eigen::ColMajor
 * @tparam E        Typename of the mesh element
 * @tparam O        Typename of the input data
 * @param vec       The data
 * @return ConstEigenVectorMap_T<T, k, Options, Alignment>
 */
template <typename T, std::size_t k, int Options = Eigen::RowMajor, typename E, typename O>
auto EigenMap(const MeshData<E, O>& vec) -> ConstEigenVectorMap_T<T, k, Options, EigenTraits<E, O>::Alignment> {
  static_assert(std::is_standard_layout<T>::value && std::is_trivially_copyable<O>::value, "O must be a POD type.");
  static_assert(sizeof(O) == k * sizeof(T), "sizeof(O) must be a k multiple of sizeof(T)");
  return EigenMap<T, k, Options, EigenTraits<E, O>::Alignment>(reinterpret_cast<const T*>(vec.raw().data()),
                                                               vec.size());
}


/**
 * @brief Generate an Eigen Map to an aligned raw buffer.
 *
 * This function will force cast the data from O to type T if O is a POD type.
 *
 * @tparam T    Typename of the output data
 * @tparam k    Number of columns in the output
 * @tparam E    Typename of the mesh element
 * @tparam O    Typename of the input data
 * @param vec   The data
 * @return EigenVectorMap_T<T, k, Options, Alignment>    Eigen map object to the buffer
 */
template <typename T, std::size_t k, typename E, typename O>
auto FlattenedEigenMap(MeshData<E, O>& vec) -> EigenVectorMap_T<T, 1, Eigen::ColMajor, EigenTraits<E, O>::Alignment> {
  using Return_T = EigenVectorMap_T<T, 1, Eigen::ColMajor, EigenTraits<E, O>::Alignment>;
  static_assert(std::is_standard_layout<T>::value && std::is_trivially_copyable<O>::value, "O must be a POD type.");
  static_assert(sizeof(O) == k * sizeof(T), "sizeof(O) must be a k multiple of sizeof(T)");
  return Return_T(reinterpret_cast<T*>(vec.raw().data()), k * vec.size());
}

/**
 * @brief Generate an Eigen Map to an aligned raw buffer.
 *
 * This function will force cast the data from O to type T if O is a POD type.
 *
 * @tparam T    Typename of the output data
 * @tparam k    Number of columns in the output
 * @tparam E    Typename of the mesh element
 * @tparam O    Typename of the input data
 * @param vec   The data
 * @return ConstEigenVectorMap_T<T, k, Options, Alignment>    Eigen map object to the buffer
 */
template <typename T, std::size_t k, typename E, typename O>
auto FlattenedEigenMap(const MeshData<E, O>& vec)
    -> ConstEigenVectorMap_T<T, 1, Eigen::ColMajor, EigenTraits<E, O>::Alignment> {
  using Return_T = ConstEigenVectorMap_T<T, 1, Eigen::ColMajor, EigenTraits<E, O>::Alignment>;
  static_assert(std::is_standard_layout<T>::value && std::is_trivially_copyable<O>::value, "O must be a POD type.");
  static_assert(sizeof(O) == k * sizeof(T), "sizeof(O) must be a k multiple of sizeof(T)");
  return Return_T(reinterpret_cast<const T*>(vec.raw().data()), k * vec.size());
}
} // namespace geometrycentral
