#pragma once

#include "geometrycentral/numerical/linear_algebra_utilities.h"

namespace geometrycentral {
namespace surface {

template <typename T>
ManifoldSurfaceMesh::ManifoldSurfaceMesh(const Eigen::MatrixBase<T>& faces)
    : ManifoldSurfaceMesh(unpackMatrixToStdVector<size_t>(faces.template cast<size_t>())) {}


} // namespace surface
} // namespace geometrycentral
