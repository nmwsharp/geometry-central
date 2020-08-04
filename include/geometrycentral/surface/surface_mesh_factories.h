#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include <memory>
#include <tuple>


// These methods construct the connectivity and geometry of a mesh simultaneously.

namespace geometrycentral {
namespace surface {

// Assumes manifoldness, errors our if not
std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
makeManifoldSurfaceMeshAndGeometry(const std::vector<std::vector<size_t>>& polygons, const std::vector<Vector3> vertexPositions);


// Like above, but with known twin connectivity
std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
makeManifoldSurfaceMeshAndGeometry(const std::vector<std::vector<size_t>>& polygons,
                        const std::vector<std::vector<std::tuple<size_t, size_t>>>& twins,
                        const std::vector<Vector3> vertexPositions);


// Same a above, but constructs a potentially-nonmanifold surface mesh
std::tuple<std::unique_ptr<SurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
makeSurfaceMeshAndGeometry(const std::vector<std::vector<size_t>>& polygons,
                               const std::vector<Vector3> vertexPositions);


// Like above, but with known twin connectivity
std::tuple<std::unique_ptr<SurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
makeSurfaceMeshAndGeometry(const std::vector<std::vector<size_t>>& polygons,
                               const std::vector<std::vector<std::tuple<size_t, size_t>>>& twins,
                               const std::vector<Vector3> vertexPositions);


// Make a manifold mesh from Eigen matrices
template <typename Scalar_V, typename Scalar_F>
std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
makeManifoldSurfaceMeshAndGeometry(const Eigen::MatrixBase<Scalar_V>& vMat, const Eigen::MatrixBase<Scalar_F>& fMat);

// Make a general mesh from Eigen matrices
template <typename Scalar_V, typename Scalar_F>
std::tuple<std::unique_ptr<SurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
makeSurfaceMeshAndGeometry(const Eigen::MatrixBase<Scalar_V>& vMat, const Eigen::MatrixBase<Scalar_F>& fMat);

}}

#include "geometrycentral/surface/surface_mesh_factories.ipp"
