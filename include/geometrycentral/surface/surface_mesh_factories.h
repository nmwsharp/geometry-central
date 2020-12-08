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
makeManifoldSurfaceMeshAndGeometry(const std::vector<std::vector<size_t>>& polygons,
                                   const std::vector<Vector3> vertexPositions);


// Like above, but with known twin connectivity (full, general version)
std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>,
           std::unique_ptr<CornerData<Vector2>>>
makeManifoldSurfaceMeshAndGeometry(const std::vector<std::vector<size_t>>& polygons,
                                   const std::vector<std::vector<std::tuple<size_t, size_t>>>& twins,
                                   const std::vector<Vector3> vertexPositions,
                                   const std::vector<std::vector<Vector2>>& paramCoordinates);

// Make a manifold surface mesh with both geometry and UV coordinates
std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>,
           std::unique_ptr<CornerData<Vector2>>>
makeParameterizedManifoldSurfaceMeshAndGeometry(const std::vector<std::vector<size_t>>& polygons,
                                                const std::vector<Vector3> vertexPositions,
                                                const std::vector<std::vector<Vector2>>& paramCoordinates);

// Same a above, but constructs a potentially-nonmanifold surface mesh
std::tuple<std::unique_ptr<SurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
makeSurfaceMeshAndGeometry(const std::vector<std::vector<size_t>>& polygons,
                           const std::vector<Vector3> vertexPositions);


// Like above, but with known twin connectivity (full, general version)
std::tuple<std::unique_ptr<SurfaceMesh>, std::unique_ptr<VertexPositionGeometry>, std::unique_ptr<CornerData<Vector2>>>
makeSurfaceMeshAndGeometry(const std::vector<std::vector<size_t>>& polygons,
                           const std::vector<std::vector<std::tuple<size_t, size_t>>>& twins,
                           const std::vector<Vector3> vertexPositions,
                           const std::vector<std::vector<Vector2>>& paramCoordinates);

// Like above, but with UV coordinates
std::tuple<std::unique_ptr<SurfaceMesh>, std::unique_ptr<VertexPositionGeometry>, std::unique_ptr<CornerData<Vector2>>>
makeParameterizedSurfaceMeshAndGeometry(const std::vector<std::vector<size_t>>& polygons,
                                        const std::vector<Vector3> vertexPositions,
                                        const std::vector<std::vector<Vector2>>& paramCoordinates);


// Make a manifold mesh from Eigen matrices
template <typename Scalar_V, typename Scalar_F>
std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
makeManifoldSurfaceMeshAndGeometry(const Eigen::MatrixBase<Scalar_V>& vMat, const Eigen::MatrixBase<Scalar_F>& fMat);

// Make a general mesh from Eigen matrices
template <typename Scalar_V, typename Scalar_F>
std::tuple<std::unique_ptr<SurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
makeSurfaceMeshAndGeometry(const Eigen::MatrixBase<Scalar_V>& vMat, const Eigen::MatrixBase<Scalar_F>& fMat);

} // namespace surface
} // namespace geometrycentral

#include "geometrycentral/surface/surface_mesh_factories.ipp"
