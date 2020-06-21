#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include <memory>
#include <tuple>


// These methods construct the connectivity and geometry of a mesh simultaneously.

namespace geometrycentral {
namespace surface {

// Constructs a halfedge mesh and associated geometry from 0-indexed list of polygons and corresponding vertex
// positions.
//  - compressIndices: if true, will search the polygons for any unused vertices and re-index to exclude them when
//  constructing the mesh. This is necessary if there may be such unused elements, because the halfedge mesh cannot
//  understand/represent them.
std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
makeHalfedgeAndGeometry(const std::vector<std::vector<size_t>>& polygons, const std::vector<Vector3> vertexPositions);


// Like above, but with known twin connectivity
std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
makeHalfedgeAndGeometry(const std::vector<std::vector<size_t>>& polygons,
                        const std::vector<std::vector<std::tuple<size_t, size_t>>>& twins,
                        const std::vector<Vector3> vertexPositions);


// Same a above, but constructs a nonmanifold surface mesh
std::tuple<std::unique_ptr<SurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
makeNonmanifoldHalfedgeAndGeometry(const std::vector<std::vector<size_t>>& polygons,
                                   const std::vector<Vector3> vertexPositions);


// Like above, but with known twin connectivity
std::tuple<std::unique_ptr<SurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
makeNonmanifoldHalfedgeAndGeometry(const std::vector<std::vector<size_t>>& polygons,
                                   const std::vector<std::vector<std::tuple<size_t, size_t>>>& twins,
                                   const std::vector<Vector3> vertexPositions);

} // namespace surface
} // namespace geometrycentral
