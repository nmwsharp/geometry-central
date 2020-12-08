#include "geometrycentral/surface/halfedge_factories.h"

// NOTE these are DEPRECATED, and exist only for compatability. Prefer the versions from surface_mesh_factories.h

namespace geometrycentral {
namespace surface {


// (these just forward to the newly-named versions)

std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
makeHalfedgeAndGeometry(const std::vector<std::vector<size_t>>& polygons, const std::vector<Vector3> vertexPositions) {
  return makeHalfedgeAndGeometry(polygons, {}, vertexPositions);
}

std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
makeHalfedgeAndGeometry(const std::vector<std::vector<size_t>>& polygons,
                        const std::vector<std::vector<std::tuple<size_t, size_t>>>& twins,
                        const std::vector<Vector3> vertexPositions) {
  auto lvals = makeManifoldSurfaceMeshAndGeometry(polygons, twins, vertexPositions, {});

  return std::tuple<std::unique_ptr<ManifoldSurfaceMesh>,
                    std::unique_ptr<VertexPositionGeometry>>(std::move(std::get<0>(lvals)),  // mesh
                                                             std::move(std::get<1>(lvals))); // geometry
}


std::tuple<std::unique_ptr<SurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
makeGeneralHalfedgeAndGeometry(const std::vector<std::vector<size_t>>& polygons,
                               const std::vector<Vector3> vertexPositions) {
  return makeGeneralHalfedgeAndGeometry(polygons, {}, vertexPositions);
}


std::tuple<std::unique_ptr<SurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
makeGeneralHalfedgeAndGeometry(const std::vector<std::vector<size_t>>& polygons,
                               const std::vector<std::vector<std::tuple<size_t, size_t>>>& twins,
                               const std::vector<Vector3> vertexPositions) {

  auto lvals = makeSurfaceMeshAndGeometry(polygons, twins, vertexPositions, {});

  return std::tuple<std::unique_ptr<SurfaceMesh>,
                    std::unique_ptr<VertexPositionGeometry>>(std::move(std::get<0>(lvals)),  // mesh
                                                             std::move(std::get<1>(lvals))); // geometry
}


} // namespace surface
} // namespace geometrycentral
