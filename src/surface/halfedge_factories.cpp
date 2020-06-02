#include "geometrycentral/surface/halfedge_factories.h"


namespace geometrycentral {
namespace surface {


std::tuple<std::unique_ptr<HalfedgeMesh>, std::unique_ptr<VertexPositionGeometry>>
makeHalfedgeAndGeometry(const std::vector<std::vector<size_t>>& polygons, const std::vector<Vector3> vertexPositions) {

  return makeHalfedgeAndGeometry(polygons, {}, vertexPositions);
}

std::tuple<std::unique_ptr<HalfedgeMesh>, std::unique_ptr<VertexPositionGeometry>>
makeHalfedgeAndGeometry(const std::vector<std::vector<size_t>>& polygons,
                        const std::vector<std::vector<std::tuple<size_t, size_t>>>& twins,
                        const std::vector<Vector3> vertexPositions) {

  // Construct
  std::unique_ptr<HalfedgeMesh> mesh;
  if (twins.empty()) {
    mesh.reset(new HalfedgeMesh(polygons));
  } else {
    mesh.reset(new HalfedgeMesh(polygons, twins, true));
  }
  std::unique_ptr<VertexPositionGeometry> geometry(new VertexPositionGeometry(*mesh));
  for (Vertex v : mesh->vertices()) {
    // Use the low-level indexers here since we're constructing
    (*geometry).inputVertexPositions[v] = vertexPositions[v.getIndex()];
  }

  return std::make_tuple(std::move(mesh), std::move(geometry));
}

} // namespace surface
} // namespace geometrycentral
