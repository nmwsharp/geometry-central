#include "geometrycentral/surface/surface_mesh_factories.h"


namespace geometrycentral {
namespace surface {


std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
makeManifoldSurfaceMeshAndGeometry(const std::vector<std::vector<size_t>>& polygons,
                                   const std::vector<Vector3> vertexPositions) {
  auto lvals = makeManifoldSurfaceMeshAndGeometry(polygons, {}, vertexPositions, {});

  return std::tuple<std::unique_ptr<ManifoldSurfaceMesh>,
                    std::unique_ptr<VertexPositionGeometry>>(std::move(std::get<0>(lvals)),  // mesh
                                                             std::move(std::get<1>(lvals))); // geometry
}

std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>,
           std::unique_ptr<CornerData<Vector2>>>
makeManifoldSurfaceMeshAndGeometry(const std::vector<std::vector<size_t>>& polygons,
                                   const std::vector<std::vector<std::tuple<size_t, size_t>>>& twins,
                                   const std::vector<Vector3> vertexPositions,
                                   const std::vector<std::vector<Vector2>>& paramCoordinates) {

  // Construct
  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  if (twins.empty()) {
    mesh.reset(new ManifoldSurfaceMesh(polygons));
  } else {
    mesh.reset(new ManifoldSurfaceMesh(polygons, twins));
  }
  std::unique_ptr<VertexPositionGeometry> geometry(new VertexPositionGeometry(*mesh));
  for (Vertex v : mesh->vertices()) {
    // Use the low-level indexers here since we're constructing
    (*geometry).vertexPositions[v] = vertexPositions[v.getIndex()];
  }
  std::unique_ptr<CornerData<Vector2>> parameterization(new CornerData<Vector2>(*mesh));
  if (paramCoordinates.size() == mesh->nFaces()) {
    for (size_t i = 0; i < mesh->nFaces(); i++) {
      Halfedge h = mesh->face(i).halfedge();
      for (size_t j = 0; j < paramCoordinates[i].size(); j++) {
        (*parameterization)[h.corner()] = paramCoordinates[i][j];
        h = h.next();
      }
    }
  }

  return std::make_tuple(std::move(mesh), std::move(geometry), std::move(parameterization));
}

std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>,
           std::unique_ptr<CornerData<Vector2>>>
makeParameterizedManifoldSurfaceMeshAndGeometry(const std::vector<std::vector<size_t>>& polygons,
                                                const std::vector<Vector3> vertexPositions,
                                                const std::vector<std::vector<Vector2>>& paramCoordinates) {

  return makeManifoldSurfaceMeshAndGeometry(polygons, {}, vertexPositions, paramCoordinates);
}

std::tuple<std::unique_ptr<SurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
makeSurfaceMeshAndGeometry(const std::vector<std::vector<size_t>>& polygons,
                           const std::vector<Vector3> vertexPositions) {

  auto lvals = makeSurfaceMeshAndGeometry(polygons, {}, vertexPositions, {});
  return std::tuple<std::unique_ptr<SurfaceMesh>,
                    std::unique_ptr<VertexPositionGeometry>>(std::move(std::get<0>(lvals)),  // mesh
                                                             std::move(std::get<1>(lvals))); // geometry
}

// Like above, but with UV coordinates
std::tuple<std::unique_ptr<SurfaceMesh>, std::unique_ptr<VertexPositionGeometry>, std::unique_ptr<CornerData<Vector2>>>
makeParameterizedSurfaceMeshAndGeometry(const std::vector<std::vector<size_t>>& polygons,
                                        const std::vector<Vector3> vertexPositions,
                                        const std::vector<std::vector<Vector2>>& paramCoordinates) {

  return makeSurfaceMeshAndGeometry(polygons, {}, vertexPositions, paramCoordinates);
}


std::tuple<std::unique_ptr<SurfaceMesh>, std::unique_ptr<VertexPositionGeometry>, std::unique_ptr<CornerData<Vector2>>>
makeSurfaceMeshAndGeometry(const std::vector<std::vector<size_t>>& polygons,
                           const std::vector<std::vector<std::tuple<size_t, size_t>>>& twins,
                           const std::vector<Vector3> vertexPositions,
                           const std::vector<std::vector<Vector2>>& paramCoordinates) {

  // Construct
  std::unique_ptr<SurfaceMesh> mesh;
  if (twins.empty()) {
    mesh.reset(new SurfaceMesh(polygons));
  } else {
    mesh.reset(new SurfaceMesh(polygons, twins));
  }
  std::unique_ptr<VertexPositionGeometry> geometry(new VertexPositionGeometry(*mesh));
  for (Vertex v : mesh->vertices()) {
    // Use the low-level indexers here since we're constructing
    (*geometry).vertexPositions[v] = vertexPositions[v.getIndex()];
  }

  std::unique_ptr<CornerData<Vector2>> parameterization(new CornerData<Vector2>(*mesh));
  if (paramCoordinates.size() == mesh->nFaces()) {
    for (size_t i = 0; i < mesh->nFaces(); i++) {
      Halfedge h = mesh->face(i).halfedge();
      for (size_t j = 0; j < paramCoordinates[i].size(); j++) {
        (*parameterization)[h.corner()] = paramCoordinates[i][j];
        h = h.next();
      }
    }
  }

  return std::make_tuple(std::move(mesh), std::move(geometry), std::move(parameterization));
}


} // namespace surface
} // namespace geometrycentral
