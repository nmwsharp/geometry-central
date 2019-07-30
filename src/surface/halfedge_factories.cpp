#include "geometrycentral/surface/halfedge_factories.h"


namespace geometrycentral {
namespace surface {
std::tuple<std::unique_ptr<HalfedgeMesh>, std::unique_ptr<VertexPositionGeometry>>
makeHalfedgeAndGeometry(const std::vector<std::vector<size_t>>& polygons, const std::vector<Vector3> vertexPositions,
                        bool compressIndices, bool verbose) {

  if (compressIndices) {

    // Check which indices are used
    size_t nV = vertexPositions.size();
    std::vector<char> vertexUsed(nV, false);
    for (auto poly : polygons) {
      for (auto i : poly) {
        GC_SAFETY_ASSERT(i < nV,
                         "polygon list has index " + std::to_string(i) + " >= num vertices " + std::to_string(nV));
        vertexUsed[i] = true;
      }
    }


    // Re-index
    std::vector<size_t> newInd(nV, INVALID_IND);
    std::vector<Vector3> newVertexPositions(nV);
    size_t nNewV = 0;
    for (size_t iOldV = 0; iOldV < nV; iOldV++) {
      if (!vertexUsed[iOldV]) continue;
      size_t iNewV = nNewV++;
      newInd[iOldV] = iNewV;
      newVertexPositions[iNewV] = vertexPositions[iOldV];
    }

    // Translate the polygon listing
    std::vector<std::vector<size_t>> newPolygons = polygons;
    for (auto& poly : newPolygons) {
      for (auto& i : poly) {
        i = newInd[i];
      }
    }


    // Construct
    std::unique_ptr<HalfedgeMesh> mesh(new HalfedgeMesh(newPolygons, verbose));
    std::unique_ptr<VertexPositionGeometry> geometry(new VertexPositionGeometry(*mesh));
    for (Vertex v : mesh->vertices()) {
      // Use the low-level indexers here since we're constructing
      (*geometry).inputVertexPositions[v] = newVertexPositions[v.getIndex()];
    }

    return std::make_tuple(std::move(mesh), std::move(geometry));

  } else {

    // Construct
    std::unique_ptr<HalfedgeMesh> mesh(new HalfedgeMesh(polygons, verbose));
    std::unique_ptr<VertexPositionGeometry> geometry(new VertexPositionGeometry(*mesh));
    for (Vertex v : mesh->vertices()) {
      // Use the low-level indexers here since we're constructing
      (*geometry).inputVertexPositions[v] = vertexPositions[v.getIndex()];
    }

    return std::make_tuple(std::move(mesh), std::move(geometry));
  }
}

} // namespace surface
} // namespace geometrycentral
