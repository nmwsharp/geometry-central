#include "geometrycentral/surface/parameterize.h"

#include "geometrycentral/numerical/linear_solvers.h"
#include "geometrycentral/surface/edge_length_geometry.h"
#include "geometrycentral/surface/uniformize.h"
#include "geometrycentral/utilities/elementary_geometry.h"

#include <algorithm>
#include <queue>
#include <unordered_map>
#include <utility>
#include <vector>

namespace geometrycentral {
namespace surface {

VertexData<Vector2> parameterizeDisk(ManifoldSurfaceMesh& origMesh, IntrinsicGeometryInterface& origGeom) {
  // Check that it's a (punctured) disk

  /*
   if ((long long int)origMesh.eulerCharacteristic() - (long long int)(2 * origMesh.nBoundaryLoops()) != -2) {
    long long int val =
        ((long long int)origMesh.eulerCharacteristic() - (long long int)(2 * origMesh.nBoundaryLoops()));
    throw std::runtime_error("parameterizeDisk(): input origMesh must be a (possibly punctured) disk, chi - 2b = " +
                             std::to_string(val));
  }
  */

  // Get uniformized edge lengths
  // Copy the mesh, since we will flip its edges
  std::unique_ptr<ManifoldSurfaceMesh> meshPtr = origMesh.copy();
  ManifoldSurfaceMesh& mesh = *meshPtr;
  origGeom.requireEdgeLengths();
  EdgeData<double> copyLens = origGeom.edgeLengths.reinterpretTo(mesh);
  EdgeLengthGeometry geometry(mesh, copyLens);
  origGeom.unrequireEdgeLengths();
  EdgeData<double> uLens = uniformizeDisk(mesh, geometry, true);

  // Layout
  VertexData<Vector2> coords(mesh);
  VertexData<char> haveCoords(mesh, false);

  // == Layout until finished

  // <n_neighbors, -n_round, vert>
  // prefer verts with more neighbors, breaking ties by earlier verts
  using WeightedVertex = std::tuple<int, int, Vertex>;
  std::priority_queue<WeightedVertex> pq;
  int layoutRound = 0;

  // Helper to add vertices to the layout queue if they have enough neighbors to layout
  auto considerVertex = [&](Vertex v) {
    if (haveCoords[v]) return;

    int neighCount = 0;
    for (Halfedge he : v.outgoingHalfedges()) {
      if (!he.isInterior()) continue;
      Halfedge heOpp = he.next();
      if (haveCoords[heOpp.vertex()] && haveCoords[heOpp.twin().vertex()]) neighCount++;
    }

    if (neighCount >= 1) pq.emplace(neighCount, -layoutRound, v);
  };

  // initialize the layout
  Halfedge he0 = mesh.halfedge(0);
  Vertex v0 = he0.vertex();
  Vertex v1 = he0.twin().vertex();
  coords[v0] = Vector2{0., 0.};
  coords[v1] = Vector2{uLens[he0.edge()], 0.};
  haveCoords[v0] = true;
  haveCoords[v1] = true;
  considerVertex(he0.next().next().vertex());
  considerVertex(he0.twin().next().next().vertex());

  // layout until finished
  while (!pq.empty()) {

    int topCount = std::get<0>(pq.top());
    Vertex topVert = std::get<2>(pq.top());
    pq.pop();

    if (haveCoords[topVert]) continue;

    // std::cout << "laying out " << topVert << " with " << topCount << " neighbors" << std::endl;

    // Compute a new position for the vertex, as an average of laid out position from all neighbors
    Vector2 avgPos{0., 0.};
    int avgCount = 0;
    for (Halfedge he : topVert.outgoingHalfedges()) {
      if (!he.isInterior()) continue;
      Halfedge heOpp = he.next();
      Vertex vA = heOpp.vertex();
      Vertex vB = heOpp.twin().vertex();
      if (haveCoords[vA] && haveCoords[vB]) {

        Vector2 pA = coords[vA];
        Vector2 pB = coords[vB];
        double lCA = uLens[he.edge()];
        double lBC = uLens[heOpp.next().edge()];

        Vector2 newPos = layoutTriangleVertex(pA, pB, lBC, lCA);
        avgPos += newPos;
        avgCount++;
      }
    }

    // set the new positions
    coords[topVert] = avgPos / avgCount;
    haveCoords[topVert] = true;

    // add neighbors
    for (Vertex vn : topVert.adjacentVertices()) {
      considerVertex(vn);
    }

    layoutRound++;
  }

  return coords.reinterpretTo(origMesh);
}

} // namespace surface
} // namespace geometrycentral
