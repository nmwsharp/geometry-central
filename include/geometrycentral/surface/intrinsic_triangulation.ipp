#pragma once

#include "geometrycentral/surface/intrinsic_triangulation.h"

namespace geometrycentral {
namespace surface {

inline bool IntrinsicTriangulation::isFixed(Edge e) const {
  if (e.isBoundary()) return true;
  if (markedEdges.size() > 0 && markedEdges[e]) return true;
  return false;
}

inline bool IntrinsicTriangulation::isOnFixedEdge(Vertex v) const {
  for (Edge e : v.adjacentEdges()) {
    if (isFixed(e)) return true;
  }
  return false;
}


template <typename T>
VertexData<T> IntrinsicTriangulation::sampleFromInput(const VertexData<T>& dataOnInput) {
  VertexData<T> output(mesh);
  for (Vertex v : mesh.vertices()) {
    output[v] = vertexLocations[v].interpolate(dataOnInput);
  }
  return output;
}

template <typename T>
VertexData<T> IntrinsicTriangulation::restrictToInput(const VertexData<T>& dataOnIntrinsic) {
  VertexData<T> output(inputMesh);
  for (Vertex v : mesh.vertices()) {
    if (vertexLocations[v].type == SurfacePointType::Vertex) {
      output[vertexLocations[v].vertex] = dataOnIntrinsic[v];
    }
  }
  return output;
}

inline std::array<Vector2, 4> IntrinsicTriangulation::layoutDiamond(Halfedge iHe) {

  // Conventions:
  //  - iHe points from vertex 2 to vertex 0, other vertices are numbered ccw
  //  - iHe is incident on face A, other is face B
  //  - halfedges within face are numbered CCW as A0, A1, A2 (etc),
  //    starting with iHe and twin(iHe)
  //  - When we lay out the triangle, p3 is at the origin and
  //    edge 3-0 is along the X-axis
  //  - flips is always ccw, so iHe points from vertex 3 --> 1 after

  // Gather index values
  Halfedge iHeA0 = iHe;
  Halfedge iHeA1 = iHeA0.next();
  Halfedge iHeA2 = iHeA1.next();
  Halfedge iHeB0 = iHe.twin();
  Halfedge iHeB1 = iHeB0.next();
  Halfedge iHeB2 = iHeB1.next();

  // Gather length values
  double l01 = edgeLengths[iHeA1.edge()];
  double l12 = edgeLengths[iHeA2.edge()];
  double l23 = edgeLengths[iHeB1.edge()];
  double l30 = edgeLengths[iHeB2.edge()];
  double l02 = edgeLengths[iHeA0.edge()];

  // Lay out the vertices of the diamond
  Vector2 p3{0., 0.};
  Vector2 p0{l30, 0.};
  Vector2 p2 = layoutTriangleVertex(p3, p0, l02, l23); // involves more arithmetic than strictly necessary
  Vector2 p1 = layoutTriangleVertex(p2, p0, l01, l12);

  return {p0, p1, p2, p3};
}

inline double IntrinsicTriangulation::shortestEdge(Face f) const {
  Halfedge he = f.halfedge();
  double lA = edgeLengths[he.edge()];
  he = he.next();
  double lB = edgeLengths[he.edge()];
  he = he.next();
  double lC = edgeLengths[he.edge()];
  return std::fmin(std::fmin(lA, lB), lC);
}

inline std::array<Vector2, 3> IntrinsicTriangulation::vertexCoordinatesInTriangle(Face face) {
  return {Vector2{0., 0.}, halfedgeVectorsInFace[face.halfedge()],
          -halfedgeVectorsInFace[face.halfedge().next().next()]};
}

} // namespace surface
} // namespace geometrycentral
