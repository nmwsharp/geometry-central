
namespace geometrycentral {
namespace surface {

inline double SignpostIntrinsicTriangulation::standardizeAngle(Vertex vert, double angle) const {
  if (vert.isBoundary()) {
    // can't wrap around at vertices
    return angle;
  }
  return std::fmod(angle, intrinsicVertexAngleSums[vert]);
}

inline Vector2 SignpostIntrinsicTriangulation::halfedgeVector(Halfedge he) const {
  double edgeAngle = intrinsicHalfedgeDirections[he];
  double scaleFac = 1.0 / vertexAngleScaling(he.vertex());
  Vector2 traceVec = Vector2::fromAngle(edgeAngle * scaleFac) * intrinsicEdgeLengths[he.edge()];
  return traceVec;
}

inline double SignpostIntrinsicTriangulation::vertexAngleScaling(Vertex v) const {
  return intrinsicVertexAngleSums[v] / (v.isBoundary() ? M_PI : 2. * M_PI);
}

inline double SignpostIntrinsicTriangulation::cornerAngle(Corner c) const {
  Halfedge heA = c.halfedge();
  Halfedge heOpp = heA.next();
  Halfedge heB = heOpp.next();

  double lOpp = intrinsicEdgeLengths[heOpp.edge()];
  double lA = intrinsicEdgeLengths[heA.edge()];
  double lB = intrinsicEdgeLengths[heB.edge()];

  double q = (lA * lA + lB * lB - lOpp * lOpp) / (2. * lA * lB);
  q = clamp(q, -1.0, 1.0);
  double angle = std::acos(q);

  return angle;
}

inline double SignpostIntrinsicTriangulation::halfedgeCotanWeight(Halfedge heI) const {
  if (heI.isInterior()) {
    Halfedge he = heI;
    double l_ij = intrinsicEdgeLengths[he.edge()];
    he = he.next();
    double l_jk = intrinsicEdgeLengths[he.edge()];
    he = he.next();
    double l_ki = intrinsicEdgeLengths[he.edge()];
    he = he.next();
    GC_SAFETY_ASSERT(he == heI, "faces mush be triangular");
    double areaV = area(he.face());
    double cotValue = (-l_ij * l_ij + l_jk * l_jk + l_ki * l_ki) / (4. * areaV);
    return cotValue / 2;
  } else {
    return 0.;
  }
}

inline double SignpostIntrinsicTriangulation::edgeCotanWeight(Edge e) const {
  return halfedgeCotanWeight(e.halfedge()) + halfedgeCotanWeight(e.halfedge().twin());
}

inline std::array<Vector2, 4> SignpostIntrinsicTriangulation::layoutDiamond(Halfedge iHe) {

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
  double l01 = intrinsicEdgeLengths[iHeA1.edge()];
  double l12 = intrinsicEdgeLengths[iHeA2.edge()];
  double l23 = intrinsicEdgeLengths[iHeB1.edge()];
  double l30 = intrinsicEdgeLengths[iHeB2.edge()];
  double l02 = intrinsicEdgeLengths[iHeA0.edge()];

  // Lay out the vertices of the diamond
  Vector2 p3{0., 0.};
  Vector2 p0{l30, 0.};
  Vector2 p2 = layoutTriangleVertex(p3, p0, l02, l23); // involves more arithmetic than strictly necessary
  Vector2 p1 = layoutTriangleVertex(p2, p0, l01, l12);

  return {p0, p1, p2, p3};
}

inline double SignpostIntrinsicTriangulation::shortestEdge(Face f) const {
  Halfedge he = f.halfedge();
  double lA = intrinsicEdgeLengths[he.edge()];
  he = he.next();
  double lB = intrinsicEdgeLengths[he.edge()];
  he = he.next();
  double lC = intrinsicEdgeLengths[he.edge()];
  return std::fmin(std::fmin(lA, lB), lC);
}

inline double SignpostIntrinsicTriangulation::area(Face f) const {
  Halfedge he = f.halfedge();
  double a = intrinsicEdgeLengths[he.edge()];
  he = he.next();
  double b = intrinsicEdgeLengths[he.edge()];
  he = he.next();
  double c = intrinsicEdgeLengths[he.edge()];
  return triangleArea(a, b, c);
}

inline double SignpostIntrinsicTriangulation::circumradius(Face f) const {
  Halfedge he = f.halfedge();
  double a = intrinsicEdgeLengths[he.edge()];
  he = he.next();
  double b = intrinsicEdgeLengths[he.edge()];
  he = he.next();
  double c = intrinsicEdgeLengths[he.edge()];

  double A = triangleArea(a, b, c);
  return a * b * c / (4. * A);
}

inline std::array<Vector2, 3> SignpostIntrinsicTriangulation::vertexCoordinatesInTriangle(Face face) {
  return {Vector2{0., 0.}, halfedgeVectorsInFace[face.halfedge()],
          -halfedgeVectorsInFace[face.halfedge().next().next()]};
}

inline bool SignpostIntrinsicTriangulation::isFixed(Edge e) {
  if (e.isBoundary()) return true;
  if (markedEdges.size() > 0 && markedEdges[e]) return true;
  return false;
}

inline bool SignpostIntrinsicTriangulation::isOnFixedEdge(Vertex v) {
  for (Edge e : v.adjacentEdges()) {
    if (isFixed(e)) return true;
  }
  return false;
}


template <typename T>
VertexData<T> SignpostIntrinsicTriangulation::sampleFromInput(const VertexData<T>& dataOnInput) {
  VertexData<T> output(mesh);
  for (Vertex v : mesh.vertices()) {
    output[v] = vertexLocations[v].interpolate(dataOnInput);
  }
  return output;
}

template <typename T>
VertexData<T> SignpostIntrinsicTriangulation::restrictToInput(const VertexData<T>& dataOnIntrinsic) {
  VertexData<T> output(inputMesh);
  for (Vertex v : mesh.vertices()) {
    if(vertexLocations[v].type == SurfacePointType::Vertex) {
      output[vertexLocations[v].vertex] = dataOnIntrinsic[v];
    }
  }
  return output;
}

} // namespace surface
} // namespace geometrycentral
