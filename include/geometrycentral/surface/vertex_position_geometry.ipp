namespace geometrycentral {
namespace surface {


template <typename T>
VertexPositionGeometry::VertexPositionGeometry(SurfaceMesh& mesh_, const Eigen::MatrixBase<T>& vMat)
    : VertexPositionGeometry(mesh_) {

  // sanity checks on input matrix dimensions
  GC_SAFETY_ASSERT(vMat.cols() == 3, "input must be a V x 3 matrix -- cols() should == 3");
  GC_SAFETY_ASSERT(static_cast<size_t>(vMat.rows()) == mesh_.nVertices(),
                   "input must be a V x 3 matrix -- rows() should == nVertices()");

  size_t iV = 0;
  for (Vertex v : mesh_.vertices()) {
    double x = vMat(iV, 0);
    double y = vMat(iV, 1);
    double z = vMat(iV, 2);
    inputVertexPositions[v] = Vector3{x, y, z};
    iV++;
  }
}

inline double VertexPositionGeometry::edgeLength(Edge e) const {
  Halfedge he = e.halfedge();
  Vector3 pA = inputVertexPositions[he.vertex()];
  Vector3 pB = inputVertexPositions[he.next().vertex()];
  return norm(pA - pB);
}

// Face areas
inline double VertexPositionGeometry::faceArea(Face f) const {
  // WARNING: Logic duplicated between cached and immediate version
  Halfedge he = f.halfedge();
  Vector3 pA = inputVertexPositions[he.vertex()];
  he = he.next();
  Vector3 pB = inputVertexPositions[he.vertex()];
  he = he.next();
  Vector3 pC = inputVertexPositions[he.vertex()];

  GC_SAFETY_ASSERT(he.next() == f.halfedge(), "faces mush be triangular");

  double area = 0.5 * norm(cross(pB - pA, pC - pA));
  return area;
}

// Corner angles
inline double VertexPositionGeometry::cornerAngle(Corner c) const {
  // WARNING: Logic duplicated between cached and immediate version
  Halfedge he = c.halfedge();
  Vector3 pA = inputVertexPositions[he.vertex()];
  he = he.next();
  Vector3 pB = inputVertexPositions[he.vertex()];
  he = he.next();
  Vector3 pC = inputVertexPositions[he.vertex()];

  GC_SAFETY_ASSERT(he.next() == c.halfedge(), "faces mush be triangular");

  double q = dot(unit(pB - pA), unit(pC - pA));
  q = clamp(q, -1.0, 1.0);
  double angle = std::acos(q);
  return angle;
}

inline double VertexPositionGeometry::halfedgeCotanWeight(Halfedge heI) const {
  // WARNING: Logic duplicated between cached and immediate version
  if (heI.isInterior()) {
    Halfedge he = heI;
    Vector3 pB = inputVertexPositions[he.vertex()];
    he = he.next();
    Vector3 pC = inputVertexPositions[he.vertex()];
    he = he.next();
    Vector3 pA = inputVertexPositions[he.vertex()];
    GC_SAFETY_ASSERT(he.next() == heI, "faces mush be triangular");

    Vector3 vecR = pB - pA;
    Vector3 vecL = pC - pA;

    double cotValue = dot(vecR, vecL) / norm(cross(vecR, vecL));
    return cotValue / 2;
  } else {
    return 0.;
  }
}

inline double VertexPositionGeometry::edgeCotanWeight(Edge e) const {
  double sum = 0;
  for (Halfedge he : e.adjacentInteriorHalfedges()) {
    sum += halfedgeCotanWeight(he);
  }
  return sum;
}

// Face normal
inline Vector3 VertexPositionGeometry::faceNormal(Face f) const {
  // For general polygons, take the sum of the cross products at each corner
  Vector3 normalSum = Vector3::zero();
  for (Halfedge heF : f.adjacentHalfedges()) {

    // Gather vertex positions for next three vertices
    Halfedge he = heF;
    Vector3 pA = inputVertexPositions[he.vertex()];
    he = he.next();
    Vector3 pB = inputVertexPositions[he.vertex()];
    he = he.next();
    Vector3 pC = inputVertexPositions[he.vertex()];

    normalSum += cross(pB - pA, pC - pA);

    // In the special case of a triangle, there is no need to to repeat at all three corners; the result will be the
    // same
    if (he.next() == heF) break;
  }

  Vector3 normal = unit(normalSum);
  return normal;
}
  

inline Vector3 VertexPositionGeometry::halfedgeVector(Halfedge he) const {
  return inputVertexPositions[he.tipVertex()] - inputVertexPositions[he.tailVertex()];
}

} // namespace surface
} // namespace geometrycentral
