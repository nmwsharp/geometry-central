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
    vertexPositions[v] = Vector3{x, y, z};
    iV++;
  }
}

inline double VertexPositionGeometry::edgeLength(Edge e) const {
  Halfedge he = e.halfedge();
  Vector3 pA = vertexPositions[he.vertex()];
  Vector3 pB = vertexPositions[he.next().vertex()];
  return norm(pA - pB);
}

// Face areas
inline double VertexPositionGeometry::faceArea(Face f) const {
  // WARNING: Logic duplicated between cached and immediate version
  Halfedge he = f.halfedge();
  Vector3 pA = vertexPositions[he.vertex()];
  he = he.next();
  Vector3 pB = vertexPositions[he.vertex()];
  he = he.next();
  Vector3 pC = vertexPositions[he.vertex()];

  GC_SAFETY_ASSERT(he.next() == f.halfedge(), "faces must be triangular");

  double area = 0.5 * norm(cross(pB - pA, pC - pA));
  return area;
}

// Vertex dual areas
inline double VertexPositionGeometry::vertexDualArea(Vertex v) const {
  // WARNING: Logic duplicated between cached and immediate version
  double area = 0.;
  for (Face f : v.adjacentFaces()) {
    area += faceArea(f);
  }
  return area / 3.;
}

// Corner angles
inline double VertexPositionGeometry::cornerAngle(Corner c) const {
  // WARNING: Logic duplicated between cached and immediate version
  Halfedge he = c.halfedge();
  Vector3 pA = vertexPositions[he.vertex()];
  he = he.next();
  Vector3 pB = vertexPositions[he.vertex()];
  he = he.next();
  Vector3 pC = vertexPositions[he.vertex()];

  GC_SAFETY_ASSERT(he.next() == c.halfedge(), "faces must be triangular");

  double q = dot(unit(pB - pA), unit(pC - pA));
  q = clamp(q, -1.0, 1.0);
  double angle = std::acos(q);
  return angle;
}

inline double VertexPositionGeometry::halfedgeCotanWeight(Halfedge heI) const {
  // WARNING: Logic duplicated between cached and immediate version
  if (heI.isInterior()) {
    Halfedge he = heI;
    Vector3 pB = vertexPositions[he.vertex()];
    he = he.next();
    Vector3 pC = vertexPositions[he.vertex()];
    he = he.next();
    Vector3 pA = vertexPositions[he.vertex()];
    GC_SAFETY_ASSERT(he.next() == heI, "faces must be triangular");

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
    Vector3 pA = vertexPositions[he.vertex()];
    he = he.next();
    Vector3 pB = vertexPositions[he.vertex()];
    he = he.next();
    Vector3 pC = vertexPositions[he.vertex()];

    normalSum += cross(pB - pA, pC - pA);

    // In the special case of a triangle, there is no need to to repeat at all three corners; the result will be the
    // same
    if (he.next() == heF) break;
  }

  Vector3 normal = unit(normalSum);
  return normal;
}


inline Vector3 VertexPositionGeometry::halfedgeVector(Halfedge he) const {
  return vertexPositions[he.tipVertex()] - vertexPositions[he.tailVertex()];
}

inline double VertexPositionGeometry::edgeDihedralAngle(Edge e) const {
  // WARNING: Logic duplicated between cached and immediate version
  if (e.isBoundary() || !e.isManifold()) {
    return 0.;
  }

  Vector3 N1 = faceNormal(e.halfedge().face());
  Vector3 N2 = faceNormal(e.halfedge().sibling().face());
  Vector3 pTail = vertexPositions[e.halfedge().vertex()];
  Vector3 pTip = vertexPositions[e.halfedge().next().vertex()];
  Vector3 edgeDir = unit(pTip - pTail);

  return atan2(dot(edgeDir, cross(N1, N2)), dot(N1, N2));
}

inline double VertexPositionGeometry::vertexMeanCurvature(Vertex v) const {
  // WARNING: Logic duplicated between cached and immediate version
  double meanCurvature = 0.;
  for (Halfedge he : v.outgoingHalfedges()) {
    double len = edgeLength(he.edge());
    double alpha = edgeDihedralAngle(he.edge());
    meanCurvature += alpha * len / 2.;
  }
  return meanCurvature / 2.;
}

inline double VertexPositionGeometry::vertexGaussianCurvature(Vertex v) const {
  // WARNING: Logic duplicated between cached and immediate version

  // the triangles neighboring any boundary vertex can be flattened into
  // the plane without any stretching/distortion; hence, a boundary
  // vertex has no Gaussian curvature
  if (v.isBoundary()) return 0.;

  double gaussianCurvature = 2. * PI;
  for (Corner c : v.adjacentCorners()) {
    gaussianCurvature -= cornerAngle(c);
  }
  return gaussianCurvature;
}

inline double VertexPositionGeometry::vertexMaxPrincipalCurvature(Vertex v) const {
  // WARNING: Logic duplicated between cached and immediate version
  return vertexPrincipalCurvature(2, v);
}

inline double VertexPositionGeometry::vertexMinPrincipalCurvature(Vertex v) const {
  // WARNING: Logic duplicated between cached and immediate version
  return vertexPrincipalCurvature(1, v);
}

inline double VertexPositionGeometry::vertexPrincipalCurvature(int whichCurvature, Vertex v) const {
  // WARNING: Logic duplicated between cached and immediate version
  double A = vertexDualArea(v);
  double H = vertexMeanCurvature(v) / A;
  double K = vertexGaussianCurvature(v) / A;

  // The two principal curvatures are given by
  //    H +/- sqrt( H^2 - K )
  double c = std::sqrt(std::max(0., H * H - K));
  double k1 = H - c;
  double k2 = H + c;

  if (whichCurvature == 1)
    return std::min(k1, k2);
  else
    return std::max(k1, k2);
}

inline Vector3 VertexPositionGeometry::vertexDualMeanCurvatureNormal(Vertex v) const {
  // WARNING: Logic duplicated between cached and immediate version
  Vector3 hN = Vector3::zero();
  for (Halfedge he : v.outgoingHalfedges()) {
    double w = edgeCotanWeight(he.edge());
    hN += w * (vertexPositions[v] - vertexPositions[he.tipVertex()]) / 2.;
  }
  return hN;
}

} // namespace surface
} // namespace geometrycentral
