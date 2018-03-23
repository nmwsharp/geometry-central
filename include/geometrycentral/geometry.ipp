#include <algorithm>
#include <cmath>

namespace geometrycentral {

// Generic methods =============================================================

template <class T>
inline HalfedgeMesh* Geometry<T>::getMesh(void) {
  return &mesh;
}

template <class T>
Geometry<T>* Geometry<T>::copyUsingTransfer(HalfedgeMeshDataTransfer& transfer) {

  VertexData<T> newCoords = transfer.transfer(p);
  Geometry<T>* newGeom = new Geometry<T>(*transfer.newMesh);

  for (VertexPtr v : transfer.newMesh->vertices()) {
    newGeom->position(v) = newCoords[v];
  }

  return newGeom;
}

// Vertex attributes (primal)

template <class T>
inline T& Geometry<T>::position(VertexPtr v) {
  return p[v];
}

template <class T>
inline T Geometry<T>::position(VertexPtr v) const {
  return p[v];
}

template <class T>
inline double Geometry<T>::volume(VertexPtr v) {
  return 1.;
}

template <class T>
inline double Geometry<T>::dualArea(VertexPtr v) {
  double sum = 0;
  for (FacePtr f : v.adjacentFaces()) {
    sum += area(f);
  }
  return sum / 3.0;
}

template <class T>
inline double Geometry<T>::angleDefect(VertexPtr v) {
  double sum = 0.;

  for (HalfedgePtr h : v.outgoingHalfedges()) {
    if (h.isReal()) sum += angle(h.next());
  }

  return 2. * M_PI - sum;
}

template <class T>
inline Vector3 Geometry<T>::normal(VertexPtr v) {
  Vector3 N{0., 0., 0.};

  for (FacePtr f : v.adjacentFaces()) {
    N += area(f) * normal(f);
  }

  return unit(N);
}

template <class T>
inline Vector3 Geometry<T>::boundaryNormal(VertexPtr v) {
  if (!v.isBoundary()) {
    return Vector3{0., 0., 0.};
  }

  Vector3 N1{0., 0., 0.};
  Vector3 N2{0., 0., 0.};
  for (HalfedgePtr he : v.incomingHalfedges()) {
    if (!he.isReal()) {
      N1 = vector(he).rotate_around(normal(v), PI / 2.0);
    }
    if (!he.twin().isReal()) {
      N2 = vector(he.twin()).rotate_around(normal(v), PI / 2.0);
    }
  }

  return unit(N1 + N2);
}

template <class T>
inline Complex Geometry<T>::tangentVectorToComplexAngle(VertexPtr v, const Vector3& inVec) {
  // Assumes vector is tangent
  Vector3 N = normal(v);
  Vector3 refEdge = unit(projectToTangentSpace(v, vector(v.halfedge())));
  Vector3 V = unit(inVec);

  double realPart = dot(refEdge, V);
  double imagPart = dot(refEdge.rotate_around(N, PI / 2.0), V);

  return Complex(realPart, imagPart);
}

template <class T>
inline Vector3 Geometry<T>::complexAngleToTangentVector(VertexPtr v, Complex inAngle) {
  // Assumes vector is tangent
  Vector3 N = normal(v);
  Vector3 refEdge = unit(projectToTangentSpace(v, vector(v->halfedge)));
  double theta = std::arg(inAngle);

  return refEdge.rotate_around(N, theta) * std::abs(inAngle);
}

template <class T>
inline Vector3 Geometry<T>::projectToTangentSpace(VertexPtr v, const Vector3& inVec) {
  Vector3 N = normal(v);
  return inVec - dot(inVec, N) * N;
}

template <class T>
Complex Geometry<T>::principalDirection(VertexPtr v) {
  // NOTE: This logic is duplicated here and in the cached method

  Complex principalDir(0, 0);

  for (HalfedgePtr he : v.outgoingHalfedges()) {
    double len = length(he->edge);
    double alpha = dihedralAngle(he->edge);
    double theta = angularCoordinate(he);

    Complex r2(cos(2 * theta), sin(2 * theta));

    principalDir += -len * std::abs(alpha) * r2;
  }

  return principalDir / 4.0;
}


// Edge attributes (primal)

template <class T>
inline T Geometry<T>::midpoint(EdgePtr e) {
  return .5 * (p[e.halfedge().vertex()] + p[e.halfedge().twin().vertex()]);
}

template <class T>
inline double Geometry<T>::length(EdgePtr e) {
  HalfedgePtr h = e.halfedge();

  T p0 = p[h.vertex()];
  h = h.twin();
  T p1 = p[h.vertex()];

  return norm(p1 - p0);
}

template <class T>
inline double Geometry<T>::cotanWeight(EdgePtr e) {
  HalfedgePtr h = e.halfedge();
  return 0.5 * (cotan(h) + cotan(h.twin()));
}

template <class T>
inline double Geometry<T>::dihedralAngle(EdgePtr e) {
  if (e.isBoundary()) {
    return 0;
  }

  Vector3 N1 = normal(e.halfedge().face());
  Vector3 N2 = normal(e.halfedge().twin().face());
  Vector3 edgeV = unit(-vector(e.halfedge()));

  return atan2(dot(edgeV, cross(N1, N2)), dot(N1, N2));
}


// Face attributes (primal)

template <class T>
inline double Geometry<T>::area(FacePtr f) {
  if (!f.isReal()) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  return norm(areaVector(f));
}

template <class T>
inline Vector3 Geometry<T>::normal(FacePtr f) {
  if (!f.isReal()) {
    return Vector3::undefined();
  }

  return unit(areaVector(f));
}

template <class T>
inline Vector3 Geometry<T>::areaVector(FacePtr f) {
  if (!f.isReal()) {
    return Vector3::undefined();
  }

  Vector3 AN{0., 0., 0.};
  for (HalfedgePtr h : f.adjacentHalfedges()) {
    Vector3 pi = position(h.vertex());
    Vector3 pj = position(h.twin().vertex());
    AN += cross(pi, pj);
  }

  return AN / 2.;
}

template <class T>
inline T Geometry<T>::barycenter(FacePtr f) {
  if (!f.isReal()) {
    return Vector3::undefined();
  }

  T sum = T::zero();
  double k = 0.;
  for (VertexPtr v : f.adjacentVertices()) {
    sum += p[v];
    k += 1.;
  }

  return sum / k;
}

template <class T>
inline T Geometry<T>::circumcenter(FacePtr f) {
  if (!f.isReal()) {
    return Vector3::undefined();
  }

  HalfedgePtr h = f.halfedge();
  const Vector3& a = position(h.vertex());
  h = h.next();
  const Vector3& b = position(h.vertex());
  h = h.next();
  const Vector3& c = position(h.vertex());

  if (h.next() != f.halfedge()) {
    throw std::domain_error("Circumcenter only defined for triangles");
  }

  return a + (norm2(c - a) * cross(cross(b - a, c - a), b - a) + norm2(b - a) * cross(c - a, cross(b - a, c - a))) /
                 (2. * norm2(cross(b - a, c - a)));
}


// Halfedge attributes (primal)

template <class T>
inline double Geometry<T>::angle(HalfedgePtr h) {
  T t1 = vector(h.next().next());
  T t2 = -vector(h.next());

  return geometrycentral::angle(t1, t2);
}

template <class T>
inline double Geometry<T>::angle(CornerPtr c) {
  return angle(c.halfedge());
}

template <class T>
inline T Geometry<T>::vector(HalfedgePtr h) {
  VertexPtr v0 = h.vertex();
  VertexPtr v1 = h.twin().vertex();

  return p[v1] - p[v0];
}

template <class T>
inline double Geometry<T>::cotan(HalfedgePtr h) {
  if (!h.isReal()) {
    return 0.;
  }

  T v1 = vector(h.next().next());
  T v2 = -vector(h.next());

  return dot(v1, v2) / norm(cross(v1, v2));
}

template <class T>
inline double Geometry<T>::angularCoordinate(HalfedgePtr h) {
  double coord = 0;
  double angleSum = 0;

  for (HalfedgePtr he : h.vertex().outgoingHalfedges()) {
    if (he == h) {
      coord = angleSum;
    }
    angleSum += angle(he.twin().next().next());
  }

  coord *= 2 * PI / angleSum; // normalize to 2PI
  return 2 * PI - coord;      // we want the CCW value, but we just measured the CW
                              // value (due to the direction orbited by twin->next)
}

// Spherical specializations ===================================================

template <>
inline double Geometry<Spherical>::length(EdgePtr e) {
  UnitVector3 u0 = p[e.halfedge().vertex()];
  UnitVector3 u1 = p[e.halfedge().twin().vertex()];

  return geometrycentral::angle(u0, u1);
}

template <>
inline double Geometry<Spherical>::cotanWeight(EdgePtr e) {
  throw std::domain_error("Edge cotangent weights not meaningful/useful for spherical geometry.");
}

template <>
inline double Geometry<Spherical>::area(FacePtr f) {
  if (!f.isReal()) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  double sum = 0.;

  for (HalfedgePtr h : f.adjacentHalfedges()) {
    sum += angle(h);
  }

  return sum - M_PI;
}

template <>
inline Vector3 Geometry<Spherical>::normal(FacePtr f) {
  HalfedgePtr h = f.halfedge();

  Vector3 pi = p[h.vertex()];
  h = h.next();
  Vector3 pj = p[h.vertex()];
  h = h.next();
  Vector3 pk = p[h.vertex()];

  double thetaIJ = geometrycentral::angle(pi, pj);
  double thetaJK = geometrycentral::angle(pj, pk);
  double thetaKI = geometrycentral::angle(pk, pi);

  return thetaIJ * unit(cross(pj, pi)) + thetaJK * unit(cross(pk, pj)) + thetaKI * unit(cross(pi, pk));
}

template <>
inline UnitVector3 Geometry<Spherical>::barycenter(FacePtr f) {
  const double eps = 1e-7;
  const double delta = 1e-5;
  Vector3 u{0., 0., 0.}; // mean tangent vector
  double phi;

  // Use normalized Euclidean average as initial guess for barycenter
  Vector3 c{0., 0., 0.};
  for (HalfedgePtr h : f.adjacentHalfedges()) {
    c += p[h.vertex()];
  }

  // Iteratively average in Lie algebra until convergence
  do {
    for (HalfedgePtr h : f.adjacentHalfedges()) {
      Vector3 q = p[h.vertex()];

      double theta = geometrycentral::angle(c, q);
      Vector3 w = cross(c, q);
      double m = norm(w);
      if (m > eps) {
        w /= m;
        Vector3 v = theta * unit(cross(w, c));
        u += v;
      }
    }

    u /= 3.;

    phi = norm(u);

    if (phi > eps) {
      c = cos(phi) * c + sin(phi) * unit(u);
    }
  } while (phi > delta);

  return c;
}

template <>
inline UnitVector3 Geometry<Spherical>::vector(HalfedgePtr h) {
  // Since the sphere is not a vector space, we cannot just return
  // the difference of points on the sphere.  Instead, we do the next
  // best thing and return a unit vector tangent to the first vertex and
  // pointing toward the second vertex.  Note that if this vector is
  // subsequently scaled by the angle between the two points, then the
  // exponential map at the first point will map this vector to the
  // second point.
  UnitVector3 u0 = p[h.vertex()];
  UnitVector3 u1 = p[h.twin().vertex()];

  return UnitVector3(cross(cross(u0, u1), u0));
}

template <>
inline double Geometry<Spherical>::cotan(HalfedgePtr h) {
  throw std::domain_error("Halfedge cotangents not meaningful/useful for spherical geometry.");
}

template <>
inline UnitVector3 Geometry<Spherical>::center(void) {
  throw std::domain_error("Center of mass not implemented for spherical geometry.");
}

// Global attributes ===========================================================

template <class T>
inline double Geometry<T>::totalArea(void) {
  double A = 0.;

  for (FacePtr f : mesh.faces()) {
    A += area(f);
  }

  return A;
}

template <class T>
inline T Geometry<T>::center(void) {
  double A = 0.;
  T c = T::zero();

  for (FacePtr f : mesh.faces()) {
    double Af = area(f);
    A += Af;
    c += Af * barycenter(f);
  }

  return c / A;
}

template <class T>
inline void Geometry<T>::boundingBox(T& bboxMin, T& bboxMax) {
  bboxMin = T::infinity();
  bboxMax = -T::infinity();

  for (VertexPtr v : mesh.vertices()) {
    bboxMin = componentwiseMin(bboxMin, p[v]);
    bboxMax = componentwiseMax(bboxMax, p[v]);
  }
}

template <class T>
inline T Geometry<T>::extent(void) {
  T bboxMin, bboxMax;
  boundingBox(bboxMin, bboxMax);
  return bboxMax - bboxMin;
}

template <class T>
inline double Geometry<T>::lengthScale(void) {
  return norm(extent());
}

// Convenience methods for caching current attributes ==========================

template <class T>
void Geometry<T>::getVertexPositions(VertexData<T>& vertexPosition) {
  vertexPosition = p;
}

template <class T>
void Geometry<T>::getVertexNormals(VertexData<Vector3>& vertexNormal) {
  vertexNormal = VertexData<Vector3>(&mesh);
  for (VertexPtr v : mesh.vertices()) {
    vertexNormal[v] = normal(v);
  }
}

template <class T>
void Geometry<T>::getVertexAngleDefects(VertexData<double>& vertexAngleDefect) {
  vertexAngleDefect = VertexData<double>(&mesh);
  for (VertexPtr v : mesh.vertices()) {
    vertexAngleDefect[v] = angleDefect(v);
  }
}

template <class T>
void Geometry<T>::getPrincipalDirections(VertexData<Complex>& principalDirections) {
  HalfedgeData<double> angularCoordinates;
  getAngularCoordinates(angularCoordinates);
  getPrincipalDirections(principalDirections, angularCoordinates);
}

template <class T>
void Geometry<T>::getPrincipalDirections(VertexData<Complex>& principalDirections,
                                         HalfedgeData<double>& angularCoordinates) {
  principalDirections = VertexData<Complex>(&mesh);

  // NOTE: This logic is duplicated here and in the per-vertex method (to use a
  // cached coordinate array)

  for (VertexPtr v : mesh.vertices()) {
    Complex principalDir(0, 0);

    for (HalfedgePtr he : v.outgoingHalfedges()) {
      double len = length(he.edge());
      double alpha = dihedralAngle(he.edge());
      double theta = angularCoordinates[he];

      Complex r2(cos(2 * theta), sin(2 * theta));

      principalDir += -len * std::abs(alpha) * r2;
    }

    principalDirections[v] = principalDir / 4.0;
  }
}


template <class T>
void Geometry<T>::getEdgeLengths(EdgeData<double>& edgeLength) {
  edgeLength = EdgeData<double>(&mesh);
  for (EdgePtr e : mesh.edges()) {
    edgeLength[e] = length(e);
  }
}

template <class T>
void Geometry<T>::getEdgeCotanWeights(EdgeData<double>& edgeCotanWeight) {
  edgeCotanWeight = EdgeData<double>(&mesh);
  for (EdgePtr e : mesh.edges()) {
    edgeCotanWeight[e] = cotanWeight(e);
  }
}

template <class T>
void Geometry<T>::getFaceAreas(FaceData<double>& faceArea) {
  faceArea = FaceData<double>(&mesh);
  for (FacePtr f : mesh.faces()) {
    faceArea[f] = area(f);
  }
}


template <class T>
void Geometry<T>::getFaceNormals(FaceData<Vector3>& faceNormal) {
  faceNormal = FaceData<Vector3>(&mesh);
  for (FacePtr f : mesh.faces()) {
    faceNormal[f] = normal(f);
  }
}

template <class T>
void Geometry<T>::getFaceBarycenters(FaceData<T>& faceBarycenter) {
  faceBarycenter = FaceData<T>(&mesh);
  for (FacePtr f : mesh.faces()) {
    faceBarycenter[f] = barycenter(f);
  }
}

template <class T>
void Geometry<T>::getHalfedgeVectors(HalfedgeData<T>& halfedgeVector) {
  halfedgeVector = HalfedgeData<T>(&mesh);
  for (HalfedgePtr h : mesh.allHalfedges()) {
    halfedgeVector[h] = vector(h);
  }
}

template <class T>
void Geometry<T>::getHalfedgeAngles(HalfedgeData<double>& halfedgeAngle) {
  halfedgeAngle = HalfedgeData<double>(&mesh);
  for (HalfedgePtr h : mesh.allHalfedges()) {
    halfedgeAngle[h] = angle(h);
  }
}

template <class T>
void Geometry<T>::getCornerAngles(CornerData<double>& cornerAngle) {
  cornerAngle = CornerData<double>(&mesh);
  for (CornerPtr c : mesh.corners()) {
    cornerAngle[c] = angle(c);
  }
}

template <class T>
void Geometry<T>::getHalfedgeCotans(HalfedgeData<double>& halfedgeCotan) {
  halfedgeCotan = HalfedgeData<double>(&mesh);
  for (HalfedgePtr h : mesh.realHalfedges()) {
    halfedgeCotan[h] = cotan(h);
  }
}

template <class T>
void Geometry<T>::getAngularCoordinates(HalfedgeData<double>& angularCoordinates) {
  angularCoordinates = HalfedgeData<double>(&mesh);
  for (HalfedgePtr h : mesh.allHalfedges()) {
    angularCoordinates[h] = angularCoordinate(h);
  }
}


} // namespace geometrycentral
