#pragma once

namespace geometrycentral {
namespace surface {

// == Constructors
inline BarycentricVector::BarycentricVector() : type(BarycentricVectorType::Face) {}
inline BarycentricVector::BarycentricVector(Face f)
    : type(BarycentricVectorType::Face), face(f), faceCoords(Vector3::zero()) {}
inline BarycentricVector::BarycentricVector(Vertex v) : type(BarycentricVectorType::Vertex), vertex(v) {}
inline BarycentricVector::BarycentricVector(Edge e, Vector2 edgeCoords_)
    : type(BarycentricVectorType::Edge), edge(e), edgeCoords(edgeCoords_) {}
inline BarycentricVector::BarycentricVector(Face f, Vector3 faceCoords_)
    : type(BarycentricVectorType::Face), face(f), faceCoords(faceCoords_) {}
inline BarycentricVector::BarycentricVector(SurfacePoint pA, SurfacePoint pB) {

  // Detect if pA, pB share face, edge, or vertex.
  if (pA.type == SurfacePointType::Vertex && pB.type == SurfacePointType::Vertex) {
    if (pA.vertex == pB.vertex) {
      type = BarycentricVectorType::Vertex;
      vertex = pA.vertex;
      return;
    }
  }

  Edge e = sharedEdge(pA, pB);
  if (e != Edge()) {
    type = BarycentricVectorType::Edge;
    SurfacePoint pA_in_edge = pA.inEdge(e);
    SurfacePoint pB_in_edge = pB.inEdge(e);
    double a = pA_in_edge.tEdge;
    double b = pB_in_edge.tEdge;
    edge = e;
    edgeCoords = Vector2{b - a, a - b};
    return;
  }

  Face f = sharedFace(pA, pB);
  if (f == Face()) throw std::invalid_argument("Input SurfacePoints of a BarycentricVector must share a face.");
  type = BarycentricVectorType::Face;
  SurfacePoint pA_in_face = pA.inFace(f);
  SurfacePoint pB_in_face = pB.inFace(f);
  face = f;
  faceCoords = pB_in_face.faceCoords - pA_in_face.faceCoords;
}

// == Methods

inline BarycentricVector BarycentricVector::inSomeFace() const {

  switch (type) {
  case BarycentricVectorType::Face:
    return *this;
    break;
  case BarycentricVectorType::Edge: {
    Halfedge targetHe = edge.halfedge();
    Face inFace = targetHe.face();
    Halfedge he = inFace.halfedge(); // tail should be the "first" vertex around this face
    if (targetHe == he) {
      return BarycentricVector(inFace, Vector3{edgeCoords[0], edgeCoords[1], 0.});
    }
    he = he.next();
    if (targetHe == he) {
      return BarycentricVector(inFace, Vector3{0., edgeCoords[0], edgeCoords[1]});
    }
    return BarycentricVector(inFace, Vector3{edgeCoords[1], 0., edgeCoords[0]});
    break;
  }
  case BarycentricVectorType::Vertex: {
    Halfedge he = vertex.halfedge();
    Face inFace = he.face();
    return BarycentricVector(inFace, Vector3::zero());
    break;
  }
  }

  throw std::logic_error("bad switch");
  return *this;
}

inline BarycentricVector BarycentricVector::inFace(Face targetFace) const {

  switch (type) {
  case BarycentricVectorType::Face:
    if (face == targetFace) {
      return *this;
    };
    break;
  case BarycentricVectorType::Edge: {
    Halfedge he = targetFace.halfedge();
    if (edge == he.edge()) {
      Vector3 faceCoords = (he == edge.halfedge()) ? Vector3{edgeCoords[0], edgeCoords[1], 0.}
                                                   : Vector3{edgeCoords[1], edgeCoords[0], 0.};
      return BarycentricVector(targetFace, faceCoords);
    }
    he = he.next();
    if (edge == he.edge()) {
      Vector3 faceCoords = (he == edge.halfedge()) ? Vector3{0., edgeCoords[0], edgeCoords[1]}
                                                   : Vector3{0., edgeCoords[1], edgeCoords[0]};
      return BarycentricVector(targetFace, faceCoords);
    }
    he = he.next();
    Vector3 faceCoords =
        (he == edge.halfedge()) ? Vector3{edgeCoords[1], 0., edgeCoords[0]} : Vector3{edgeCoords[0], 0., edgeCoords[1]};
    return BarycentricVector(targetFace, faceCoords);
    break;
  }
  case BarycentricVectorType::Vertex:
    for (Vertex v : targetFace.adjacentVertices()) {
      if (v == vertex) return BarycentricVector(targetFace, Vector3::zero());
    }
    break;
  }

  throw std::logic_error("BarycentricVector " + std::to_string(*this) + " not adjacent to target face " +
                         std::to_string(targetFace));
  return *this;
}

inline BarycentricVector BarycentricVector::inEdge(Edge targetEdge) const {

  switch (type) {
  case BarycentricVectorType::Face:
    break;
  case BarycentricVectorType::Edge:
    if (edge == targetEdge) {
      return *this;
    }
    break;
  case BarycentricVectorType::Vertex:
    for (Vertex v : targetEdge.adjacentVertices()) {
      if (v == vertex) return BarycentricVector(targetEdge, Vector2::zero());
    }
    break;
  }

  throw std::logic_error("BarycentricVector " + std::to_string(*this) + " not adjacent to target edge " +
                         std::to_string(targetEdge));
  return *this;
}

inline void BarycentricVector::validate() const {

  double eps = 1e-5;

  switch (type) {
  case BarycentricVectorType::Face:
    if (face == Face()) throw std::logic_error("BarycentricVector with Type::Face has invalid face ref");
    if (!isfinite(faceCoords)) throw std::logic_error("BarycentricVector with Type::Face has non-finite coords");
    if (abs(sum(faceCoords)) > eps)
      throw std::logic_error("BarycentricVector with Type::Face has coords that do not sum to 0 " +
                             std::to_string(faceCoords));
    break;
  case BarycentricVectorType::Edge:
    if (edge == Edge()) throw std::logic_error("BarycentricVector with Type::Edge has invalid edge ref");
    if (!isfinite(edgeCoords)) throw std::logic_error("BarycentricVector with Type::Edge has non-finite coords");
    if (abs(sum(edgeCoords)) > eps)
      throw std::logic_error("BarycentricVector with Type::Edge has coords that do not sum to 0 " +
                             std::to_string(edgeCoords));
    break;
  case BarycentricVectorType::Vertex:
    if (vertex == Vertex()) throw std::logic_error("BarycentricVector with Type::Vertex has invalid vertex ref");
    break;
  }
}

inline void BarycentricVector::normalizeDisplacement() {
  switch (type) {
  case BarycentricVectorType::Face:
    faceCoords -= Vector3::constant(sum(faceCoords) / 3.);
    break;
  case BarycentricVectorType::Edge:
    edgeCoords -= Vector2::constant(sum(edgeCoords) / 2.);
    break;
  default:
    break;
  }
}

inline BarycentricVector BarycentricVector::normalizedDisplacement() const {
  switch (type) {
  case BarycentricVectorType::Face:
    return BarycentricVector(face, faceCoords - Vector3::constant(sum(faceCoords) / 3.));
    break;
  case BarycentricVectorType::Edge:
    return BarycentricVector(edge, edgeCoords - Vector2::constant(sum(edgeCoords) / 2.));
    break;
  default:
    break;
  }
  return *this;
}

inline BarycentricVector BarycentricVector::rotated90(IntrinsicGeometryInterface& geom) const {

  switch (type) {
  case BarycentricVectorType::Face: {
    double ui = faceCoords[0];
    double uj = faceCoords[1];
    double uk = faceCoords[2];
    geom.requireEdgeLengths();
    double l_ij = geom.edgeLengths[face.halfedge().edge()];
    double l_jk = geom.edgeLengths[face.halfedge().next().edge()];
    double l_ki = geom.edgeLengths[face.halfedge().next().next().edge()];
    geom.unrequireEdgeLengths();
    // coefficients of matrix D taking barycentric coords to 3D local coords, squared and multiplied by 2
    double ai = l_ij * l_ij - l_jk * l_jk + l_ki * l_ki;
    double aj = l_jk * l_jk - l_ki * l_ki + l_ij * l_ij;
    double ak = l_ki * l_ki - l_ij * l_ij + l_jk * l_jk;
    double s = 0.5 * (l_ij + l_jk + l_ki);
    double A = std::sqrt(s * (s - l_ij) * (s - l_jk) * (s - l_ki));
    Vector3 newCoords = {ak * uk - aj * uj, ai * ui - ak * uk, aj * uj - ai * ui};
    newCoords /= (4. * A);
    return BarycentricVector(face, newCoords);
    break;
  }
  case BarycentricVectorType::Edge:
    throw std::logic_error("Cannot rotate BarycentricVector of Type::Edge");
    break;
  default:
    // If vertex-type, do nothing
    break;
  }
  return *this;
}

// == Overloaded operators

inline BarycentricVector BarycentricVector::operator+(const BarycentricVector& w) const {

  // If the two vectors share a vertex, just return either operand.
  if (type == BarycentricVectorType::Vertex && w.type == BarycentricVectorType::Vertex && vertex == w.vertex) {
    return *this;
  }
  // If the two vectors share an edge, return a vector also on this edge.
  Edge e = sharedEdge(*this, w);
  if (e != Edge()) {
    BarycentricVector uEdge = this->inEdge(e);
    BarycentricVector wEdge = w.inEdge(e);
    return BarycentricVector(e, uEdge.edgeCoords + wEdge.edgeCoords);
  }
  // If the two vectors share a face, return a vector in this face.
  Face f = sharedFace(*this, w);
  if (f != Face()) {
    BarycentricVector uFace = this->inFace(f);
    BarycentricVector wFace = w.inFace(f);
    return BarycentricVector(f, uFace.faceCoords + wFace.faceCoords);
  }

  throw std::logic_error("BarycentricVectors must share a face to be added");
  return *this;
}

inline BarycentricVector BarycentricVector::operator-(const BarycentricVector& w) const {

  // If the two vectors share a vertex, just return either operand.
  if (type == BarycentricVectorType::Vertex && w.type == BarycentricVectorType::Vertex && vertex == w.vertex) {
    return *this;
  }
  // If the two vectors share an edge, return a vector also on this edge.
  Edge e = sharedEdge(*this, w);
  if (e != Edge()) {
    BarycentricVector uEdge = this->inEdge(e);
    BarycentricVector wEdge = w.inEdge(e);
    return BarycentricVector(e, uEdge.edgeCoords - wEdge.edgeCoords);
  }
  // If the two vectors share a face, return a vector in this face.
  Face f = sharedFace(*this, w);
  if (f != Face()) {
    BarycentricVector uFace = this->inFace(f);
    BarycentricVector wFace = w.inFace(f);
    return BarycentricVector(f, uFace.faceCoords - wFace.faceCoords);
  }

  throw std::logic_error("BarycentricVectors must share a face to be subtracted");
  return *this;
}

inline BarycentricVector BarycentricVector::operator*(double s) const {
  switch (type) {
  case BarycentricVectorType::Face:
    return BarycentricVector(face, faceCoords * s);
    break;
  case BarycentricVectorType::Edge:
    return BarycentricVector(edge, edgeCoords * s);
    break;
  default:
    return *this;
    break;
  }
}

inline BarycentricVector BarycentricVector::operator/(double s) const {
  switch (type) {
  case BarycentricVectorType::Face:
    return BarycentricVector(face, faceCoords / s);
    break;
  case BarycentricVectorType::Edge:
    return BarycentricVector(edge, edgeCoords / s);
    break;
  default:
    break;
  }
  return *this;
}

inline const BarycentricVector BarycentricVector::operator-() const {
  switch (type) {
  case BarycentricVectorType::Face:
    return BarycentricVector(face, -faceCoords);
    break;
  case BarycentricVectorType::Edge:
    return BarycentricVector(edge, -edgeCoords);
    break;
  default:
    break;
  }
  return *this;
}

template <typename T>
inline BarycentricVector operator*(const T s, const BarycentricVector& v) {
  switch (v.type) {
  case BarycentricVectorType::Face:
    return BarycentricVector(v.face, s * v.faceCoords);
    break;
  case BarycentricVectorType::Edge:
    return BarycentricVector(v.edge, s * v.edgeCoords);
    break;
  default:
    break;
  }
  return v;
}

inline BarycentricVector& BarycentricVector::operator+=(const BarycentricVector& w) {

  // If the two vectors share a vertex, do nothing.
  if (type == BarycentricVectorType::Vertex && w.type == BarycentricVectorType::Vertex && vertex == w.vertex) {
    return *this;
  }
  // If the two vectors share an edge, change *this to be a vector also on this edge.
  Edge e = sharedEdge(*this, w);
  if (e != Edge()) {
    BarycentricVector uEdge = this->inEdge(e);
    BarycentricVector wEdge = w.inEdge(e);
    *this = BarycentricVector(e, uEdge.edgeCoords + wEdge.edgeCoords);
    return *this;
  }
  // If the two vectors share a face, change *this to be a vector in this face.
  Face f = sharedFace(*this, w);
  if (f != Face()) {
    BarycentricVector uFace = this->inFace(f);
    BarycentricVector wFace = w.inFace(f);
    *this = BarycentricVector(f, uFace.faceCoords + wFace.faceCoords);
    return *this;
  }

  throw std::logic_error("BarycentricVectors must share a face to be added");
  return *this;
}

inline BarycentricVector& BarycentricVector::operator-=(const BarycentricVector& w) {

  // If the two vectors share a vertex, do nothing.
  if (type == BarycentricVectorType::Vertex && w.type == BarycentricVectorType::Vertex && vertex == w.vertex) {
    return *this;
  }
  // If the two vectors share an edge, change *this to be a vector also on this edge.
  Edge e = sharedEdge(*this, w);
  if (e != Edge()) {
    BarycentricVector uEdge = this->inEdge(e);
    BarycentricVector wEdge = w.inEdge(e);
    *this = BarycentricVector(e, uEdge.edgeCoords - wEdge.edgeCoords);
    return *this;
  }
  // If the two vectors share a face, change *this to be a vector in this face.
  Face f = sharedFace(*this, w);
  if (f != Face()) {
    BarycentricVector uFace = this->inFace(f);
    BarycentricVector wFace = w.inFace(f);
    *this = BarycentricVector(f, uFace.faceCoords - wFace.faceCoords);
    return *this;
  }

  throw std::logic_error("BarycentricVectors must share a face to be subtracted");
  return *this;
}

inline BarycentricVector& BarycentricVector::operator*=(const double& s) {
  switch (type) {
  case BarycentricVectorType::Face:
    faceCoords *= s;
    break;
  case BarycentricVectorType::Edge:
    edgeCoords *= s;
    break;
  default:
    break;
  }
  return *this;
}

inline BarycentricVector& BarycentricVector::operator/=(const double& s) {
  switch (type) {
  case BarycentricVectorType::Face:
    faceCoords /= s;
    break;
  case BarycentricVectorType::Edge:
    edgeCoords /= s;
    break;
  default:
    break;
  }
  return *this;
}

inline bool BarycentricVector::operator==(const BarycentricVector& other) const {

  if (type != other.type) return false;

  switch (type) {
  case BarycentricVectorType::Face:
    return (face == other.face) && (faceCoords == other.faceCoords);
    break;
  case BarycentricVectorType::Edge:
    return (edge == other.edge) && (edgeCoords == other.edgeCoords);
    break;
  case BarycentricVectorType::Vertex:
    return (vertex == other.vertex);
    break;
  }
  return false; // should never be reached
}

inline bool BarycentricVector::operator!=(const BarycentricVector& other) const { return !(*this == other); }

inline double BarycentricVector::norm(IntrinsicGeometryInterface& geom) const {

  geom.requireEdgeLengths();
  double val = 0.;

  switch (type) {
  case BarycentricVectorType::Face: {
    double ui = faceCoords[0];
    double uj = faceCoords[1];
    double uk = faceCoords[2];
    double l_ij = geom.edgeLengths[face.halfedge().edge()];
    double l_jk = geom.edgeLengths[face.halfedge().next().edge()];
    double l_ki = geom.edgeLengths[face.halfedge().next().next().edge()];
    val = std::sqrt(-((l_ij * l_ij * ui * uj) + (l_jk * l_jk * uj * uk) + (l_ki * l_ki * uk * ui)));
    break;
  }
  case BarycentricVectorType::Edge: {
    double l = geom.edgeLengths[edge];
    val = std::sqrt(-l * l * edgeCoords[0] * edgeCoords[1]);
    break;
  }
  case BarycentricVectorType::Vertex:
    break;
  }
  geom.unrequireEdgeLengths();
  return val;
}

inline double BarycentricVector::norm2(IntrinsicGeometryInterface& geom) const {

  geom.requireEdgeLengths();
  double val = 0.;

  switch (type) {
  case BarycentricVectorType::Face: { // Could do a loop with circular shifts, but I'd rather the code be easiser to
                                      // read.
    double ui = faceCoords[0];
    double uj = faceCoords[1];
    double uk = faceCoords[2];
    double l_ij = geom.edgeLengths[face.halfedge().edge()];
    double l_jk = geom.edgeLengths[face.halfedge().next().edge()];
    double l_ki = geom.edgeLengths[face.halfedge().next().next().edge()];
    val = -((l_ij * l_ij * ui * uj) + (l_jk * l_jk * uj * uk) + (l_ki * l_ki * uk * ui));
    break;
  }
  case BarycentricVectorType::Edge: {
    double l = geom.edgeLengths[edge];
    val = -l * l * edgeCoords[0] * edgeCoords[1];
    break;
  }
  case BarycentricVectorType::Vertex:
    break;
  }
  geom.unrequireEdgeLengths();
  return val;
}

inline Face sharedFace(const BarycentricVector& u, const BarycentricVector& w) {

  switch (u.type) {
  case BarycentricVectorType::Face: {
    switch (w.type) {
    case BarycentricVectorType::Face:
      if (u.face == w.face) return u.face;
      break;
    case BarycentricVectorType::Edge:
      for (Edge e : u.face.adjacentEdges()) {
        if (e == w.edge) return u.face;
      }
      break;
    case BarycentricVectorType::Vertex:
      for (Vertex v : u.face.adjacentVertices()) {
        if (v == w.vertex) return u.face;
        break;
      }
    }
    break;
  }
  case BarycentricVectorType::Edge: {
    switch (w.type) {
    case BarycentricVectorType::Face:
      for (Edge e : w.face.adjacentEdges()) {
        if (e == u.edge) return w.face;
        break;
      }
    case BarycentricVectorType::Edge:
      for (Face f : u.edge.adjacentFaces()) {
        for (Edge e : f.adjacentEdges()) {
          if (e == w.edge) return f;
        }
      }
      break;
    case BarycentricVectorType::Vertex:
      for (Face f : u.edge.adjacentFaces()) {
        for (Vertex v : f.adjacentVertices()) {
          if (v == w.vertex) return f;
        }
      }
      break;
    }
    break;
  }
  case BarycentricVectorType::Vertex: {
    switch (w.type) {
    case BarycentricVectorType::Face:
      for (Vertex v : w.face.adjacentVertices()) {
        if (v == u.vertex) return w.face;
      }
      break;
    case BarycentricVectorType::Edge:
      for (Face f : w.edge.adjacentFaces()) {
        for (Vertex v : f.adjacentVertices()) {
          if (v == u.vertex) return f;
        }
      }
      break;
    case BarycentricVectorType::Vertex:
      for (Face f : u.vertex.adjacentFaces()) {
        for (Vertex v : f.adjacentVertices()) {
          if (v == w.vertex) return f;
        }
      }
      break;
    }
    break;
  }
  }
  // No shared face
  return Face();
}

inline Edge sharedEdge(const BarycentricVector& u, const BarycentricVector& w) {

  if (u.type == BarycentricVectorType::Face || w.type == BarycentricVectorType::Face) return Edge();

  if (u.type == BarycentricVectorType::Edge) {
    Edge e = u.edge;
    if (w.type == BarycentricVectorType::Edge) {
      if (e == w.edge) return e;
    } else {
      for (Vertex v : e.adjacentVertices()) {
        if (v == w.vertex) return e;
      }
    }
  } else {
    if (w.type == BarycentricVectorType::Edge) {
      Edge e = w.edge;
      for (Vertex v : e.adjacentVertices()) {
        if (u.vertex == v) return e;
      }
    } else {
      for (Halfedge he : u.vertex.outgoingHalfedges()) {
        if (he.tipVertex() == w.vertex) return he.edge();
      }
    }
  }
  // No shared edge
  return Edge();
}

inline double norm(IntrinsicGeometryInterface& geom, const BarycentricVector& v) { return v.norm(geom); }
inline double norm2(IntrinsicGeometryInterface& geom, const BarycentricVector& v) { return v.norm2(geom); }

// Multiplication
// inline BarycentricVector cross(IntrinsicGeometryInterface& geom, const BarycentricVector& u, const
// BarycentricVector& v) {}

inline double dot(IntrinsicGeometryInterface& geom, const BarycentricVector& u, const BarycentricVector& v) {

  // First convert to shared face.
  Face f = sharedFace(u, v);
  if (f == Face())
    throw std::logic_error("Cannot compute inner product of BarycentricVectors that do not share a face.");

  BarycentricVector uf = u.inFace(f);
  BarycentricVector vf = v.inFace(f);

  double ui = uf.faceCoords[0];
  double uj = uf.faceCoords[1];
  double uk = uf.faceCoords[2];

  double vi = vf.faceCoords[0];
  double vj = vf.faceCoords[1];
  double vk = vf.faceCoords[2];

  geom.requireEdgeLengths();
  double l_ij = geom.edgeLengths[f.halfedge().edge()];
  double l_jk = geom.edgeLengths[f.halfedge().next().edge()];
  double l_ki = geom.edgeLengths[f.halfedge().next().next().edge()];
  geom.unrequireEdgeLengths();

  double term1 = (uj * vi + ui * vj) * l_ij * l_ij;
  double term2 = (uk * vj + uj * vk) * l_jk * l_jk;
  double term3 = (uk * vi + ui * vk) * l_ki * l_ki;
  return -0.5 * (term1 + term2 + term3);
}

// inline double angle(const BarycentricVector& u, const BarycentricVector& v) {}

inline std::ostream& operator<<(std::ostream& output, const BarycentricVector& v) {

  switch (v.type) {
  case BarycentricVectorType::Face:
    output << "[BarycentricVector: type=Face, face= " << v.face << " faceCoords= " << v.faceCoords << "]";
    break;
  case BarycentricVectorType::Edge:
    output << "[BarycentricVector: type=Edge, edge= " << v.edge << " edgeCoords= " << v.edgeCoords << "]";
    break;
  case BarycentricVectorType::Vertex:
    output << "[BarycentricVector: type=Vertex, vertex= " << v.vertex << "]";
    break;
  }
  return output;
}

} // namespace surface
} // namespace geometrycentral

namespace std {
inline std::string to_string(geometrycentral::surface::BarycentricVector v) {
  ostringstream output;
  output << v;
  return output.str();
}
} // namespace std