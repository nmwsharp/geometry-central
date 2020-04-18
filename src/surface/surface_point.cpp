#include "geometrycentral/surface/surface_point.h"


namespace geometrycentral {
namespace surface {

// == Constructors
SurfacePoint::SurfacePoint() : type(SurfacePointType::Vertex) {}
SurfacePoint::SurfacePoint(Vertex v) : type(SurfacePointType::Vertex), vertex(v) {}
SurfacePoint::SurfacePoint(Edge e, double tEdge_) : type(SurfacePointType::Edge), edge(e), tEdge(tEdge_) {}
SurfacePoint::SurfacePoint(Halfedge he, double tHalfedge)
    : type(SurfacePointType::Edge), edge(he.edge()), tEdge(he == he.edge().halfedge() ? tHalfedge : (1. - tHalfedge)) {}
SurfacePoint::SurfacePoint(Face f, Vector3 faceCoords_)
    : type(SurfacePointType::Face), face(f), faceCoords(faceCoords_) {}


// == Methods

std::ostream& operator<<(std::ostream& output, const SurfacePoint& p) {
  switch (p.type) {
  case SurfacePointType::Vertex: {
    output << "[SurfacePoint: type=Vertex, vertex= " << p.vertex << "]";
    break;
  }
  case SurfacePointType::Edge: {
    output << "[SurfacePoint: type=Edge, edge= " << p.edge << " tEdge= " << p.tEdge << "]";
    break;
  }
  case SurfacePointType::Face: {
    output << "[SurfacePoint: type=Face, face= " << p.face << " faceCoords= " << p.faceCoords << "]";
    break;
  }
  }

  return output;
}


SurfacePoint SurfacePoint::inSomeFace() const {

  switch (type) {
  case SurfacePointType::Vertex: {

    Halfedge he = vertex.halfedge();
    Face targetFace = he.face();
    Halfedge targetHe = targetFace.halfedge();

    // Find the appropriate barycentric coordinate and return
    if (he == targetHe) {
      return SurfacePoint(targetFace, Vector3{1., 0., 0.});
    }
    he = he.next();
    if (he == targetHe) {
      return SurfacePoint(targetFace, Vector3{0., 0., 1.});
    }
    return SurfacePoint(targetFace, Vector3{0., 1., 0.});

    break;
  }
  case SurfacePointType::Edge: {

    Halfedge he = edge.halfedge();
    Face targetFace = he.face();
    Halfedge targetHe = targetFace.halfedge();

    // Find the appropriate barycentric coordinate and return
    if (he == targetHe) {
      return SurfacePoint(targetFace, Vector3{1. - tEdge, tEdge, 0.});
    }
    he = he.next();
    if (he == targetHe) {
      return SurfacePoint(targetFace, Vector3{tEdge, 0., 1. - tEdge});
    }
    return SurfacePoint(targetFace, Vector3{0., 1. - tEdge, tEdge});

    break;
  }
  case SurfacePointType::Face: {
    return *this;
    break;
  }
  }

  throw std::logic_error("bad switch");
  return *this;
}


SurfacePoint SurfacePoint::inFace(Face targetFace) const {

  switch (type) {
  case SurfacePointType::Vertex: {

    Halfedge he = targetFace.halfedge();

    // Find the appropriate barycentric coordinate and return
    if (he.vertex() == vertex) {
      return SurfacePoint(targetFace, Vector3{1., 0., 0.});
    }
    he = he.next();
    if (he.vertex() == vertex) {
      return SurfacePoint(targetFace, Vector3{0., 1., 0.});
    }
    he = he.next();
    if (he.vertex() == vertex) {
      return SurfacePoint(targetFace, Vector3{0., 0., 1.});
    }

    break;
  }

  case SurfacePointType::Edge: {

    double thisT = tEdge;
    for (Halfedge targetHe : {edge.halfedge(), edge.halfedge().twin()}) {

      int i = 0;
      for (Halfedge he : targetFace.adjacentHalfedges()) {
        if (he == targetHe) {
          // Find the appropriate barycentric coordinate and return

          Vector3 bary = Vector3::zero();
          bary[i] = 1.0 - thisT;
          bary[(i + 1) % 3] = thisT;

          return SurfacePoint(targetFace, bary);
        }
        i++;
      }

      // Flip the point to be along the other halfedge
      thisT = 1. - thisT;
    }

    break;
  }

  case SurfacePointType::Face: {
    if (face == targetFace) {
      return *this;
    };
    break;
  }
  }

  throw std::logic_error("SurfacePoint " + std::to_string(*this) + " not adjacent to target face " +
                         std::to_string(targetFace));
  return *this;
}


Vertex SurfacePoint::nearestVertex() const {

  switch (type) {
  case SurfacePointType::Vertex: {
    return vertex;
    break;
  }
  case SurfacePointType::Edge: {
    if (tEdge < .5) return edge.halfedge().vertex();
    return edge.halfedge().twin().vertex();
    break;
  }
  case SurfacePointType::Face: {
    if (faceCoords.x >= faceCoords.y && faceCoords.x >= faceCoords.z) {
      return face.halfedge().vertex();
    }
    if (faceCoords.y >= faceCoords.x && faceCoords.y >= faceCoords.z) {
      return face.halfedge().next().vertex();
    }
    return face.halfedge().next().next().vertex();
    break;
  }
  }

  throw std::logic_error("bad switch");
  return vertex;
}


void SurfacePoint::validate() const {
  const double EPS = 1e-4;

  switch (type) {
  case SurfacePointType::Vertex: {
    if (vertex == Vertex()) throw std::logic_error("surface point with Type::Vertex has invalid vertex ref");
    break;
  }
  case SurfacePointType::Edge: {
    if (edge == Edge()) throw std::logic_error("surface point with Type::Edge has invalid edge ref");
    if (!std::isfinite(tEdge)) throw std::logic_error("surface point with Type::Edge has non-finite tEdge");
    if (tEdge < -EPS || tEdge > (1. + EPS))
      throw std::logic_error("surface point with Type::Edge has tEdge outside of [0,1] = " + std::to_string(tEdge));
    break;
  }
  case SurfacePointType::Face: {
    if (face == Face()) throw std::logic_error("surface point with Type::Face has invalid face ref");
    if (!isfinite(faceCoords)) throw std::logic_error("surface point with Type::Face has non-finite face coords");
    if (faceCoords.x < -EPS || faceCoords.y < -EPS || faceCoords.z < -EPS)
      throw std::logic_error("surface point with Type::Face has negative bary coord " + std::to_string(faceCoords));

    if ((faceCoords.x + faceCoords.y + faceCoords.z) > (1. + EPS))
      throw std::logic_error("surface point with Type::Face has bary coord that sum to > 1 " +
                             std::to_string(faceCoords));
    break;
  }
  }
}

bool SurfacePoint::operator==(const SurfacePoint& other) const {
  if (type != other.type) return false;
  switch (type) {
  case SurfacePointType::Vertex: {
    return vertex == other.vertex;
    break;
  }
  case SurfacePointType::Edge: {
    return edge == other.edge && tEdge == other.tEdge;
    break;
  }
  case SurfacePointType::Face: {
    return face == other.face && faceCoords == other.faceCoords;
    break;
  }
  }
  return false; // should never be reached
}

bool SurfacePoint::operator!=(const SurfacePoint& other) const { return !(*this == other); }


namespace { // helpers

// Helpers for the surface point check
// Next 9 methods are all combinations of {Vertex,Edge,Face} adjacency

bool checkAdjacent(Vertex vA, Vertex vB) {
  bool adjacent = false;
  for (Vertex vN : vA.adjacentVertices()) {
    if (vN == vB) adjacent = true;
  }
  return adjacent;
}

bool checkAdjacent(Vertex vA, Edge eB) {
  bool adjacent = false;
  for (Halfedge he : vA.outgoingHalfedges()) {
    if (eB == he.next().edge()) adjacent = true;
    if (eB == he.edge()) adjacent = true;
  }
  return adjacent;
}

bool checkAdjacent(Vertex vA, Face fB) {
  bool adjacent = false;
  for (Face fN : vA.adjacentFaces()) {
    if (fN == fB) adjacent = true;
  }
  return adjacent;
}

bool checkAdjacent(Edge eA, Vertex vB) { return checkAdjacent(vB, eA); }

bool checkAdjacent(Edge eA, Edge eB) {
  bool adjacent = false;

  // Must have a shared face
  adjacent |= (eA.halfedge().face() == eB.halfedge().face());
  adjacent |= (eA.halfedge().twin().face() == eB.halfedge().face());
  adjacent |= (eA.halfedge().face() == eB.halfedge().twin().face());
  adjacent |= (eA.halfedge().twin().face() == eB.halfedge().twin().face());

  return adjacent;
}

bool checkAdjacent(Edge eA, Face fB) {
  bool adjacent = false;
  for (Edge eN : fB.adjacentEdges()) {
    if (eN == eA) adjacent = true;
  }
  return adjacent;
}


bool checkAdjacent(Face fA, Vertex vB) { return checkAdjacent(vB, fA); }
bool checkAdjacent(Face fA, Edge eB) { return checkAdjacent(eB, fA); }
bool checkAdjacent(Face fA, Face fB) { return fA == fB; }


} // namespace


bool checkAdjacent(const SurfacePoint& pA, const SurfacePoint& pB) {

  switch (pA.type) {
  case SurfacePointType::Vertex:
    switch (pB.type) {
    case SurfacePointType::Vertex:
      return checkAdjacent(pA.vertex, pB.vertex);
      break;
    case SurfacePointType::Edge:
      return checkAdjacent(pA.vertex, pB.edge);
      break;
    case SurfacePointType::Face:
      return checkAdjacent(pA.vertex, pB.face);
      break;
    }
    break;
  case SurfacePointType::Edge:
    switch (pB.type) {
    case SurfacePointType::Vertex:
      return checkAdjacent(pA.edge, pB.vertex);
      break;
    case SurfacePointType::Edge:
      return checkAdjacent(pA.edge, pB.edge);
      break;
    case SurfacePointType::Face:
      return checkAdjacent(pA.edge, pB.face);
      break;
    }
    break;
  case SurfacePointType::Face:
    switch (pB.type) {
    case SurfacePointType::Vertex:
      return checkAdjacent(pA.face, pB.vertex);
      break;
    case SurfacePointType::Edge:
      return checkAdjacent(pA.face, pB.edge);
      break;
    case SurfacePointType::Face:
      return checkAdjacent(pA.face, pB.face);
      break;
    }
    break;
  }

  return false;
}


bool onSameElement(const SurfacePoint& pA, const SurfacePoint& pB) {

  if (pA.type != pB.type) return false;

  switch (pA.type) {
  case SurfacePointType::Vertex:
    return pA.vertex == pB.vertex;
    break;
  case SurfacePointType::Edge:
    return pA.edge == pB.edge;
    break;
  case SurfacePointType::Face:
    return pA.face == pB.face;
    break;
  }

  // unreachable
  throw std::runtime_error("unreachable");
  return false;
}

Face sharedFace(const SurfacePoint& pA, const SurfacePoint& pB) {

  switch (pA.type) {

  case SurfacePointType::Vertex:
    for (Face f : pA.vertex.adjacentFaces()) {
      if (checkAdjacent(SurfacePoint(f, Vector3::zero()), pB)) return f;
    }
    break;

  case SurfacePointType::Edge:

    if (checkAdjacent(SurfacePoint(pA.edge.halfedge().face(), Vector3::zero()), pB)) {
      return pA.edge.halfedge().face();
    }
    if (checkAdjacent(SurfacePoint(pA.edge.halfedge().twin().face(), Vector3::zero()), pB)) {
      return pA.edge.halfedge().twin().face();
    }
    break;

  case SurfacePointType::Face:
    if (checkAdjacent(pA, pB)) return pA.face;
    break;
  }

  // no shared face
  return Face();
}

} // namespace surface
} // namespace geometrycentral


namespace std {
std::string to_string(geometrycentral::surface::SurfacePoint p) {
  ostringstream output;
  output << p;
  return output.str();
}
} // namespace std
