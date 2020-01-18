#include "geometrycentral/surface/surface_point.h"


namespace geometrycentral {
namespace surface {

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


} // namespace surface
} // namespace geometrycentral
