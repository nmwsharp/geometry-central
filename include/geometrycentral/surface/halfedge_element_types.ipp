#pragma once

// Implementations for halfedge_mesh_types.ipp

// Make the element types hashable (this _should_ be doable for just the parent class, but I couldn't sort out how)
namespace std {
// clang-format off
template <> struct hash<geometrycentral::surface::Vertex>         { std::size_t operator()(const geometrycentral::surface::Vertex& e)         const { return std::hash<size_t>{}(e.getIndex()); } };
template <> struct hash<geometrycentral::surface::Halfedge>       { std::size_t operator()(const geometrycentral::surface::Halfedge& e)       const { return std::hash<size_t>{}(e.getIndex()); } };
template <> struct hash<geometrycentral::surface::Corner>         { std::size_t operator()(const geometrycentral::surface::Corner& e)         const { return std::hash<size_t>{}(e.getIndex()); } };
template <> struct hash<geometrycentral::surface::Edge>           { std::size_t operator()(const geometrycentral::surface::Edge& e)           const { return std::hash<size_t>{}(e.getIndex()); } };
template <> struct hash<geometrycentral::surface::Face>           { std::size_t operator()(const geometrycentral::surface::Face& e)           const { return std::hash<size_t>{}(e.getIndex()); } };
template <> struct hash<geometrycentral::surface::BoundaryLoop>   { std::size_t operator()(const geometrycentral::surface::BoundaryLoop& e)   const { return std::hash<size_t>{}(e.getIndex()); } };
// clang-format on
} // namespace std

namespace geometrycentral {
namespace surface {

// clang-format off

// ==========================================================
// ================        Vertex          ==================
// ==========================================================

// Constructors
// (see note in header, these should be inherited but aren't due to compiler issues)
inline Vertex::Vertex() {}
inline Vertex::Vertex(SurfaceMesh* mesh_, size_t ind_) : Element(mesh_,ind_) {}
//inline Vertex::Vertex(const DynamicElement<Vertex>& e) : Element(e.getMesh(), e.getIndex()) {}

// Navigators
inline Halfedge Vertex::halfedge() const    { return Halfedge(mesh, mesh->vHalfedge(ind)); }
inline Corner Vertex::corner() const        { return halfedge().corner(); }
inline bool Vertex::isDead() const          { return mesh->vertexIsDead(ind); }

// Properties
inline bool Vertex::isBoundary() const {
  if(mesh->usesImplicitTwin()) {
    return !halfedge().twin().isInterior();
  } else {
    for (Edge e : adjacentEdges()) {
      if (e.isBoundary()) {return true;}
    }
    return false;
  }
}
inline bool Vertex::isManifoldAndOriented() const {
  // TODO this routine is actually pretty nontrivial, it probably deserves some more thought
  // strategy: bootstrap off of the adjacency iterator, which already has functionality for nonmanifold vertices
  if(mesh->usesImplicitTwin()) return true;
  size_t d = degree(); 
  size_t kManif = 0; 
  Halfedge firstHe = halfedge();
  Halfedge currHe = firstHe;
  do {
    kManif++;
    currHe = currHe.twin().next();
    if(currHe.vertex() != *this) return false;
    if(kManif > d) return false;
  } while(currHe != firstHe);

  return kManif == d;
}
inline bool Vertex::isManifold() const { 
  if(mesh->usesImplicitTwin()) return true;
 
  // test if all incoming edges are manifold
  for(Edge e : adjacentEdges()) {
    if(!e.isManifold()) return false;
  }
  
  // test if there is a single edge-connected component
  std::vector<Face> toProc{halfedge().face()};
  std::unordered_set<Face> seen;
  while(!toProc.empty()) {
    Face f = toProc.back();
    toProc.pop_back();
    if(seen.find(f) != seen.end()) continue;
    seen.insert(f);
    for(Halfedge he : f.adjacentHalfedges()) {
      if(he.vertex()  == *this || he.next().vertex() ==*this) {
        for(Face fn : he.edge().adjacentFaces()) {
          if(seen.find(fn) == seen.end()) toProc.push_back(fn);
        }
      }
    }
  }
  for (Face f : adjacentFaces()) { 
    if(seen.find(f) == seen.end()) return false;
  }

  return true;
}
inline size_t Vertex::degree() const {
  size_t k = 0;
  for (Edge e: adjacentEdges()) { k++; }
  return k;
}
inline size_t Vertex::faceDegree() const {
  size_t k = 0;
  for (Face f : adjacentFaces()) { k++; }
  return k;
}

// Navigation iterators 
inline NavigationSetBase<VertexIncomingHalfedgeNavigator> Vertex::incomingHalfedges() const { 
  return NavigationSetBase<VertexIncomingHalfedgeNavigator>(halfedge().prevOrbitFace()); 
}
inline NavigationSetBase<VertexOutgoingHalfedgeNavigator> Vertex::outgoingHalfedges() const { 
  return NavigationSetBase<VertexOutgoingHalfedgeNavigator>(halfedge()); 
}
inline NavigationSetBase<VertexAdjacentVertexNavigator> Vertex::adjacentVertices() const { 
  return NavigationSetBase<VertexAdjacentVertexNavigator>(VertexNeighborIteratorState(halfedge(), mesh->usesImplicitTwin())); 
}
inline NavigationSetBase<VertexAdjacentFaceNavigator> Vertex::adjacentFaces() const { 
  return NavigationSetBase<VertexAdjacentFaceNavigator>(halfedge()); 
}
inline NavigationSetBase<VertexAdjacentEdgeNavigator> Vertex::adjacentEdges() const { 
  return NavigationSetBase<VertexAdjacentEdgeNavigator>(VertexNeighborIteratorState(halfedge(), mesh->usesImplicitTwin())); 
}
inline NavigationSetBase<VertexAdjacentCornerNavigator> Vertex::adjacentCorners() const {
  return NavigationSetBase<VertexAdjacentCornerNavigator>(halfedge());
}

// == Range iterators
inline bool VertexRangeF::elementOkay(const SurfaceMesh& mesh, size_t ind) {
  return !mesh.vertexIsDead(ind);
}

// == Navigation iterators

// Iterating around vertices on a nonmanifold/nonoriented mesh is hard. 
// Our data structure offers iteration around the outgoing halfedges from each vertex, which takes care of outoing halfedges. Furthermore one thing we know is that for every face in which the vertex appears, there is at least one incoming and one outgoing halfedge. We can use this to find the faces and outgoing halfedges.  For edges and vertices, 
// In both the mannifold and nonmanifold case, if a vertex appears in a face multiple times, (aka its a Delta-complex), then these iterators will return elements multiple times.

inline VertexNeighborIteratorState::VertexNeighborIteratorState(Halfedge currHe_, bool useImplicitTwin_) : useImplicitTwin(useImplicitTwin_), currHe(currHe_), firstHe(currHe_) {}

// clang-format on
inline void VertexNeighborIteratorState::advance() {
  if (useImplicitTwin) {
    currHe = currHe.nextOutgoingNeighbor(); // twin().next()
  } else {
    if (!processingIncoming) {
      // this happens first
      currHe = currHe.nextOutgoingNeighbor();
      if (currHe == firstHe) { // switch to processing incoming if needed
        processingIncoming = true;
        currHe = firstHe.prevOrbitFace();
        firstHe = currHe;
      }
    } else {
      // this happens second
      currHe = currHe.nextIncomingNeighbor();
      if (currHe == firstHe) { // switch back to processing outgoing if needed (returning to initial state)
        processingIncoming = false;
        currHe = firstHe.next();
        firstHe = currHe;
      }
    }
  }
}

inline bool VertexNeighborIteratorState::isHalfedgeCanonical() const {
  // TODO I _think_ that this leads to different Delta-complex behavior on implicit twin vs. without wrt yielding
  // elements multiple times when there is a self-edge...
  if (useImplicitTwin) {
    return true;
  } else {
    return currHe == currHe.edge().halfedge();
  }
}

inline bool VertexNeighborIteratorState::operator==(const VertexNeighborIteratorState& rhs) const {
  return currHe == rhs.currHe && processingIncoming == rhs.processingIncoming;
}
// clang-format off

inline void VertexAdjacentVertexNavigator::advance() { currE.advance(); }
inline bool VertexAdjacentVertexNavigator::isValid() const { return currE.isHalfedgeCanonical(); }
inline Vertex VertexAdjacentVertexNavigator::getCurrent() const {
  if(currE.useImplicitTwin || !currE.processingIncoming) {
    return currE.currHe.next().vertex();
  } else {
    return currE.currHe.vertex();
  }
}

inline void VertexIncomingHalfedgeNavigator::advance() { currE = currE.nextIncomingNeighbor(); }
inline bool VertexIncomingHalfedgeNavigator::isValid() const { return true; }
inline Halfedge VertexIncomingHalfedgeNavigator::getCurrent() const { return currE; }

inline void VertexOutgoingHalfedgeNavigator::advance() {currE = currE.nextOutgoingNeighbor(); }
inline bool VertexOutgoingHalfedgeNavigator::isValid() const { return true; }
inline Halfedge VertexOutgoingHalfedgeNavigator::getCurrent() const { return currE; }

inline void VertexAdjacentCornerNavigator::advance() { currE = currE.nextOutgoingNeighbor(); }
inline bool VertexAdjacentCornerNavigator::isValid() const { return currE.isInterior(); }
inline Corner VertexAdjacentCornerNavigator::getCurrent() const { return currE.corner(); }

inline void VertexAdjacentEdgeNavigator::advance() { currE.advance(); }
inline bool VertexAdjacentEdgeNavigator::isValid() const { return currE.isHalfedgeCanonical(); }
inline Edge VertexAdjacentEdgeNavigator::getCurrent() const { return currE.currHe.edge(); }

inline void VertexAdjacentFaceNavigator::advance() { currE = currE.nextOutgoingNeighbor(); }
inline bool VertexAdjacentFaceNavigator::isValid() const { return currE.isInterior(); }
inline Face VertexAdjacentFaceNavigator::getCurrent() const { return currE.face(); }


// ==========================================================
// ================        Halfedge        ==================
// ==========================================================

// Constructors
inline Halfedge::Halfedge() {}
inline Halfedge::Halfedge(SurfaceMesh* mesh_, size_t ind_) : Element(mesh_,ind_) {}
//inline Halfedge::Halfedge(const DynamicElement<Halfedge>& e) : Element(e.getMesh(), e.getIndex()) {}

// Navigators
inline Halfedge Halfedge::twin() const                  { return Halfedge(mesh, mesh->heSibling(ind)); }
inline Halfedge Halfedge::sibling() const               { return Halfedge(mesh, mesh->heSibling(ind)); }
inline Halfedge Halfedge::next() const                  { return Halfedge(mesh, mesh->heNext(ind)); }
inline Vertex Halfedge::vertex() const                  { return Vertex(mesh, mesh->heVertex(ind)); }
inline Vertex Halfedge::tipVertex() const               { return next().vertex(); }
inline Vertex Halfedge::tailVertex() const              { return vertex(); }
inline Halfedge Halfedge::nextOutgoingNeighbor() const  { return Halfedge(mesh, mesh->heNextOutgoingNeighbor(ind)); }
inline Halfedge Halfedge::nextIncomingNeighbor() const  { return Halfedge(mesh, mesh->heNextIncomingNeighbor(ind)); }
inline Edge Halfedge::edge() const                      { return Edge(mesh, mesh->heEdge(ind)); }
inline Face Halfedge::face() const                      { return Face(mesh, mesh->heFace(ind)); }
inline Corner Halfedge::corner() const                  { return Corner(mesh, ind); }
inline bool Halfedge::isDead() const                    { return mesh->halfedgeIsDead(ind); }


// Super-navigators
inline Halfedge Halfedge::prevOrbitFace() const  { 
  Halfedge currHe = *this;
  while(true) {
    Halfedge nextHe = currHe.next();
    if(nextHe == *this) break;
    currHe = nextHe;
  }
  return currHe;
}
inline Halfedge Halfedge::prevOrbitVertex() const  { 
  Halfedge currHe = twin();
  while(true) {
    Halfedge nextHe = currHe.next();
    if(nextHe == *this) break;
    currHe = nextHe.twin();
  }
  return currHe;
}

// Properties
inline bool Halfedge::isInterior() const { return  mesh->heIsInterior(ind); }
inline bool Halfedge::orientation() const { return  mesh->heOrientation(ind); }

// Range iterators
inline bool HalfedgeRangeF::elementOkay(const SurfaceMesh& mesh, size_t ind) {
  return !mesh.halfedgeIsDead(ind);
}
inline bool HalfedgeInteriorRangeF::elementOkay(const SurfaceMesh& mesh, size_t ind) {
  return !mesh.halfedgeIsDead(ind) && mesh.heIsInterior(ind);
}
inline bool HalfedgeExteriorRangeF::elementOkay(const SurfaceMesh& mesh, size_t ind) {
  return !mesh.halfedgeIsDead(ind) && !mesh.heIsInterior(ind);
}

// ==========================================================
// ================        Corner          ==================
// ==========================================================

// Constructors
inline Corner::Corner() {}
inline Corner::Corner(SurfaceMesh* mesh_, size_t ind_) : Element(mesh_,ind_) {}
//inline Corner::Corner(const DynamicElement<Corner>& e) : Element(e.getMesh(), e.getIndex()) {}

// Navigators
inline Halfedge Corner::halfedge() const { return Halfedge(mesh, ind); }
inline Vertex Corner::vertex() const { return halfedge().vertex(); }
inline Face Corner::face() const { return halfedge().face(); }
inline bool Corner::isDead() const    { return halfedge().isDead(); }

// Range iterators
inline bool CornerRangeF::elementOkay(const SurfaceMesh& mesh, size_t ind) {
  return !mesh.halfedgeIsDead(ind) && mesh.heIsInterior(ind);
}

// ==========================================================
// ================          Edge          ==================
// ==========================================================

// Constructors
inline Edge::Edge() {}
inline Edge::Edge(SurfaceMesh* mesh_, size_t ind_) : Element(mesh_,ind_) {}
//inline Edge::Edge(const DynamicElement<Edge>& e) : Element(e.getMesh(), e.getIndex()) {}

// Navigator
inline Halfedge Edge::halfedge() const      { return Halfedge(mesh, mesh->eHalfedge(ind)); }
inline Vertex Edge::otherVertex(Vertex v) const { 
  if(halfedge().tailVertex() == v) {
    return halfedge().tipVertex();
  } else {
    return halfedge().tailVertex();
  }
}
inline Vertex Edge::firstVertex() const     { return halfedge().tailVertex(); }
inline Vertex Edge::secondVertex() const    { return halfedge().tipVertex(); }
inline std::array<Halfedge,4> Edge::diamondBoundary() const {
   Halfedge h = halfedge();
   Halfedge hN = h.next();
   Halfedge hNN = hN.next();

   Halfedge hT = h.sibling();
   Halfedge hTN = hT.next();
   Halfedge hTNN = hTN.next();

#ifndef NGC_SAFETY_CHECKS
   if(hT == h || !hT.isInterior()) throw std::runtime_error("cannot construct diamondBoundary() of boundary edge");
   if(hT.sibling() != h) throw std::runtime_error("cannot construct diamondBoundary() of nonmanifold edge");
   if(hNN.next() != h) throw std::runtime_error("cannot construct diamondBoundary() for non-triangular face");
   if(hTNN.next() != hT) throw std::runtime_error("cannot construct diamondBoundary() for non-triangular face");
#endif

   return std::array<Halfedge, 4>{hN, hNN, hTN, hTNN};
}
inline bool Edge::isDead() const            { return mesh->edgeIsDead(ind); }

// Properties
inline bool Edge::isBoundary() const { 
  if(mesh->usesImplicitTwin()) {
    return !halfedge().twin().isInterior();
  } else {
    return halfedge().sibling() == halfedge();
  }
} 
inline bool Edge::isManifold() const { return halfedge().sibling().sibling() == halfedge() || halfedge().sibling() == halfedge(); }
inline bool Edge::isOriented() const { 
  Halfedge he = halfedge();
  if(he == he.sibling()) return true; // boundary
  if(he.sibling().sibling() != he) return false; // not manifold
  return he.orientation() != he.sibling().orientation();
}
inline size_t Edge::degree() const { 
  size_t k = 0;
  for(Halfedge he : adjacentInteriorHalfedges()) {
    k++;
  }
  return k;
}

// Range iterators
inline bool EdgeRangeF::elementOkay(const SurfaceMesh& mesh, size_t ind) {
  return !mesh.edgeIsDead(ind);
}

// Navigation iterators

inline NavigationSetBase<EdgeAdjacentHalfedgeNavigator> Edge::adjacentHalfedges() const { 
  return NavigationSetBase<EdgeAdjacentHalfedgeNavigator>(halfedge()); 
}
inline NavigationSetBase<EdgeAdjacentInteriorHalfedgeNavigator> Edge::adjacentInteriorHalfedges() const { 
  return NavigationSetBase<EdgeAdjacentInteriorHalfedgeNavigator>(halfedge()); 
}
inline NavigationSetBase<EdgeAdjacentFaceNavigator> Edge::adjacentFaces() const { 
  return NavigationSetBase<EdgeAdjacentFaceNavigator>(halfedge()); 
}

inline void EdgeAdjacentHalfedgeNavigator::advance() { currE = currE.sibling(); }
inline bool EdgeAdjacentHalfedgeNavigator::isValid() const { return true; }
inline Halfedge EdgeAdjacentHalfedgeNavigator::getCurrent() const { return currE; }

inline void EdgeAdjacentInteriorHalfedgeNavigator::advance() { currE = currE.sibling(); }
inline bool EdgeAdjacentInteriorHalfedgeNavigator::isValid() const { return currE.isInterior(); }
inline Halfedge EdgeAdjacentInteriorHalfedgeNavigator::getCurrent() const { return currE; }

inline void EdgeAdjacentFaceNavigator::advance() { currE = currE.sibling(); }
inline bool EdgeAdjacentFaceNavigator::isValid() const { return currE.isInterior(); }
inline Face EdgeAdjacentFaceNavigator::getCurrent() const { return currE.face(); }


inline std::array<Vertex, 2> Edge::adjacentVertices() const { 
  Halfedge he = halfedge();
  return std::array<Vertex, 2>{he.tailVertex(), he.tipVertex()};
}


// ==========================================================
// ================          Face          ==================
// ==========================================================

// Constructors
inline Face::Face() {}
inline Face::Face(SurfaceMesh* mesh_, size_t ind_) : Element(mesh_,ind_) {}
//inline Face::Face(const DynamicElement<Face>& e) : Element(e.getMesh(), e.getIndex()) {}

// Navigators
inline Halfedge Face::halfedge() const { return Halfedge(mesh, mesh->fHalfedge(ind)); }
inline BoundaryLoop Face::asBoundaryLoop() const { 
  GC_SAFETY_ASSERT(isBoundaryLoop(), "face must be boundary loop to call asBoundaryLoop()")
  return BoundaryLoop(mesh, mesh->faceIndToBoundaryLoopInd(ind)); 
}
inline bool Face::isDead() const    { return mesh->faceIsDead(ind); }

// Properties
inline bool Face::isTriangle() const {
  Halfedge he = halfedge();
  return he == he.next().next().next();
}
inline size_t Face::degree() const {
  size_t k = 0;
  for (Halfedge h : adjacentHalfedges()) { k++; }
  return k;
}

inline bool Face::isBoundaryLoop() const { return mesh->faceIsBoundaryLoop(ind); }

// Navigation iterators
inline NavigationSetBase<FaceAdjacentHalfedgeNavigator> Face::adjacentHalfedges() const { 
  return NavigationSetBase<FaceAdjacentHalfedgeNavigator>(halfedge()); 
}
inline NavigationSetBase<FaceAdjacentVertexNavigator> Face::adjacentVertices() const { 
  return NavigationSetBase<FaceAdjacentVertexNavigator>(halfedge()); 
}
inline NavigationSetBase<FaceAdjacentFaceNavigator> Face::adjacentFaces() const { 
  return NavigationSetBase<FaceAdjacentFaceNavigator>(std::make_pair(halfedge(), halfedge())); 
}
inline NavigationSetBase<FaceAdjacentEdgeNavigator> Face::adjacentEdges() const { 
  return NavigationSetBase<FaceAdjacentEdgeNavigator>(halfedge()); 
}
inline NavigationSetBase<FaceAdjacentCornerNavigator> Face::adjacentCorners() const { 
  return NavigationSetBase<FaceAdjacentCornerNavigator>(halfedge()); 
}


// == Range iterators
inline bool FaceRangeF::elementOkay(const SurfaceMesh& mesh, size_t ind) {
  return !mesh.faceIsDead(ind);
}

// == Navigation iterators

inline void FaceAdjacentVertexNavigator::advance() { currE = currE.next(); }
inline bool FaceAdjacentVertexNavigator::isValid() const { return true; }
inline Vertex FaceAdjacentVertexNavigator::getCurrent() const { return currE.vertex(); }

inline void FaceAdjacentHalfedgeNavigator::advance() { currE = currE.next(); }
inline bool FaceAdjacentHalfedgeNavigator::isValid() const { return true; }
inline Halfedge FaceAdjacentHalfedgeNavigator::getCurrent() const { return currE; }

inline void FaceAdjacentCornerNavigator::advance() { currE = currE.next(); }
inline bool FaceAdjacentCornerNavigator::isValid() const { return true; }
inline Corner FaceAdjacentCornerNavigator::getCurrent() const { return currE.corner(); }

inline void FaceAdjacentEdgeNavigator::advance() { currE = currE.next(); }
inline bool FaceAdjacentEdgeNavigator::isValid() const { return true; }
inline Edge FaceAdjacentEdgeNavigator::getCurrent() const { return currE.edge(); }

inline void FaceAdjacentFaceNavigator::advance() { 
  currE.second = currE.second.sibling();
  if(currE.first == currE.second) {
    currE.first = currE.first.next(); 
    currE.second = currE.first;
  }
}
inline bool FaceAdjacentFaceNavigator::isValid() const { return currE.first != currE.second && currE.second.isInterior(); }
inline Face FaceAdjacentFaceNavigator::getCurrent() const { return currE.second.face(); }


// ==========================================================
// ================     Boundary Loop      ==================
// ==========================================================

// Constructors
inline BoundaryLoop::BoundaryLoop() {}
inline BoundaryLoop::BoundaryLoop(SurfaceMesh* mesh_, size_t ind_) : Element(mesh_,ind_) {}
//inline BoundaryLoop::BoundaryLoop(const DynamicElement<BoundaryLoop>& e) : Element(e.getMesh(), e.getIndex()) {}


inline Halfedge BoundaryLoop::halfedge() const { return asFace().halfedge(); }
inline Face BoundaryLoop::asFace() const { return Face(mesh, mesh->boundaryLoopIndToFaceInd(ind)); }
inline bool BoundaryLoop::isDead() const    { return asFace().isDead(); }

inline size_t BoundaryLoop::degree() const {
  size_t k = 0;
  for (Halfedge h : adjacentHalfedges()) { k++; }
  return k;
}

inline NavigationSetBase<BoundaryLoopAdjacentHalfedgeNavigator> BoundaryLoop::adjacentHalfedges() const { 
  return NavigationSetBase<BoundaryLoopAdjacentHalfedgeNavigator>(halfedge()); 
}
inline NavigationSetBase<BoundaryLoopAdjacentVertexNavigator> BoundaryLoop::adjacentVertices() const { 
  return NavigationSetBase<BoundaryLoopAdjacentVertexNavigator>(halfedge()); 
}
inline NavigationSetBase<BoundaryLoopAdjacentEdgeNavigator> BoundaryLoop::adjacentEdges() const { 
  return NavigationSetBase<BoundaryLoopAdjacentEdgeNavigator>(halfedge()); 
}

// == Range iterators
inline bool BoundaryLoopRangeF::elementOkay(const SurfaceMesh& mesh, size_t ind) {
  return !mesh.faceIsDead(mesh.boundaryLoopIndToFaceInd(ind));
}

// == Navigation iterators

inline void BoundaryLoopAdjacentVertexNavigator::advance() { currE = currE.next(); }
inline bool BoundaryLoopAdjacentVertexNavigator::isValid() const { return true; }
inline Vertex BoundaryLoopAdjacentVertexNavigator::getCurrent() const { return currE.vertex(); }

inline void BoundaryLoopAdjacentHalfedgeNavigator::advance() { currE = currE.next(); }
inline bool BoundaryLoopAdjacentHalfedgeNavigator::isValid() const { return true; }
inline Halfedge BoundaryLoopAdjacentHalfedgeNavigator::getCurrent() const { return currE; }

inline void BoundaryLoopAdjacentEdgeNavigator::advance() { currE = currE.next(); }
inline bool BoundaryLoopAdjacentEdgeNavigator::isValid() const { return true; }
inline Edge BoundaryLoopAdjacentEdgeNavigator::getCurrent() const { return currE.edge(); }


// clang-format on

} // namespace surface
} // namespace geometrycentral
