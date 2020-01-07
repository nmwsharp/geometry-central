#pragma once

// Implementations for halfedge_mesh_types.ipp

namespace geometrycentral {
namespace surface {

// clang-format off

// ==========================================================
// ================      Base Element      ==================
// ==========================================================

// Constructors
template<typename T> 
Element<T>::Element() {}
template<typename T> 
Element<T>::Element(HalfedgeMesh* mesh_, size_t ind_) : mesh(mesh_), ind(ind_) {}
template<typename T> 
Element<T>::Element(const DynamicElement<T>& e) : mesh(e.getMesh()), ind(e.getIndex()) {}

// Comparators
template<typename T> 
inline bool Element<T>::operator==(const Element<T>& other) const { return ind == other.ind; }
template<typename T> 
inline bool Element<T>::operator!=(const Element<T>& other) const { return !(*this == other); }
template<typename T> 
inline bool Element<T>::operator>(const Element<T>& other) const { return ind > other.ind; }
template<typename T> 
inline bool Element<T>::operator>=(const Element<T>& other) const { return ind >= other.ind; }
template<typename T> 
inline bool Element<T>::operator<(const Element<T>& other) const { return ind < other.ind; }
template<typename T> 
inline bool Element<T>::operator<=(const Element<T>& other) const { return ind <= other.ind; }

template <typename T>
size_t Element<T>::getIndex() const { return ind; }

template <typename T>
HalfedgeMesh* Element<T>::getMesh() const { return mesh; }

template <typename T>
inline ::std::ostream& operator<<(::std::ostream& output, const Element<T>& e) {
  output << typeShortName<T>() << "_" << e.ind;
  return output;
}

// Dynamic element
template<typename S> 
DynamicElement<S>::DynamicElement() {}

template<typename S> 
DynamicElement<S>::DynamicElement(HalfedgeMesh* mesh_, size_t ind_) : S(mesh_, ind_) {
  registerWithMesh();
}

template<typename S> 
DynamicElement<S>::DynamicElement(const S& e) : S(e) {
  registerWithMesh();
}
  
template <typename S>
DynamicElement<S>::DynamicElement(const DynamicElement& other) : S(other.mesh, other.ind) {
  registerWithMesh();
}

template <typename S>
DynamicElement<S>::DynamicElement(DynamicElement&& other) : S(other.mesh, other.ind) {
  registerWithMesh();
}

template <typename S>
DynamicElement<S>& DynamicElement<S>::operator=(const DynamicElement<S>& other) {
  deregisterWithMesh();
  this->mesh = other.mesh;
  this->ind = other.ind;
  registerWithMesh();
  return *this;
}

template <typename S>
DynamicElement<S>& DynamicElement<S>::operator=(DynamicElement<S>&& other) noexcept {
  deregisterWithMesh();
  this->mesh = other.mesh;
  this->ind = other.ind;
  registerWithMesh();
  return *this;
}

template<typename S> 
DynamicElement<S>::~DynamicElement() {
  deregisterWithMesh();
}

template<typename S> 
void DynamicElement<S>::registerWithMesh() {
  
  // Callback function on permutation
  std::function<void(const std::vector<size_t>&)> permuteFunc = [this](const std::vector<size_t>& perm) {
    // TODO FIXME not implemented. See note in mesh compression callbacks.
  };

  // Callback function on mesh delete
  std::function<void()> deleteFunc = [this]() {
    // Ensures that we don't try to remove with iterators on deconstruct of this object
    this->mesh = nullptr;
  };

  permuteCallbackIt = getPermuteCallbackList<S>(this->mesh).insert(getPermuteCallbackList<S>(this->mesh).end(), permuteFunc);
  deleteCallbackIt = this->mesh->meshDeleteCallbackList.insert(this->mesh->meshDeleteCallbackList.end(), deleteFunc);
}

template<typename S> 
void DynamicElement<S>::deregisterWithMesh() {
  if (this->mesh == nullptr) return;
  getPermuteCallbackList<S>(this->mesh).erase(permuteCallbackIt);
  this->mesh->meshDeleteCallbackList.erase(deleteCallbackIt);
}


// Base iterators
template <typename F>
inline RangeIteratorBase<F>::RangeIteratorBase(HalfedgeMesh* mesh_, size_t iStart_, size_t iEnd_) : mesh(mesh_), iCurr(iStart_), iEnd(iEnd_) {
  if (iCurr != iEnd && !F::elementOkay(*mesh, iCurr)) {
    this->operator++();
  }
}

template <typename F>
inline const RangeIteratorBase<F>& RangeIteratorBase<F>::operator++() {
  iCurr++;
  while (iCurr != iEnd && !F::elementOkay(*mesh, iCurr)) {
    iCurr++;
  }
  return *this;
}

template <typename F>
inline bool RangeIteratorBase<F>::operator==(const RangeIteratorBase<F>& other) const {
	return iCurr == other.iCurr;
}

template <typename F>
inline bool RangeIteratorBase<F>::operator!=(const RangeIteratorBase<F>& other) const {
	return !(*this == other);
}

template <typename F>
inline typename F::Etype RangeIteratorBase<F>::operator*() const { return typename F::Etype(mesh, iCurr); }

template <typename F>
RangeSetBase<F>::RangeSetBase(HalfedgeMesh* mesh_, size_t iStart_, size_t iEnd_) : mesh(mesh_), iStart(iStart_), iEnd(iEnd_) {}  

template <typename F>
inline RangeIteratorBase<F> RangeSetBase<F>::begin() const { return RangeIteratorBase<F>(mesh, iStart, iEnd); }

template <typename F>
inline RangeIteratorBase<F> RangeSetBase<F>::end() const { return RangeIteratorBase<F>(mesh, iEnd, iEnd); }

// ==========================================================
// ================        Vertex          ==================
// ==========================================================

// Constructors
// (see note in header, these should be inherited but aren't due to compiler issues)
inline Vertex::Vertex() {}
inline Vertex::Vertex(HalfedgeMesh* mesh_, size_t ind_) : Element(mesh_,ind_) {}
inline Vertex::Vertex(const DynamicElement<Vertex>& e) : Element(e.getMesh(), e.getIndex()) {}

// Navigators
inline Halfedge Vertex::halfedge() const    { return Halfedge(mesh, mesh->vHalfedge[ind]); }
inline Corner Vertex::corner() const        { return halfedge().corner(); }
inline bool Vertex::isDead() const          { return mesh->vertexIsDead(ind); }

// Properties
inline bool Vertex::isBoundary() const { return !halfedge().twin().isInterior(); }
inline size_t Vertex::degree() const {
  size_t k = 0;
  for (Halfedge h : outgoingHalfedges()) { k++; }
  return k;
}
inline size_t Vertex::faceDegree() const {
  size_t d = degree();
  if(isBoundary()) {
    return d - 1;
  } else {
    return d;
  }
}

// Navigation iterators 
inline NavigationSetBase<VertexIncomingHalfedgeNavigator> Vertex::incomingHalfedges() const { 
  return NavigationSetBase<VertexIncomingHalfedgeNavigator>(halfedge()); 
}
inline NavigationSetBase<VertexOutgoingHalfedgeNavigator> Vertex::outgoingHalfedges() const { 
  return NavigationSetBase<VertexOutgoingHalfedgeNavigator>(halfedge()); 
}
inline NavigationSetBase<VertexAdjacentVertexNavigator> Vertex::adjacentVertices() const { 
  return NavigationSetBase<VertexAdjacentVertexNavigator>(halfedge()); 
}
inline NavigationSetBase<VertexAdjacentFaceNavigator> Vertex::adjacentFaces() const { 
  return NavigationSetBase<VertexAdjacentFaceNavigator>(halfedge()); 
}
inline NavigationSetBase<VertexAdjacentEdgeNavigator> Vertex::adjacentEdges() const { 
  return NavigationSetBase<VertexAdjacentEdgeNavigator>(halfedge()); 
}
inline NavigationSetBase<VertexAdjacentCornerNavigator> Vertex::adjacentCorners() const {
  return NavigationSetBase<VertexAdjacentCornerNavigator>(halfedge());
}

// Range iterators
inline bool VertexRangeF::elementOkay(const HalfedgeMesh& mesh, size_t ind) {
  return !mesh.vertexIsDead(ind);
}

// ==========================================================
// ================        Halfedge        ==================
// ==========================================================

// Constructors
inline Halfedge::Halfedge() {}
inline Halfedge::Halfedge(HalfedgeMesh* mesh_, size_t ind_) : Element(mesh_,ind_) {}
inline Halfedge::Halfedge(const DynamicElement<Halfedge>& e) : Element(e.getMesh(), e.getIndex()) {}

// Navigators
inline Halfedge Halfedge::twin() const  { return Halfedge(mesh, HalfedgeMesh::heTwin(ind)); }
inline Halfedge Halfedge::next() const  { return Halfedge(mesh, mesh->heNext[ind]); }
inline Vertex Halfedge::vertex() const  { return Vertex(mesh, mesh->heVertex[ind]); }
inline Edge Halfedge::edge() const      { return Edge(mesh, HalfedgeMesh::heEdge(ind)); }
inline Face Halfedge::face() const      { return Face(mesh, mesh->heFace[ind]); }
inline Corner Halfedge::corner() const  { return Corner(mesh, ind); }
inline bool Halfedge::isDead() const    { return mesh->halfedgeIsDead(ind); }


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

// Range iterators
inline bool HalfedgeRangeF::elementOkay(const HalfedgeMesh& mesh, size_t ind) {
  return !mesh.halfedgeIsDead(ind);
}
inline bool HalfedgeInteriorRangeF::elementOkay(const HalfedgeMesh& mesh, size_t ind) {
  return !mesh.halfedgeIsDead(ind) && mesh.heIsInterior(ind);
}
inline bool HalfedgeExteriorRangeF::elementOkay(const HalfedgeMesh& mesh, size_t ind) {
  return !mesh.halfedgeIsDead(ind) && !mesh.heIsInterior(ind);
}

// ==========================================================
// ================        Corner          ==================
// ==========================================================

// Constructors
inline Corner::Corner() {}
inline Corner::Corner(HalfedgeMesh* mesh_, size_t ind_) : Element(mesh_,ind_) {}
inline Corner::Corner(const DynamicElement<Corner>& e) : Element(e.getMesh(), e.getIndex()) {}

// Navigators
inline Halfedge Corner::halfedge() const { return Halfedge(mesh, ind); }
inline Vertex Corner::vertex() const { return halfedge().vertex(); }
inline Face Corner::face() const { return halfedge().face(); }
inline bool Corner::isDead() const    { return halfedge().isDead(); }

// Range iterators
inline bool CornerRangeF::elementOkay(const HalfedgeMesh& mesh, size_t ind) {
  return !mesh.halfedgeIsDead(ind) && mesh.heIsInterior(ind);
}

// ==========================================================
// ================          Edge          ==================
// ==========================================================

// Constructors
inline Edge::Edge() {}
inline Edge::Edge(HalfedgeMesh* mesh_, size_t ind_) : Element(mesh_,ind_) {}
inline Edge::Edge(const DynamicElement<Edge>& e) : Element(e.getMesh(), e.getIndex()) {}

// Navigators
inline Halfedge Edge::halfedge() const { return Halfedge(mesh, HalfedgeMesh::eHalfedge(ind)); }
inline bool Edge::isDead() const    { return mesh->edgeIsDead(ind); }

// Properties
inline bool Edge::isBoundary() const { return !halfedge().isInterior() || !halfedge().twin().isInterior(); }

// Range iterators
inline bool EdgeRangeF::elementOkay(const HalfedgeMesh& mesh, size_t ind) {
  return !mesh.edgeIsDead(ind);
}


// ==========================================================
// ================          Face          ==================
// ==========================================================

// Constructors
inline Face::Face() {}
inline Face::Face(HalfedgeMesh* mesh_, size_t ind_) : Element(mesh_,ind_) {}
inline Face::Face(const DynamicElement<Face>& e) : Element(e.getMesh(), e.getIndex()) {}

// Navigators
inline Halfedge Face::halfedge() const { return Halfedge(mesh, mesh->fHalfedge[ind]); }
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
  return NavigationSetBase<FaceAdjacentFaceNavigator>(halfedge()); 
}
inline NavigationSetBase<FaceAdjacentEdgeNavigator> Face::adjacentEdges() const { 
  return NavigationSetBase<FaceAdjacentEdgeNavigator>(halfedge()); 
}
inline NavigationSetBase<FaceAdjacentCornerNavigator> Face::adjacentCorners() const { 
  return NavigationSetBase<FaceAdjacentCornerNavigator>(halfedge()); 
}


// Range iterators
inline bool FaceRangeF::elementOkay(const HalfedgeMesh& mesh, size_t ind) {
  return !mesh.faceIsDead(ind);
}

// ==========================================================
// ================     Boundary Loop      ==================
// ==========================================================

// Constructors
inline BoundaryLoop::BoundaryLoop() {}
inline BoundaryLoop::BoundaryLoop(HalfedgeMesh* mesh_, size_t ind_) : Element(mesh_,ind_) {}
inline BoundaryLoop::BoundaryLoop(const DynamicElement<BoundaryLoop>& e) : Element(e.getMesh(), e.getIndex()) {}


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

// Range iterators
inline bool BoundaryLoopRangeF::elementOkay(const HalfedgeMesh& mesh, size_t ind) {
  return !mesh.faceIsDead(mesh.boundaryLoopIndToFaceInd(ind));
}


// clang-format on

} // namespace surface
} // namespace geometrycentral

namespace std {

// For lookup reasons I don't entirely understand, need to list these out explicitly, template on base does not resolve

// clang-format off
template <> struct hash<geometrycentral::surface::Vertex>         { std::size_t operator()(const geometrycentral::surface::Vertex& e)         const { return std::hash<size_t>{}(e.getIndex()); } };
template <> struct hash<geometrycentral::surface::Halfedge>       { std::size_t operator()(const geometrycentral::surface::Halfedge& e)       const { return std::hash<size_t>{}(e.getIndex()); } };
template <> struct hash<geometrycentral::surface::Corner>         { std::size_t operator()(const geometrycentral::surface::Corner& e)         const { return std::hash<size_t>{}(e.getIndex()); } };
template <> struct hash<geometrycentral::surface::Edge>           { std::size_t operator()(const geometrycentral::surface::Edge& e)           const { return std::hash<size_t>{}(e.getIndex()); } };
template <> struct hash<geometrycentral::surface::Face>           { std::size_t operator()(const geometrycentral::surface::Face& e)           const { return std::hash<size_t>{}(e.getIndex()); } };
template <> struct hash<geometrycentral::surface::BoundaryLoop>   { std::size_t operator()(const geometrycentral::surface::BoundaryLoop& e)   const { return std::hash<size_t>{}(e.getIndex()); } };

// clang-format on

} // namespace std
