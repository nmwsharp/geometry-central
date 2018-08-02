#pragma once

// Implementations for halfedge_mesh_types.ipp

namespace geometrycentral {

// Halfedge
inline HalfedgePtr::HalfedgePtr() : ptr(nullptr) {}
inline HalfedgePtr::HalfedgePtr(Halfedge* ptr_) : ptr(ptr_) {}
inline HalfedgePtr::HalfedgePtr(DynamicHalfedgePtr ptr_) : HalfedgePtr(ptr_.mesh->halfedge(ptr_.ind)) {}
inline HalfedgePtr HalfedgePtr::twin() const { return ptr->twin; }
inline HalfedgePtr HalfedgePtr::next() const { return ptr->next; }
inline HalfedgePtr HalfedgePtr::prev() const {
  Halfedge* h = ptr;
  while (h->next != ptr) {
    h = h->next;
  }
  return HalfedgePtr{h};
}
inline VertexPtr HalfedgePtr::vertex() const { return ptr->vertex; }
inline EdgePtr HalfedgePtr::edge() const { return ptr->edge; }
inline FacePtr HalfedgePtr::face() const { return ptr->face; }
inline CornerPtr HalfedgePtr::corner() const { return CornerPtr{ptr}; }
inline bool HalfedgePtr::isReal() const { return ptr->isReal; }
inline Halfedge HalfedgePtr::operator*() { return *ptr; }
inline Halfedge HalfedgePtr::operator*() const { return *ptr; }
inline Halfedge* HalfedgePtr::operator->() { return ptr; }
inline const Halfedge* HalfedgePtr::operator->() const { return ptr; }
inline bool HalfedgePtr::operator==(const HalfedgePtr& other) const {
  return ptr == other.ptr;
}
inline bool HalfedgePtr::operator!=(const HalfedgePtr& other) const {
  return !(*this == other);
}
inline bool HalfedgePtr::operator>(const HalfedgePtr& other) const {
  return ptr > other.ptr;
}
inline bool HalfedgePtr::operator>=(const HalfedgePtr& other) const {
  return ptr >= other.ptr;
}
inline bool HalfedgePtr::operator<(const HalfedgePtr& other) const {
  return ptr < other.ptr;
}
inline bool HalfedgePtr::operator<=(const HalfedgePtr& other) const {
  return ptr <= other.ptr;
}
inline bool HalfedgePtr::operator==(void* n) const {
  return ptr == n;
}
inline bool HalfedgePtr::operator!=(void* n) const {
  return ptr != n;
}
inline bool HalfedgePtr::operator>(void* n) const { return ptr > n; }
inline bool HalfedgePtr::operator>=(void* n) const {
  return ptr >= n;
}
inline bool HalfedgePtr::operator<(void* n) const { return ptr < n; }
inline bool HalfedgePtr::operator<=(void* n) const {
  return ptr <= n;
}
inline unsigned int HalfedgePtr::operator-(const HalfedgePtr& other) const {
#ifndef NDEBUG
  assert(ptr >= other.ptr &&
         "Pointer subtraction must yield a nonnegative result");
#endif
  return ptr - other.ptr;
}
inline HalfedgePtr& HalfedgePtr::operator++() {
  ptr++;
  return *this;
}
inline HalfedgePtr HalfedgePtr::operator++(int) {
  ptr++;
  return HalfedgePtr(ptr - 1);
}
inline HalfedgePtr& HalfedgePtr::operator--() {
  ptr--;
  return *this;
}
inline HalfedgePtr HalfedgePtr::operator--(int) {
  ptr--;
  return HalfedgePtr(ptr + 1);
}

inline ::std::ostream& operator<<(::std::ostream& output,
                                  const HalfedgePtr& he) {
  output << "he_" << he.ptr;
  return output;
}

inline DynamicHalfedgePtr::DynamicHalfedgePtr(Halfedge* ptr_, HalfedgeMesh* mesh_) : ind(mesh_->indexOf(ptr_)), mesh(mesh_) { }
inline DynamicHalfedgePtr::DynamicHalfedgePtr(HalfedgePtr ptr_, HalfedgeMesh* mesh_) : ind(mesh_->indexOf(ptr_.ptr)), mesh(mesh_) {
}
inline bool DynamicHalfedgePtr::operator==(const DynamicHalfedgePtr& other) const {
  return ind == other.ind;
}
inline bool DynamicHalfedgePtr::operator!=(const DynamicHalfedgePtr& other) const {
  return !(*this == other);
}
inline bool DynamicHalfedgePtr::operator>(const DynamicHalfedgePtr& other) const {
  return ind > other.ind;
}
inline bool DynamicHalfedgePtr::operator>=(const DynamicHalfedgePtr& other) const {
  return ind >= other.ind;
}
inline bool DynamicHalfedgePtr::operator<(const DynamicHalfedgePtr& other) const {
  return ind < other.ind;
}
inline bool DynamicHalfedgePtr::operator<=(const DynamicHalfedgePtr& other) const {
  return ind <= other.ind;
}
inline size_t DynamicHalfedgePtr::getInd() const {return ind;}

inline HalfedgePtrRangeIterator::HalfedgePtrRangeIterator(
    HalfedgePtr startingHalfedge, HalfedgeSetType type_, HalfedgePtr end_)
    : currHalfedge(startingHalfedge), type(type_), end(end_) {}
inline const HalfedgePtrRangeIterator& HalfedgePtrRangeIterator::operator++() {
  currHalfedge++;

  // = Respect the 'type' option
  // Note that we always need to return if we fall off the end of the list, and must not dereference the pointer in that case
 
  // All halfedges
  if(type == HalfedgeSetType::All || currHalfedge == end) {
    return *this;
  }
  // Real only
  else if(type == HalfedgeSetType::Real) {
    while(currHalfedge != end && !currHalfedge.isReal()) {
      currHalfedge++;
    }
    return *this;
  } 
  // Imaginary only
  else /* imag */ {
    while(currHalfedge != end && currHalfedge.isReal()) {
      currHalfedge++;
    }
    return *this;
  }

}
inline bool HalfedgePtrRangeIterator::operator==(
    const HalfedgePtrRangeIterator& other) const {
  return currHalfedge == other.currHalfedge;
}
inline bool HalfedgePtrRangeIterator::operator!=(
    const HalfedgePtrRangeIterator& other) const {
  return currHalfedge != other.currHalfedge;
}
inline HalfedgePtr HalfedgePtrRangeIterator::operator*() const {
  return currHalfedge;
}

inline HalfedgePtrSet::HalfedgePtrSet(HalfedgePtr beginptr_,
                                      HalfedgePtr endptr_, HalfedgeSetType type_)
    : beginptr(beginptr_), endptr(endptr_), type(type_) {}
inline HalfedgePtrRangeIterator HalfedgePtrSet::begin() {
  return HalfedgePtrRangeIterator(beginptr, type, endptr);
}
inline HalfedgePtrRangeIterator HalfedgePtrSet::end() {
  return HalfedgePtrRangeIterator(endptr, type, endptr);
}


// Corner
inline CornerPtr::CornerPtr() : ptr(nullptr) {}
inline CornerPtr::CornerPtr(Halfedge* ptr_) : ptr(ptr_) {}
inline CornerPtr CornerPtr::next() const { return halfedge().next().corner(); }
inline CornerPtr CornerPtr::prev() const { return halfedge().prev().corner(); }
inline HalfedgePtr CornerPtr::halfedge() const { return HalfedgePtr{ptr}; }
inline VertexPtr CornerPtr::vertex() const {
  return halfedge().prev().vertex();
}
inline FacePtr CornerPtr::face() const { return halfedge().face(); }
inline Halfedge CornerPtr::operator*() { return *ptr; }
inline Halfedge CornerPtr::operator*() const { return *ptr; }
inline Halfedge* CornerPtr::operator->() { return ptr; }
inline const Halfedge* CornerPtr::operator->() const { return ptr; }
inline bool CornerPtr::operator==(const CornerPtr& other) const {
  return ptr == other.ptr;
}
inline bool CornerPtr::operator!=(const CornerPtr& other) const {
  return !(*this == other);
}
inline bool CornerPtr::operator>(const CornerPtr& other) const {
  return ptr > other.ptr;
}
inline bool CornerPtr::operator>=(const CornerPtr& other) const {
  return ptr >= other.ptr;
}
inline bool CornerPtr::operator<(const CornerPtr& other) const {
  return ptr < other.ptr;
}
inline bool CornerPtr::operator<=(const CornerPtr& other) const {
  return ptr <= other.ptr;
}
inline bool CornerPtr::operator==(void* n) const { return ptr == n; }
inline bool CornerPtr::operator!=(void* n) const { return ptr != n; }
inline bool CornerPtr::operator>(void* n) const { return ptr > n; }
inline bool CornerPtr::operator>=(void* n) const { return ptr >= n; }
inline bool CornerPtr::operator<(void* n) const { return ptr < n; }
inline bool CornerPtr::operator<=(void* n) const { return ptr <= n; }
inline unsigned int CornerPtr::operator-(const CornerPtr& other) const {
#ifndef NDEBUG
  assert(ptr >= other.ptr &&
         "Pointer subtraction must yield a nonnegative result");
#endif
  return ptr - other.ptr;
}
inline CornerPtr& CornerPtr::operator++() {
  ptr++;
  return *this;
}
inline CornerPtr CornerPtr::operator++(int) {
  ptr++;
  return CornerPtr(ptr - 1);
}
inline CornerPtr& CornerPtr::operator--() {
  ptr--;
  return *this;
}
inline CornerPtr CornerPtr::operator--(int) {
  ptr--;
  return CornerPtr(ptr + 1);
}

inline CornerPtrRangeIterator::CornerPtrRangeIterator(CornerPtr startingCorner)
    : currCorner(startingCorner) {}
inline const CornerPtrRangeIterator& CornerPtrRangeIterator::operator++() {
  currCorner++;
  return *this;
}
inline bool CornerPtrRangeIterator::operator==(
    const CornerPtrRangeIterator& other) const {
  return currCorner == other.currCorner;
}
inline bool CornerPtrRangeIterator::operator!=(
    const CornerPtrRangeIterator& other) const {
  return currCorner != other.currCorner;
}
inline CornerPtr CornerPtrRangeIterator::operator*() const {
  return currCorner;
}

inline CornerPtrSet::CornerPtrSet(CornerPtr beginptr_, CornerPtr endptr_)
    : beginptr(beginptr_), endptr(endptr_) {}
inline CornerPtrRangeIterator CornerPtrSet::begin() {
  return CornerPtrRangeIterator(beginptr);
}
inline CornerPtrRangeIterator CornerPtrSet::end() {
  return CornerPtrRangeIterator(endptr);
}

// Vertex
inline VertexPtr::VertexPtr() : ptr(nullptr) {}
inline VertexPtr::VertexPtr(Vertex* ptr_) : ptr(ptr_) {}
inline VertexPtr::VertexPtr(DynamicVertexPtr ptr_) : VertexPtr(ptr_.mesh->vertex(ptr_.ind)) {}
inline HalfedgePtr VertexPtr::halfedge() const { return ptr->halfedge; }
inline CornerPtr VertexPtr::corner() const {
  HalfedgePtr h = halfedge();
  if (!h.isReal()) h = h.twin().next();
  return h.next().corner();
}
inline bool VertexPtr::isBoundary() const { return ptr->isBoundary; }
inline VertexIncomingHalfedgeSet VertexPtr::incomingHalfedges() {
  return VertexIncomingHalfedgeSet(halfedge().twin());
}
inline VertexOutgoingHalfedgeSet VertexPtr::outgoingHalfedges() {
  return VertexOutgoingHalfedgeSet(halfedge());
}
inline VertexIncomingInteriorHalfedgeSet
VertexPtr::incomingInteriorHalfedges() {
  return VertexIncomingInteriorHalfedgeSet(halfedge().twin());
}
inline VertexOutgoingInteriorHalfedgeSet
VertexPtr::outgoingInteriorHalfedges() {
  return VertexOutgoingInteriorHalfedgeSet(halfedge());
}
inline VertexAdjacentVertexSet VertexPtr::adjacentVertices() {
  return VertexAdjacentVertexSet(halfedge().twin());
}
inline VertexAdjacentFaceSet VertexPtr::adjacentFaces() {
  return VertexAdjacentFaceSet(halfedge());
}
inline VertexAdjacentEdgeSet VertexPtr::adjacentEdges() {
  return VertexAdjacentEdgeSet(halfedge());
}
inline VertexAdjacentCornerSet VertexPtr::adjacentCorners() {
  HalfedgePtr h = halfedge();
  if (!h.isReal()) h = h.twin().next();
  return VertexAdjacentCornerSet(h);
}
inline Vertex VertexPtr::operator*() { return *ptr; }
inline Vertex VertexPtr::operator*() const { return *ptr; }
inline Vertex* VertexPtr::operator->() { return ptr; }
inline const Vertex* VertexPtr::operator->() const { return ptr; }
inline bool VertexPtr::operator==(const VertexPtr& other) const {
  return ptr == other.ptr;
}
inline bool VertexPtr::operator!=(const VertexPtr& other) const {
  return !(*this == other);
}
inline bool VertexPtr::operator>(const VertexPtr& other) const {
  return ptr > other.ptr;
}
inline bool VertexPtr::operator>=(const VertexPtr& other) const {
  return ptr >= other.ptr;
}
inline bool VertexPtr::operator<(const VertexPtr& other) const {
  return ptr < other.ptr;
}
inline bool VertexPtr::operator<=(const VertexPtr& other) const {
  return ptr <= other.ptr;
}
inline bool VertexPtr::operator==(void* n) const { return ptr == n; }
inline bool VertexPtr::operator!=(void* n) const { return ptr != n; }
inline bool VertexPtr::operator>(void* n) const { return ptr > n; }
inline bool VertexPtr::operator>=(void* n) const { return ptr >= n; }
inline bool VertexPtr::operator<(void* n) const { return ptr < n; }
inline bool VertexPtr::operator<=(void* n) const { return ptr <= n; }
inline unsigned int VertexPtr::operator-(const VertexPtr& other) const {
#ifndef NDEBUG
  assert(ptr >= other.ptr &&
         "Pointer subtraction must yield a nonnegative result");
#endif
  return ptr - other.ptr;
}
inline VertexPtr& VertexPtr::operator++() {
  ptr++;
  return *this;
}
inline VertexPtr VertexPtr::operator++(int) {
  ptr++;
  return VertexPtr(ptr - 1);
}
inline VertexPtr& VertexPtr::operator--() {
  ptr--;
  return *this;
}
inline VertexPtr VertexPtr::operator--(int) {
  ptr--;
  return VertexPtr(ptr + 1);
}

inline ::std::ostream& operator<<(::std::ostream& output, const VertexPtr& v) {
  output << "v_" << v.ptr;
  return output;
}

inline DynamicVertexPtr::DynamicVertexPtr(Vertex* ptr_, HalfedgeMesh* mesh_) : ind(mesh_->indexOf(ptr_)), mesh(mesh_) { } 
inline DynamicVertexPtr::DynamicVertexPtr(VertexPtr ptr_, HalfedgeMesh* mesh_) : ind(mesh_->indexOf(ptr_.ptr)), mesh(mesh_) { }
inline bool DynamicVertexPtr::operator==(const DynamicVertexPtr& other) const {
  return ind == other.ind;
}
inline bool DynamicVertexPtr::operator!=(const DynamicVertexPtr& other) const {
  return !(*this == other);
}
inline bool DynamicVertexPtr::operator>(const DynamicVertexPtr& other) const {
  return ind > other.ind;
}
inline bool DynamicVertexPtr::operator>=(const DynamicVertexPtr& other) const {
  return ind >= other.ind;
}
inline bool DynamicVertexPtr::operator<(const DynamicVertexPtr& other) const {
  return ind < other.ind;
}
inline bool DynamicVertexPtr::operator<=(const DynamicVertexPtr& other) const {
  return ind <= other.ind;
}

inline size_t DynamicVertexPtr::getInd() const {return ind;}

inline VertexPtrRangeIterator::VertexPtrRangeIterator(VertexPtr startingVertex)
    : currVertex(startingVertex) {}
inline const VertexPtrRangeIterator& VertexPtrRangeIterator::operator++() {
  currVertex++;
  return *this;
}
inline bool VertexPtrRangeIterator::operator==(
    const VertexPtrRangeIterator& other) const {
  return currVertex == other.currVertex;
}
inline bool VertexPtrRangeIterator::operator!=(
    const VertexPtrRangeIterator& other) const {
  return currVertex != other.currVertex;
}
inline VertexPtr VertexPtrRangeIterator::operator*() const {
  return currVertex;
}

inline VertexPtrSet::VertexPtrSet(VertexPtr beginptr_, VertexPtr endptr_)
    : beginptr(beginptr_), endptr(endptr_) {}
inline VertexPtrRangeIterator VertexPtrSet::begin() {
  return VertexPtrRangeIterator(beginptr);
}
inline VertexPtrRangeIterator VertexPtrSet::end() {
  return VertexPtrRangeIterator(endptr);
}

// Edge
inline EdgePtr::EdgePtr() : ptr(nullptr) {}
inline EdgePtr::EdgePtr(Edge* ptr_) : ptr(ptr_) {}
inline EdgePtr::EdgePtr(DynamicEdgePtr ptr_) : EdgePtr(ptr_.mesh->edge(ptr_.ind)) {}
inline HalfedgePtr EdgePtr::halfedge() const { return ptr->halfedge; }
inline bool EdgePtr::flip() { return ptr->flip(); }
inline bool EdgePtr::isBoundary() const { return ptr->isBoundary; }
inline Edge EdgePtr::operator*() { return *ptr; }
inline Edge EdgePtr::operator*() const { return *ptr; }
inline Edge* EdgePtr::operator->() { return ptr; }
inline const Edge* EdgePtr::operator->() const { return ptr; }
inline bool EdgePtr::operator==(const EdgePtr& other) const {
  return ptr == other.ptr;
}
inline bool EdgePtr::operator!=(const EdgePtr& other) const {
  return !(*this == other);
}
inline bool EdgePtr::operator>(const EdgePtr& other) const {
  return ptr > other.ptr;
}
inline bool EdgePtr::operator>=(const EdgePtr& other) const {
  return ptr >= other.ptr;
}
inline bool EdgePtr::operator<(const EdgePtr& other) const {
  return ptr < other.ptr;
}
inline bool EdgePtr::operator<=(const EdgePtr& other) const {
  return ptr <= other.ptr;
}
inline bool EdgePtr::operator==(void* n) const { return ptr == n; }
inline bool EdgePtr::operator!=(void* n) const { return ptr != n; }
inline bool EdgePtr::operator>(void* n) const { return ptr > n; }
inline bool EdgePtr::operator>=(void* n) const { return ptr >= n; }
inline bool EdgePtr::operator<(void* n) const { return ptr < n; }
inline bool EdgePtr::operator<=(void* n) const { return ptr <= n; }
inline unsigned int EdgePtr::operator-(const EdgePtr& other) const {
#ifndef NDEBUG
  assert(ptr >= other.ptr &&
         "Pointer subtraction must yield a nonnegative result");
#endif
  return ptr - other.ptr;
}
inline EdgePtr& EdgePtr::operator++() {
  ptr++;
  return *this;
}
inline EdgePtr EdgePtr::operator++(int) {
  ptr++;
  return EdgePtr(ptr - 1);
}
inline EdgePtr& EdgePtr::operator--() {
  ptr--;
  return *this;
}
inline EdgePtr EdgePtr::operator--(int) {
  ptr--;
  return EdgePtr(ptr + 1);
}

inline ::std::ostream& operator<<(::std::ostream& output, const EdgePtr& e) {
  output << "e_" << e.ptr;
  return output;
}

inline DynamicEdgePtr::DynamicEdgePtr(Edge* ptr_, HalfedgeMesh* mesh_) : ind(mesh_->indexOf(ptr_)), mesh(mesh_) { }
inline DynamicEdgePtr::DynamicEdgePtr(EdgePtr ptr_, HalfedgeMesh* mesh_) : ind(mesh_->indexOf(ptr_.ptr)), mesh(mesh_) { }
inline bool DynamicEdgePtr::operator==(const DynamicEdgePtr& other) const {
  return ind == other.ind;
}
inline bool DynamicEdgePtr::operator!=(const DynamicEdgePtr& other) const {
  return !(*this == other);
}
inline bool DynamicEdgePtr::operator>(const DynamicEdgePtr& other) const {
  return ind > other.ind;
}
inline bool DynamicEdgePtr::operator>=(const DynamicEdgePtr& other) const {
  return ind >= other.ind;
}
inline bool DynamicEdgePtr::operator<(const DynamicEdgePtr& other) const {
  return ind < other.ind;
}
inline bool DynamicEdgePtr::operator<=(const DynamicEdgePtr& other) const {
  return ind <= other.ind;
}
inline size_t DynamicEdgePtr::getInd() const {return ind;}

inline EdgePtrRangeIterator::EdgePtrRangeIterator(EdgePtr startingEdge)
    : currEdge(startingEdge) {}
inline const EdgePtrRangeIterator& EdgePtrRangeIterator::operator++() {
  currEdge++;
  return *this;
}
inline bool EdgePtrRangeIterator::operator==(
    const EdgePtrRangeIterator& other) const {
  return currEdge == other.currEdge;
}
inline bool EdgePtrRangeIterator::operator!=(
    const EdgePtrRangeIterator& other) const {
  return currEdge != other.currEdge;
}
inline EdgePtr EdgePtrRangeIterator::operator*() const { return currEdge; }

inline EdgePtrSet::EdgePtrSet(EdgePtr beginptr_, EdgePtr endptr_)
    : beginptr(beginptr_), endptr(endptr_) {}
inline EdgePtrRangeIterator EdgePtrSet::begin() {
  return EdgePtrRangeIterator(beginptr);
}
inline EdgePtrRangeIterator EdgePtrSet::end() {
  return EdgePtrRangeIterator(endptr);
}

// Face
inline FacePtr::FacePtr() : ptr(nullptr) {}
inline FacePtr::FacePtr(Face* ptr_) : ptr(ptr_) {}
inline FacePtr::FacePtr(DynamicFacePtr ptr_) : FacePtr(ptr_.mesh->face(ptr_.ind)) {}
inline HalfedgePtr FacePtr::halfedge() const { return ptr->halfedge; }
inline CornerPtr FacePtr::corner() const { return halfedge().next().corner(); }
inline unsigned int FacePtr::degree() {
  unsigned int k = 0;
  for (HalfedgePtr h : adjacentHalfedges()) {
    k++;
  }
  return k;
}
inline bool FacePtr::isBoundary() const { return ptr->isBoundary; }
inline bool FacePtr::isReal() const { return ptr->isReal; }
inline FaceAdjacentHalfedgeSet FacePtr::adjacentHalfedges() {
  return FaceAdjacentHalfedgeSet(halfedge());
}
inline FaceAdjacentVertexSet FacePtr::adjacentVertices() {
  return FaceAdjacentVertexSet(halfedge());
}
inline FaceAdjacentFaceSet FacePtr::adjacentFaces() {
  return FaceAdjacentFaceSet(halfedge());
}
inline FaceAdjacentEdgeSet FacePtr::adjacentEdges() {
  return FaceAdjacentEdgeSet(halfedge());
}
inline FaceAdjacentCornerSet FacePtr::adjacentCorners() {
  return FaceAdjacentCornerSet(halfedge());
}
inline Face FacePtr::operator*() { return *ptr; }
inline Face FacePtr::operator*() const { return *ptr; }
inline Face* FacePtr::operator->() { return ptr; }
inline const Face* FacePtr::operator->() const { return ptr; }
inline bool FacePtr::operator==(const FacePtr& other) const {
  return ptr == other.ptr;
}
inline bool FacePtr::operator!=(const FacePtr& other) const {
  return !(*this == other);
}
inline bool FacePtr::operator>(const FacePtr& other) const {
  return ptr > other.ptr;
}
inline bool FacePtr::operator>=(const FacePtr& other) const {
  return ptr >= other.ptr;
}
inline bool FacePtr::operator<(const FacePtr& other) const {
  return ptr < other.ptr;
}
inline bool FacePtr::operator<=(const FacePtr& other) const {
  return ptr <= other.ptr;
}
inline bool FacePtr::operator==(void* n) const { return ptr == n; }
inline bool FacePtr::operator!=(void* n) const { return ptr != n; }
inline bool FacePtr::operator>(void* n) const { return ptr > n; }
inline bool FacePtr::operator>=(void* n) const { return ptr >= n; }
inline bool FacePtr::operator<(void* n) const { return ptr < n; }
inline bool FacePtr::operator<=(void* n) const { return ptr <= n; }
inline unsigned int FacePtr::operator-(const FacePtr& other) const {
#ifndef NDEBUG
  assert(ptr >= other.ptr &&
         "Pointer subtraction must yield a nonnegative result");
#endif
  return ptr - other.ptr;
}
inline FacePtr& FacePtr::operator++() {
  ptr++;
  return *this;
}
inline FacePtr FacePtr::operator++(int) {
  ptr++;
  return FacePtr(ptr - 1);
}
inline FacePtr& FacePtr::operator--() {
  ptr--;
  return *this;
}
inline FacePtr FacePtr::operator--(int) {
  ptr--;
  return FacePtr(ptr + 1);
}

inline ::std::ostream& operator<<(::std::ostream& output, const FacePtr& f) {
  output << "f_" << f.ptr;
  return output;
}

inline DynamicFacePtr::DynamicFacePtr(Face* ptr_, HalfedgeMesh* mesh_) : ind(mesh_->indexOf(ptr_)), mesh(mesh_) { }
inline DynamicFacePtr::DynamicFacePtr(FacePtr ptr_, HalfedgeMesh* mesh_) : ind(mesh_->indexOf(ptr_.ptr)), mesh(mesh_) { }
inline bool DynamicFacePtr::operator==(const DynamicFacePtr& other) const {
  return ind == other.ind;
}
inline bool DynamicFacePtr::operator!=(const DynamicFacePtr& other) const {
  return !(*this == other);
}
inline bool DynamicFacePtr::operator>(const DynamicFacePtr& other) const {
  return ind > other.ind;
}
inline bool DynamicFacePtr::operator>=(const DynamicFacePtr& other) const {
  return ind >= other.ind;
}
inline bool DynamicFacePtr::operator<(const DynamicFacePtr& other) const {
  return ind < other.ind;
}
inline bool DynamicFacePtr::operator<=(const DynamicFacePtr& other) const {
  return ind <= other.ind;
}
inline size_t DynamicFacePtr::getInd() const {return ind;}

inline FacePtrRangeIterator::FacePtrRangeIterator(FacePtr startingFace)
    : currFace(startingFace) {}
inline const FacePtrRangeIterator& FacePtrRangeIterator::operator++() {
  currFace++;
  return *this;
}
inline bool FacePtrRangeIterator::operator==(
    const FacePtrRangeIterator& other) const {
  return currFace == other.currFace;
}
inline bool FacePtrRangeIterator::operator!=(
    const FacePtrRangeIterator& other) const {
  return currFace != other.currFace;
}
inline FacePtr FacePtrRangeIterator::operator*() const { return currFace; }

inline FacePtrSet::FacePtrSet(FacePtr beginptr_, FacePtr endptr_)
    : beginptr(beginptr_), endptr(endptr_) {}
inline FacePtrRangeIterator FacePtrSet::begin() {
  return FacePtrRangeIterator(beginptr);
}
inline FacePtrRangeIterator FacePtrSet::end() {
  return FacePtrRangeIterator(endptr);
}

}  // namespace geometrycentral

namespace std {

// == Hash standard pointers

template <>
struct hash<geometrycentral::HalfedgePtr> {
  std::size_t operator()(const geometrycentral::HalfedgePtr& he) const {
    return std::hash<size_t>{}(he.ptr->ID);
  }
};

template <>
struct hash<geometrycentral::VertexPtr> {
  std::size_t operator()(const geometrycentral::VertexPtr& v) const {
    return std::hash<size_t>{}(v.ptr->ID);
  }
};

template <>
struct hash<geometrycentral::EdgePtr> {
  std::size_t operator()(const geometrycentral::EdgePtr& e) const {
    return std::hash<size_t>{}(e.ptr->ID);
  }
};

template <>
struct hash<geometrycentral::FacePtr> {
  std::size_t operator()(const geometrycentral::FacePtr& f) const {
    return std::hash<size_t>{}(f.ptr->ID);
  }
};

// == Hash dynamic pointers

template <>
struct hash<geometrycentral::DynamicHalfedgePtr> {
  std::size_t operator()(const geometrycentral::DynamicHalfedgePtr& he) const {
    return std::hash<size_t>{}(he.ind);
  }
};

template <>
struct hash<geometrycentral::DynamicVertexPtr> {
  std::size_t operator()(const geometrycentral::DynamicVertexPtr& v) const {
    return std::hash<size_t>{}(v.ind);
  }
};

template <>
struct hash<geometrycentral::DynamicEdgePtr> {
  std::size_t operator()(const geometrycentral::DynamicEdgePtr& e) const {
    return std::hash<size_t>{}(e.ind);
  }
};

template <>
struct hash<geometrycentral::DynamicFacePtr> {
  std::size_t operator()(const geometrycentral::DynamicFacePtr& f) const {
    return std::hash<size_t>{}(f.ind);
  }
};
}
