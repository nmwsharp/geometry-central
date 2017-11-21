#pragma once

namespace geometrycentral {

// Dual Halfedge
inline DualHalfedgePtr::DualHalfedgePtr() : ptr(nullptr) {}
inline DualHalfedgePtr::DualHalfedgePtr(Halfedge* ptr_) : ptr(ptr_) {}
inline DualHalfedgePtr DualHalfedgePtr::twin() const { return ptr->twin; }
inline DualHalfedgePtr DualHalfedgePtr::next() const { return ptr->twin->next; }
inline DualVertexPtr DualHalfedgePtr::vertex() const { return ptr->face; }
inline DualEdgePtr DualHalfedgePtr::edge() const { return ptr->edge; }
inline DualFacePtr DualHalfedgePtr::face() const { return ptr->vertex; }
inline HalfedgePtr DualHalfedgePtr::dual() const { return HalfedgePtr{ptr}; }
inline bool DualHalfedgePtr::isReal() const { return ptr->twin->isReal; }
inline Halfedge DualHalfedgePtr::operator*() { return *ptr; }
inline Halfedge DualHalfedgePtr::operator*() const { return *ptr; }
inline Halfedge* DualHalfedgePtr::operator->() { return ptr; }
inline const Halfedge* DualHalfedgePtr::operator->() const { return ptr; }
inline bool DualHalfedgePtr::operator==(const DualHalfedgePtr& other) const {
  return ptr == other.ptr;
}
inline bool DualHalfedgePtr::operator!=(const DualHalfedgePtr& other) const {
  return !(*this == other);
}
inline bool DualHalfedgePtr::operator>(const DualHalfedgePtr& other) const {
  return ptr > other.ptr;
}
inline bool DualHalfedgePtr::operator>=(const DualHalfedgePtr& other) const {
  return ptr >= other.ptr;
}
inline bool DualHalfedgePtr::operator<(const DualHalfedgePtr& other) const {
  return ptr < other.ptr;
}
inline bool DualHalfedgePtr::operator<=(const DualHalfedgePtr& other) const {
  return ptr <= other.ptr;
}
inline bool DualHalfedgePtr::operator==(std::nullptr_t n) const {
  return ptr == n;
}
inline bool DualHalfedgePtr::operator!=(std::nullptr_t n) const {
  return ptr != n;
}
inline bool DualHalfedgePtr::operator>(std::nullptr_t n) const {
  return ptr > n;
}
inline bool DualHalfedgePtr::operator>=(std::nullptr_t n) const {
  return ptr >= n;
}
inline bool DualHalfedgePtr::operator<(std::nullptr_t n) const {
  return ptr < n;
}
inline bool DualHalfedgePtr::operator<=(std::nullptr_t n) const {
  return ptr <= n;
}
inline unsigned int DualHalfedgePtr::operator-(
    const DualHalfedgePtr& other) const {
#ifndef NDEBUG
  assert(ptr >= other.ptr &&
         "Pointer subtraction must yield a nonnegative result");
#endif
  return ptr - other.ptr;
}
inline DualHalfedgePtr& DualHalfedgePtr::operator++() {
  ptr++;
  return *this;
}
inline DualHalfedgePtr DualHalfedgePtr::operator++(int) {
  ptr++;
  return DualHalfedgePtr(ptr - 1);
}
inline DualHalfedgePtr& DualHalfedgePtr::operator--() {
  ptr--;
  return *this;
}
inline DualHalfedgePtr DualHalfedgePtr::operator--(int) {
  ptr--;
  return DualHalfedgePtr(ptr + 1);
}
inline DualHalfedgePtr DualHalfedgePtr::prev() const {
  Halfedge* h = ptr;
  while (h->twin->next != ptr) {
    h = h->twin->next;
  }
  return DualHalfedgePtr{h};
}

inline DualHalfedgePtrRangeIterator::DualHalfedgePtrRangeIterator(
    DualHalfedgePtr startingHalfedge)
    : currHalfedge(startingHalfedge) {}
inline const DualHalfedgePtrRangeIterator& DualHalfedgePtrRangeIterator::
operator++() {
  currHalfedge++;
  return *this;
}
inline bool DualHalfedgePtrRangeIterator::operator==(
    const DualHalfedgePtrRangeIterator& other) const {
  return currHalfedge == other.currHalfedge;
}
inline bool DualHalfedgePtrRangeIterator::operator!=(
    const DualHalfedgePtrRangeIterator& other) const {
  return currHalfedge != other.currHalfedge;
}
inline DualHalfedgePtr DualHalfedgePtrRangeIterator::operator*() const {
  return currHalfedge;
}

inline DualHalfedgePtrSet::DualHalfedgePtrSet(DualHalfedgePtr beginptr_,
                                              DualHalfedgePtr endptr_)
    : beginptr(beginptr_), endptr(endptr_) {}
inline DualHalfedgePtrRangeIterator DualHalfedgePtrSet::begin() {
  return DualHalfedgePtrRangeIterator(beginptr);
}
inline DualHalfedgePtrRangeIterator DualHalfedgePtrSet::end() {
  return DualHalfedgePtrRangeIterator(endptr);
}

// Dual Vertex
inline DualVertexPtr::DualVertexPtr() : ptr(nullptr) {}
inline DualVertexPtr::DualVertexPtr(Face* ptr_) : ptr(ptr_) {}
inline DualHalfedgePtr DualVertexPtr::halfedge() const { return ptr->halfedge; }
inline FacePtr DualVertexPtr::dual() const { return FacePtr{ptr}; }
inline bool DualVertexPtr::isBoundary() const { return ptr->isBoundary; }
inline bool DualVertexPtr::isReal() const { return ptr->isReal; }
inline DualVertexIncomingDualHalfedgeSet DualVertexPtr::incomingHalfedges() {
  return DualVertexIncomingDualHalfedgeSet(halfedge().twin());
}
inline DualVertexOutgoingDualHalfedgeSet DualVertexPtr::outgoingHalfedges() {
  return DualVertexOutgoingDualHalfedgeSet(halfedge());
}
inline DualVertexAdjacentDualVertexSet DualVertexPtr::adjacentVertices() {
  return DualVertexAdjacentDualVertexSet(halfedge().twin());
}
inline DualVertexAdjacentDualFaceSet DualVertexPtr::adjacentFaces() {
  return DualVertexAdjacentDualFaceSet(halfedge());
}
inline DualVertexAdjacentDualEdgeSet DualVertexPtr::adjacentEdges() {
  return DualVertexAdjacentDualEdgeSet(halfedge());
}
inline Face DualVertexPtr::operator*() { return *ptr; }
inline Face DualVertexPtr::operator*() const { return *ptr; }
inline Face* DualVertexPtr::operator->() { return ptr; }
inline const Face* DualVertexPtr::operator->() const { return ptr; }
inline bool DualVertexPtr::operator==(const DualVertexPtr& other) const {
  return ptr == other.ptr;
}
inline bool DualVertexPtr::operator!=(const DualVertexPtr& other) const {
  return !(*this == other);
}
inline bool DualVertexPtr::operator>(const DualVertexPtr& other) const {
  return ptr > other.ptr;
}
inline bool DualVertexPtr::operator>=(const DualVertexPtr& other) const {
  return ptr >= other.ptr;
}
inline bool DualVertexPtr::operator<(const DualVertexPtr& other) const {
  return ptr < other.ptr;
}
inline bool DualVertexPtr::operator<=(const DualVertexPtr& other) const {
  return ptr <= other.ptr;
}
inline bool DualVertexPtr::operator==(std::nullptr_t n) const {
  return ptr == n;
}
inline bool DualVertexPtr::operator!=(std::nullptr_t n) const {
  return ptr != n;
}
inline bool DualVertexPtr::operator>(std::nullptr_t n) const { return ptr > n; }
inline bool DualVertexPtr::operator>=(std::nullptr_t n) const {
  return ptr >= n;
}
inline bool DualVertexPtr::operator<(std::nullptr_t n) const { return ptr < n; }
inline bool DualVertexPtr::operator<=(std::nullptr_t n) const {
  return ptr <= n;
}
inline unsigned int DualVertexPtr::operator-(const DualVertexPtr& other) const {
#ifndef NDEBUG
  assert(ptr >= other.ptr &&
         "Pointer subtraction must yield a nonnegative result");
#endif
  return ptr - other.ptr;
}
inline DualVertexPtr& DualVertexPtr::operator++() {
  ptr++;
  return *this;
}
inline DualVertexPtr DualVertexPtr::operator++(int) {
  ptr++;
  return DualVertexPtr(ptr - 1);
}
inline DualVertexPtr& DualVertexPtr::operator--() {
  ptr--;
  return *this;
}
inline DualVertexPtr DualVertexPtr::operator--(int) {
  ptr--;
  return DualVertexPtr(ptr + 1);
}

inline DualVertexPtrRangeIterator::DualVertexPtrRangeIterator(
    DualVertexPtr startingVertex)
    : currVertex(startingVertex) {}
inline const DualVertexPtrRangeIterator& DualVertexPtrRangeIterator::
operator++() {
  currVertex++;
  return *this;
}
inline bool DualVertexPtrRangeIterator::operator==(
    const DualVertexPtrRangeIterator& other) const {
  return currVertex == other.currVertex;
}
inline bool DualVertexPtrRangeIterator::operator!=(
    const DualVertexPtrRangeIterator& other) const {
  return currVertex != other.currVertex;
}
inline DualVertexPtr DualVertexPtrRangeIterator::operator*() const {
  return currVertex;
}

inline DualVertexPtrSet::DualVertexPtrSet(DualVertexPtr beginptr_,
                                          DualVertexPtr endptr_)
    : beginptr(beginptr_), endptr(endptr_) {}
inline DualVertexPtrRangeIterator DualVertexPtrSet::begin() {
  return DualVertexPtrRangeIterator(beginptr);
}
inline DualVertexPtrRangeIterator DualVertexPtrSet::end() {
  return DualVertexPtrRangeIterator(endptr);
}

// Dual Edge
inline DualEdgePtr::DualEdgePtr() : ptr(nullptr) {}
inline DualEdgePtr::DualEdgePtr(Edge* ptr_) : ptr(ptr_) {}
inline DualHalfedgePtr DualEdgePtr::halfedge() const { return ptr->halfedge; }
inline EdgePtr DualEdgePtr::dual() const { return EdgePtr{ptr}; }
inline bool DualEdgePtr::isBoundary() const { return ptr->isBoundary; }
inline Edge DualEdgePtr::operator*() { return *ptr; }
inline Edge DualEdgePtr::operator*() const { return *ptr; }
inline Edge* DualEdgePtr::operator->() { return ptr; }
inline const Edge* DualEdgePtr::operator->() const { return ptr; }
inline bool DualEdgePtr::operator==(const DualEdgePtr& other) const {
  return ptr == other.ptr;
}
inline bool DualEdgePtr::operator!=(const DualEdgePtr& other) const {
  return !(*this == other);
}
inline bool DualEdgePtr::operator>(const DualEdgePtr& other) const {
  return ptr > other.ptr;
}
inline bool DualEdgePtr::operator>=(const DualEdgePtr& other) const {
  return ptr >= other.ptr;
}
inline bool DualEdgePtr::operator<(const DualEdgePtr& other) const {
  return ptr < other.ptr;
}
inline bool DualEdgePtr::operator<=(const DualEdgePtr& other) const {
  return ptr <= other.ptr;
}
inline bool DualEdgePtr::operator==(std::nullptr_t n) const { return ptr == n; }
inline bool DualEdgePtr::operator!=(std::nullptr_t n) const { return ptr != n; }
inline bool DualEdgePtr::operator>(std::nullptr_t n) const { return ptr > n; }
inline bool DualEdgePtr::operator>=(std::nullptr_t n) const { return ptr >= n; }
inline bool DualEdgePtr::operator<(std::nullptr_t n) const { return ptr < n; }
inline bool DualEdgePtr::operator<=(std::nullptr_t n) const { return ptr <= n; }
inline unsigned int DualEdgePtr::operator-(const DualEdgePtr& other) const {
#ifndef NDEBUG
  assert(ptr >= other.ptr &&
         "Pointer subtraction must yield a nonnegative result");
#endif
  return ptr - other.ptr;
}
inline DualEdgePtr& DualEdgePtr::operator++() {
  ptr++;
  return *this;
}
inline DualEdgePtr DualEdgePtr::operator++(int) {
  ptr++;
  return DualEdgePtr(ptr - 1);
}
inline DualEdgePtr& DualEdgePtr::operator--() {
  ptr--;
  return *this;
}
inline DualEdgePtr DualEdgePtr::operator--(int) {
  ptr--;
  return DualEdgePtr(ptr + 1);
}

inline DualEdgePtrRangeIterator::DualEdgePtrRangeIterator(
    DualEdgePtr startingEdge)
    : currEdge(startingEdge) {}
inline const DualEdgePtrRangeIterator& DualEdgePtrRangeIterator::operator++() {
  currEdge++;
  return *this;
}
inline bool DualEdgePtrRangeIterator::operator==(
    const DualEdgePtrRangeIterator& other) const {
  return currEdge == other.currEdge;
}
inline bool DualEdgePtrRangeIterator::operator!=(
    const DualEdgePtrRangeIterator& other) const {
  return currEdge != other.currEdge;
}
inline DualEdgePtr DualEdgePtrRangeIterator::operator*() const {
  return currEdge;
}

inline DualEdgePtrSet::DualEdgePtrSet(DualEdgePtr beginptr_,
                                      DualEdgePtr endptr_)
    : beginptr(beginptr_), endptr(endptr_) {}
inline DualEdgePtrRangeIterator DualEdgePtrSet::begin() {
  return DualEdgePtrRangeIterator(beginptr);
}
inline DualEdgePtrRangeIterator DualEdgePtrSet::end() {
  return DualEdgePtrRangeIterator(endptr);
}

// Dual Face
inline DualFacePtr::DualFacePtr() : ptr(nullptr) {}
inline DualFacePtr::DualFacePtr(Vertex* ptr_) : ptr(ptr_) {}
inline DualHalfedgePtr DualFacePtr::halfedge() const { return ptr->halfedge; }
inline VertexPtr DualFacePtr::dual() const { return VertexPtr{ptr}; }
inline bool DualFacePtr::isBoundary() const { return ptr->isBoundary; }
inline DualFaceAdjacentDualHalfedgeSet DualFacePtr::adjacentHalfedges() {
  return DualFaceAdjacentDualHalfedgeSet(halfedge());
}
inline DualFaceAdjacentDualVertexSet DualFacePtr::adjacentVertices() {
  return DualFaceAdjacentDualVertexSet(halfedge());
}
inline DualFaceAdjacentDualFaceSet DualFacePtr::adjacentFaces() {
  return DualFaceAdjacentDualFaceSet(halfedge());
}
inline DualFaceAdjacentDualEdgeSet DualFacePtr::adjacentEdges() {
  return DualFaceAdjacentDualEdgeSet(halfedge());
}
inline Vertex DualFacePtr::operator*() { return *ptr; }
inline Vertex DualFacePtr::operator*() const { return *ptr; }
inline Vertex* DualFacePtr::operator->() { return ptr; }
inline const Vertex* DualFacePtr::operator->() const { return ptr; }
inline bool DualFacePtr::operator==(const DualFacePtr& other) const {
  return ptr == other.ptr;
}
inline bool DualFacePtr::operator!=(const DualFacePtr& other) const {
  return !(*this == other);
}
inline bool DualFacePtr::operator>(const DualFacePtr& other) const {
  return ptr > other.ptr;
}
inline bool DualFacePtr::operator>=(const DualFacePtr& other) const {
  return ptr >= other.ptr;
}
inline bool DualFacePtr::operator<(const DualFacePtr& other) const {
  return ptr < other.ptr;
}
inline bool DualFacePtr::operator<=(const DualFacePtr& other) const {
  return ptr <= other.ptr;
}
inline bool DualFacePtr::operator==(std::nullptr_t n) const { return ptr == n; }
inline bool DualFacePtr::operator!=(std::nullptr_t n) const { return ptr != n; }
inline bool DualFacePtr::operator>(std::nullptr_t n) const { return ptr > n; }
inline bool DualFacePtr::operator>=(std::nullptr_t n) const { return ptr >= n; }
inline bool DualFacePtr::operator<(std::nullptr_t n) const { return ptr < n; }
inline bool DualFacePtr::operator<=(std::nullptr_t n) const { return ptr <= n; }
inline unsigned int DualFacePtr::operator-(const DualFacePtr& other) const {
#ifndef NDEBUG
  assert(ptr >= other.ptr &&
         "Pointer subtraction must yield a nonnegative result");
#endif
  return ptr - other.ptr;
}
inline DualFacePtr& DualFacePtr::operator++() {
  ptr++;
  return *this;
}
inline DualFacePtr DualFacePtr::operator++(int) {
  ptr++;
  return DualFacePtr(ptr - 1);
}
inline DualFacePtr& DualFacePtr::operator--() {
  ptr--;
  return *this;
}
inline DualFacePtr DualFacePtr::operator--(int) {
  ptr--;
  return DualFacePtr(ptr + 1);
}
inline unsigned int DualFacePtr::degree() {
  unsigned int k = 0;
  for (DualHalfedgePtr h : adjacentHalfedges()) {
    k++;
  }
  return k;
}

inline DualFacePtrRangeIterator::DualFacePtrRangeIterator(
    DualFacePtr startingFace)
    : currFace(startingFace) {}
inline const DualFacePtrRangeIterator& DualFacePtrRangeIterator::operator++() {
  currFace++;
  return *this;
}
inline bool DualFacePtrRangeIterator::operator==(
    const DualFacePtrRangeIterator& other) const {
  return currFace == other.currFace;
}
inline bool DualFacePtrRangeIterator::operator!=(
    const DualFacePtrRangeIterator& other) const {
  return currFace != other.currFace;
}
inline DualFacePtr DualFacePtrRangeIterator::operator*() const {
  return currFace;
}

inline DualFacePtrSet::DualFacePtrSet(DualFacePtr beginptr_,
                                      DualFacePtr endptr_)
    : beginptr(beginptr_), endptr(endptr_) {}
inline DualFacePtrRangeIterator DualFacePtrSet::begin() {
  return DualFacePtrRangeIterator(beginptr);
}
inline DualFacePtrRangeIterator DualFacePtrSet::end() {
  return DualFacePtrRangeIterator(endptr);
}

} // namespace geometrycentral