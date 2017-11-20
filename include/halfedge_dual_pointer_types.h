#pragma once

#include <cstddef>

// === Types and inline methods for the halfedge mesh pointer and datatypes

// === Forward declare dual iterator set types (pointers may need to return
// these to support range-based for loops)
class DualVertexIncomingDualHalfedgeSet;
class DualVertexIncomingInteriorDualHalfedgeSet;
class DualVertexOutgoingDualHalfedgeSet;
class DualVertexOutgoingInteriorDualHalfedgeSet;
class DualVertexAdjacentDualVertexSet;
class DualVertexAdjacentDualFaceSet;
class DualVertexAdjacentDualEdgeSet;
class DualFaceAdjacentDualHalfedgeSet;
class DualFaceAdjacentDualVertexSet;
class DualFaceAdjacentDualEdgeSet;
class DualFaceAdjacentDualFaceSet;

// === Pointer types for mesh elements
class DualHalfedgePtr;
class DualVertexPtr;
class DualEdgePtr;
class DualFacePtr;

// Halfedge
class DualHalfedgePtr {
 public:
  DualHalfedgePtr();  // defaults to nullptr
  DualHalfedgePtr(Halfedge* ptr);

  // Connectivity
  DualHalfedgePtr twin() const;
  DualHalfedgePtr next() const;
  DualHalfedgePtr prev() const;
  DualVertexPtr vertex() const;
  DualEdgePtr edge() const;
  DualFacePtr face() const;

  // Duality
  HalfedgePtr dual() const;

  // Properties
  bool isReal() const;

  // Accessors
  Halfedge operator*();
  Halfedge operator*() const;
  Halfedge* operator->();
  const Halfedge* operator->() const;

  // Comparators
  bool operator==(const DualHalfedgePtr& other) const;
  bool operator!=(const DualHalfedgePtr& other) const;
  bool operator>(const DualHalfedgePtr& other) const;
  bool operator>=(const DualHalfedgePtr& other) const;
  bool operator<(const DualHalfedgePtr& other) const;
  bool operator<=(const DualHalfedgePtr& other) const;

  // Null comparators
  bool operator==(std::nullptr_t n) const;
  bool operator!=(std::nullptr_t n) const;
  bool operator>(std::nullptr_t n) const;
  bool operator>=(std::nullptr_t n) const;
  bool operator<(std::nullptr_t n) const;
  bool operator<=(std::nullptr_t n) const;

  // Arithmetic
  unsigned int operator-(const DualHalfedgePtr& other) const;
  DualHalfedgePtr& operator++();
  DualHalfedgePtr operator++(int);
  DualHalfedgePtr& operator--();
  DualHalfedgePtr operator--(int);

 protected:
  Halfedge* ptr = nullptr;

  friend class HalfedgeMesh;
};
class DualHalfedgePtrRangeIterator {
 public:
  DualHalfedgePtrRangeIterator(DualHalfedgePtr startingHalfedge);
  const DualHalfedgePtrRangeIterator& operator++();
  bool operator==(const DualHalfedgePtrRangeIterator& other) const;
  bool operator!=(const DualHalfedgePtrRangeIterator& other) const;
  DualHalfedgePtr operator*() const;

 private:
  DualHalfedgePtr currHalfedge;
};
class DualHalfedgePtrSet {
 public:
  DualHalfedgePtrSet(DualHalfedgePtr beginptr_, DualHalfedgePtr endptr_);
  DualHalfedgePtrRangeIterator begin();
  DualHalfedgePtrRangeIterator end();

 private:
  DualHalfedgePtr beginptr, endptr;
};

// Vertex
class DualVertexPtr {
 public:
  DualVertexPtr();  // defaults to nullptr
  DualVertexPtr(Face* ptr);

  // Connectivity
  DualHalfedgePtr halfedge() const;

  // Duality
  FacePtr dual() const;

  // Properties
  bool isBoundary() const;
  unsigned int degree();
  bool isReal() const;

  // Iterators
  DualVertexIncomingDualHalfedgeSet incomingHalfedges();
  DualVertexOutgoingDualHalfedgeSet outgoingHalfedges();
  DualVertexAdjacentDualVertexSet adjacentVertices();
  DualVertexAdjacentDualFaceSet adjacentFaces();
  DualVertexAdjacentDualEdgeSet adjacentEdges();

  // Accessors
  Face operator*();
  Face operator*() const;
  Face* operator->();
  const Face* operator->() const;

  // Comparators
  bool operator==(const DualVertexPtr& other) const;
  bool operator!=(const DualVertexPtr& other) const;
  bool operator>(const DualVertexPtr& other) const;
  bool operator>=(const DualVertexPtr& other) const;
  bool operator<(const DualVertexPtr& other) const;
  bool operator<=(const DualVertexPtr& other) const;

  // Null comparators
  bool operator==(std::nullptr_t n) const;
  bool operator!=(std::nullptr_t n) const;
  bool operator>(std::nullptr_t n) const;
  bool operator>=(std::nullptr_t n) const;
  bool operator<(std::nullptr_t n) const;
  bool operator<=(std::nullptr_t n) const;

  // Arithmetic
  unsigned int operator-(const DualVertexPtr& other) const;
  DualVertexPtr& operator++();
  DualVertexPtr operator++(int);
  DualVertexPtr& operator--();
  DualVertexPtr operator--(int);

 protected:
  Face* ptr = nullptr;

  friend class HalfedgeMesh;
};
class DualVertexPtrRangeIterator {
 public:
  DualVertexPtrRangeIterator(DualVertexPtr startingVertex);
  const DualVertexPtrRangeIterator& operator++();
  bool operator==(const DualVertexPtrRangeIterator& other) const;
  bool operator!=(const DualVertexPtrRangeIterator& other) const;
  DualVertexPtr operator*() const;

 private:
  DualVertexPtr currVertex;
};
class DualVertexPtrSet {
 public:
  DualVertexPtrSet(DualVertexPtr beginptr_, DualVertexPtr endptr_);
  DualVertexPtrRangeIterator begin();
  DualVertexPtrRangeIterator end();

 private:
  DualVertexPtr beginptr, endptr;
};

// Edge
class DualEdgePtr {
 public:
  DualEdgePtr();  // defaults to nullptr
  DualEdgePtr(Edge* ptr);

  // Connectivity
  DualHalfedgePtr halfedge() const;

  // Duality
  EdgePtr dual() const;

  // Properties
  bool isBoundary() const;

  // Accessors
  Edge operator*();
  Edge operator*() const;
  Edge* operator->();
  const Edge* operator->() const;

  // Comparators
  bool operator==(const DualEdgePtr& other) const;
  bool operator!=(const DualEdgePtr& other) const;
  bool operator>(const DualEdgePtr& other) const;
  bool operator>=(const DualEdgePtr& other) const;
  bool operator<(const DualEdgePtr& other) const;
  bool operator<=(const DualEdgePtr& other) const;

  // Null comparators
  bool operator==(std::nullptr_t n) const;
  bool operator!=(std::nullptr_t n) const;
  bool operator>(std::nullptr_t n) const;
  bool operator>=(std::nullptr_t n) const;
  bool operator<(std::nullptr_t n) const;
  bool operator<=(std::nullptr_t n) const;

  // Arithmetic
  unsigned int operator-(const DualEdgePtr& other) const;
  DualEdgePtr& operator++();
  DualEdgePtr operator++(int);
  DualEdgePtr& operator--();
  DualEdgePtr operator--(int);

 protected:
  Edge* ptr = nullptr;

  friend class HalfedgeMesh;
};
class DualEdgePtrRangeIterator {
 public:
  DualEdgePtrRangeIterator(DualEdgePtr startingEdge);
  const DualEdgePtrRangeIterator& operator++();
  bool operator==(const DualEdgePtrRangeIterator& other) const;
  bool operator!=(const DualEdgePtrRangeIterator& other) const;
  DualEdgePtr operator*() const;

 private:
  DualEdgePtr currEdge;
};
class DualEdgePtrSet {
 public:
  DualEdgePtrSet(DualEdgePtr beginptr_, DualEdgePtr endptr_);
  DualEdgePtrRangeIterator begin();
  DualEdgePtrRangeIterator end();

 private:
  DualEdgePtr beginptr, endptr;
};

// Face
class DualFacePtr {
 public:
  DualFacePtr();  // defaults to nullptr
  DualFacePtr(Vertex* ptr);

  // Connectivity
  DualHalfedgePtr halfedge() const;

  // Duality
  VertexPtr dual() const;

  // Properties
  unsigned int degree();
  bool isBoundary() const;

  // Iterators
  DualFaceAdjacentDualHalfedgeSet adjacentHalfedges();
  DualFaceAdjacentDualVertexSet adjacentVertices();
  DualFaceAdjacentDualFaceSet adjacentFaces();
  DualFaceAdjacentDualEdgeSet adjacentEdges();

  // Accessors
  Vertex operator*();
  Vertex operator*() const;
  Vertex* operator->();
  const Vertex* operator->() const;

  // Comparators
  bool operator==(const DualFacePtr& other) const;
  bool operator!=(const DualFacePtr& other) const;
  bool operator>(const DualFacePtr& other) const;
  bool operator>=(const DualFacePtr& other) const;
  bool operator<(const DualFacePtr& other) const;
  bool operator<=(const DualFacePtr& other) const;

  // Null comparators
  bool operator==(std::nullptr_t n) const;
  bool operator!=(std::nullptr_t n) const;
  bool operator>(std::nullptr_t n) const;
  bool operator>=(std::nullptr_t n) const;
  bool operator<(std::nullptr_t n) const;
  bool operator<=(std::nullptr_t n) const;

  // Arithmetic
  unsigned int operator-(const DualFacePtr& other) const;
  DualFacePtr& operator++();
  DualFacePtr operator++(int);
  DualFacePtr& operator--();
  DualFacePtr operator--(int);

 protected:
  Vertex* ptr = nullptr;

  friend class HalfedgeMesh;
};
class DualFacePtrRangeIterator {
 public:
  DualFacePtrRangeIterator(DualFacePtr startingFace);
  const DualFacePtrRangeIterator& operator++();
  bool operator==(const DualFacePtrRangeIterator& other) const;
  bool operator!=(const DualFacePtrRangeIterator& other) const;
  DualFacePtr operator*() const;

 private:
  DualFacePtr currFace;
};
class DualFacePtrSet {
 public:
  DualFacePtrSet(DualFacePtr beginptr_, DualFacePtr endptr_);
  DualFacePtrRangeIterator begin();
  DualFacePtrRangeIterator end();

 private:
  DualFacePtr beginptr, endptr;
};
