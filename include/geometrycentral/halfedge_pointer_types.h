#pragma once

#include <cstddef>
#include <functional>
#include <iostream>

namespace geometrycentral {

// === Types and inline methods for the halfedge mesh pointer and datatypes

// === Forward declare iterator set types (pointers may need to return these to
// support range-based for loops)
class VertexIncomingHalfedgeSet;
class VertexIncomingInteriorHalfedgeSet;
class VertexOutgoingHalfedgeSet;
class VertexOutgoingInteriorHalfedgeSet;
class VertexAdjacentVertexSet;
class VertexAdjacentFaceSet;
class VertexAdjacentEdgeSet;
class VertexAdjacentCornerSet;
class FaceAdjacentHalfedgeSet;
class FaceAdjacentVertexSet;
class FaceAdjacentEdgeSet;
class FaceAdjacentFaceSet;
class FaceAdjacentCornerSet;

// === Pointer types for mesh elements
class HalfedgePtr;
class CornerPtr;
class VertexPtr;
class EdgePtr;
class FacePtr;
class BoundaryLoopPtr;
class DualHalfedgePtr;
class DualVertexPtr;
class DualEdgePtr;
class DualFacePtr;

// Halfedge
class HalfedgePtr {
 public:
  HalfedgePtr();  // defaults to nullptr
  HalfedgePtr(Halfedge* ptr);

  // Connectivity
  HalfedgePtr twin() const;
  HalfedgePtr next() const;
  HalfedgePtr prev() const;
  VertexPtr vertex() const;
  EdgePtr edge() const;
  FacePtr face() const;
  CornerPtr corner() const;

  // Duality
  DualHalfedgePtr dual() const;

  // Properties
  bool isReal() const;

  // Accessors
  Halfedge operator*();
  Halfedge operator*() const;
  Halfedge* operator->();
  const Halfedge* operator->() const;

  // Comparators
  bool operator==(const HalfedgePtr& other) const;
  bool operator!=(const HalfedgePtr& other) const;
  bool operator>(const HalfedgePtr& other) const;
  bool operator>=(const HalfedgePtr& other) const;
  bool operator<(const HalfedgePtr& other) const;
  bool operator<=(const HalfedgePtr& other) const;

  // Null comparators
  bool operator==(std::nullptr_t n) const;
  bool operator!=(std::nullptr_t n) const;
  bool operator>(std::nullptr_t n) const;
  bool operator>=(std::nullptr_t n) const;
  bool operator<(std::nullptr_t n) const;
  bool operator<=(std::nullptr_t n) const;

  // Arithmetic
  unsigned int operator-(const HalfedgePtr& other) const;
  HalfedgePtr& operator++();
  HalfedgePtr operator++(int);
  HalfedgePtr& operator--();
  HalfedgePtr operator--(int);

 protected:
  Halfedge* ptr = nullptr;

  friend class HalfedgeMesh;
  friend std::ostream& operator<<(std::ostream& output, const HalfedgePtr& he);
  friend struct std::hash<HalfedgePtr>;
};
std::ostream& operator<<(std::ostream& output, const HalfedgePtr& he);

class HalfedgePtrRangeIterator {
 public:
  HalfedgePtrRangeIterator(HalfedgePtr startingHalfedge);
  const HalfedgePtrRangeIterator& operator++();
  bool operator==(const HalfedgePtrRangeIterator& other) const;
  bool operator!=(const HalfedgePtrRangeIterator& other) const;
  HalfedgePtr operator*() const;

 private:
  HalfedgePtr currHalfedge;
};
class HalfedgePtrSet {
 public:
  HalfedgePtrSet(HalfedgePtr beginptr_, HalfedgePtr endptr_);
  HalfedgePtrRangeIterator begin();
  HalfedgePtrRangeIterator end();

 private:
  HalfedgePtr beginptr, endptr;
};

class CutPtrRangeIterator {
 public:
  CutPtrRangeIterator(HalfedgePtr startingHalfedge, bool justStarted_);
  const CutPtrRangeIterator& operator++();
  bool operator==(const CutPtrRangeIterator& other) const;
  bool operator!=(const CutPtrRangeIterator& other) const;
  HalfedgePtr operator*() const;

 private:
  HalfedgePtr currHalfedge;
  bool justStarted;
};
class CutPtrSet {
 public:
  CutPtrSet(HalfedgePtr he);
  CutPtrRangeIterator begin();
  CutPtrRangeIterator end();

 private:
  HalfedgePtr firstHe;
};

// Corner
class CornerPtr {
 public:
  CornerPtr();  // defaults to nullptr
  CornerPtr(Halfedge* ptr);

  // Connectivity
  CornerPtr next() const;
  CornerPtr prev() const;
  HalfedgePtr halfedge() const;
  VertexPtr vertex() const;
  FacePtr face() const;

  // Accessors
  Halfedge operator*();
  Halfedge operator*() const;
  Halfedge* operator->();
  const Halfedge* operator->() const;

  // Comparators
  bool operator==(const CornerPtr& other) const;
  bool operator!=(const CornerPtr& other) const;
  bool operator>(const CornerPtr& other) const;
  bool operator>=(const CornerPtr& other) const;
  bool operator<(const CornerPtr& other) const;
  bool operator<=(const CornerPtr& other) const;

  // Null comparators
  bool operator==(std::nullptr_t n) const;
  bool operator!=(std::nullptr_t n) const;
  bool operator>(std::nullptr_t n) const;
  bool operator>=(std::nullptr_t n) const;
  bool operator<(std::nullptr_t n) const;
  bool operator<=(std::nullptr_t n) const;

  // Arithmetic
  unsigned int operator-(const CornerPtr& other) const;
  CornerPtr& operator++();
  CornerPtr operator++(int);
  CornerPtr& operator--();
  CornerPtr operator--(int);

 protected:
  Halfedge* ptr = nullptr;

  friend class HalfedgeMesh;
};
class CornerPtrRangeIterator {
 public:
  CornerPtrRangeIterator(CornerPtr startingCorner);
  const CornerPtrRangeIterator& operator++();
  bool operator==(const CornerPtrRangeIterator& other) const;
  bool operator!=(const CornerPtrRangeIterator& other) const;
  CornerPtr operator*() const;

 private:
  CornerPtr currCorner;
};
class CornerPtrSet {
 public:
  CornerPtrSet(CornerPtr beginptr_, CornerPtr endptr_);
  CornerPtrRangeIterator begin();
  CornerPtrRangeIterator end();

 private:
  CornerPtr beginptr, endptr;
};

// Vertex
class VertexPtr {
 public:
  VertexPtr();  // defaults to nullptr
  VertexPtr(Vertex* ptr);

  // Connectivity
  HalfedgePtr halfedge() const;
  CornerPtr corner() const;

  // Duality
  DualFacePtr dual() const;

  // Properties
  bool isBoundary() const;
  unsigned int degree();

  // Iterators
  VertexIncomingHalfedgeSet incomingHalfedges();
  VertexOutgoingHalfedgeSet outgoingHalfedges();
  VertexIncomingInteriorHalfedgeSet incomingInteriorHalfedges();
  VertexOutgoingInteriorHalfedgeSet outgoingInteriorHalfedges();
  VertexAdjacentVertexSet adjacentVertices();
  VertexAdjacentFaceSet adjacentFaces();
  VertexAdjacentEdgeSet adjacentEdges();
  VertexAdjacentCornerSet adjacentCorners();

  // Accessors
  Vertex operator*();
  Vertex operator*() const;
  Vertex* operator->();
  const Vertex* operator->() const;

  // Comparators
  bool operator==(const VertexPtr& other) const;
  bool operator!=(const VertexPtr& other) const;
  bool operator>(const VertexPtr& other) const;
  bool operator>=(const VertexPtr& other) const;
  bool operator<(const VertexPtr& other) const;
  bool operator<=(const VertexPtr& other) const;

  // Null comparators
  bool operator==(std::nullptr_t n) const;
  bool operator!=(std::nullptr_t n) const;
  bool operator>(std::nullptr_t n) const;
  bool operator>=(std::nullptr_t n) const;
  bool operator<(std::nullptr_t n) const;
  bool operator<=(std::nullptr_t n) const;

  // Arithmetic
  unsigned int operator-(const VertexPtr& other) const;
  VertexPtr& operator++();
  VertexPtr operator++(int);
  VertexPtr& operator--();
  VertexPtr operator--(int);

 protected:
  Vertex* ptr = nullptr;

  friend class HalfedgeMesh;
  friend std::ostream& operator<<(std::ostream& output, const VertexPtr& v);
  friend struct std::hash<VertexPtr>;
};
std::ostream& operator<<(std::ostream& output, const VertexPtr& v);

class VertexPtrRangeIterator {
 public:
  VertexPtrRangeIterator(VertexPtr startingVertex);
  const VertexPtrRangeIterator& operator++();
  bool operator==(const VertexPtrRangeIterator& other) const;
  bool operator!=(const VertexPtrRangeIterator& other) const;
  VertexPtr operator*() const;

 private:
  VertexPtr currVertex;
};
class VertexPtrSet {
 public:
  VertexPtrSet(VertexPtr beginptr_, VertexPtr endptr_);
  VertexPtrRangeIterator begin();
  VertexPtrRangeIterator end();

 private:
  VertexPtr beginptr, endptr;
};

// Edge
class EdgePtr {
 public:
  EdgePtr();  // defaults to nullptr
  EdgePtr(Edge* ptr);

  // Connectivity
  HalfedgePtr halfedge() const;

  // Duality
  DualEdgePtr dual() const;

  // Properties
  bool isBoundary() const;
  bool isCut() const;
  void markCut(bool isCut);

  // Remeshing
  bool flip();  // flips a triangle pair, returning false if the faces aren't
                // triangles or the edge is on the boundary

  // Accessors
  Edge operator*();
  Edge operator*() const;
  Edge* operator->();
  const Edge* operator->() const;

  // Comparators
  bool operator==(const EdgePtr& other) const;
  bool operator!=(const EdgePtr& other) const;
  bool operator>(const EdgePtr& other) const;
  bool operator>=(const EdgePtr& other) const;
  bool operator<(const EdgePtr& other) const;
  bool operator<=(const EdgePtr& other) const;

  // Null comparators
  bool operator==(std::nullptr_t n) const;
  bool operator!=(std::nullptr_t n) const;
  bool operator>(std::nullptr_t n) const;
  bool operator>=(std::nullptr_t n) const;
  bool operator<(std::nullptr_t n) const;
  bool operator<=(std::nullptr_t n) const;

  // Arithmetic
  unsigned int operator-(const EdgePtr& other) const;
  EdgePtr& operator++();
  EdgePtr operator++(int);
  EdgePtr& operator--();
  EdgePtr operator--(int);

 protected:
  Edge* ptr = nullptr;

  friend class HalfedgeMesh;
  friend std::ostream& operator<<(std::ostream& output, const EdgePtr& e);
  friend struct std::hash<EdgePtr>;
};
std::ostream& operator<<(std::ostream& output, const EdgePtr& e);

class EdgePtrRangeIterator {
 public:
  EdgePtrRangeIterator(EdgePtr startingEdge);
  const EdgePtrRangeIterator& operator++();
  bool operator==(const EdgePtrRangeIterator& other) const;
  bool operator!=(const EdgePtrRangeIterator& other) const;
  EdgePtr operator*() const;

 private:
  EdgePtr currEdge;
};
class EdgePtrSet {
 public:
  EdgePtrSet(EdgePtr beginptr_, EdgePtr endptr_);
  EdgePtrRangeIterator begin();
  EdgePtrRangeIterator end();

 private:
  EdgePtr beginptr, endptr;
};

// NOTE: Triangle is merely a helper class for the
// method FacePtr::triangulation(); it is NOT the
// standard representation for faces of a HalfedgeMesh.
struct Triangle {
  VertexPtr vertex[3];
  VertexPtr& operator[](size_t index) { return vertex[index]; }
};

// Face
class FacePtr {
 public:
  FacePtr();  // defaults to nullptr
  FacePtr(Face* ptr);

  // Connectivity
  HalfedgePtr halfedge() const;
  CornerPtr corner() const;

  // Duality
  DualVertexPtr dual() const;

  // Utility
  std::vector<Triangle> triangulation();

  // Properties
  unsigned int degree();
  bool isBoundary() const;
  bool isReal() const;

  // Iterators
  FaceAdjacentHalfedgeSet adjacentHalfedges();
  FaceAdjacentVertexSet adjacentVertices();
  FaceAdjacentFaceSet adjacentFaces();
  FaceAdjacentEdgeSet adjacentEdges();
  FaceAdjacentCornerSet adjacentCorners();

  // Accessors
  Face operator*();
  Face operator*() const;
  Face* operator->();
  const Face* operator->() const;

  // Comparators
  bool operator==(const FacePtr& other) const;
  bool operator!=(const FacePtr& other) const;
  bool operator>(const FacePtr& other) const;
  bool operator>=(const FacePtr& other) const;
  bool operator<(const FacePtr& other) const;
  bool operator<=(const FacePtr& other) const;

  // Null comparators
  bool operator==(std::nullptr_t n) const;
  bool operator!=(std::nullptr_t n) const;
  bool operator>(std::nullptr_t n) const;
  bool operator>=(std::nullptr_t n) const;
  bool operator<(std::nullptr_t n) const;
  bool operator<=(std::nullptr_t n) const;

  // Arithmetic
  unsigned int operator-(const FacePtr& other) const;
  FacePtr& operator++();
  FacePtr operator++(int);
  FacePtr& operator--();
  FacePtr operator--(int);

 protected:
  Face* ptr = nullptr;

  friend class HalfedgeMesh;
  friend std::ostream& operator<<(std::ostream& output, const FacePtr& f);
  friend struct std::hash<FacePtr>;
};
std::ostream& operator<<(std::ostream& output, const FacePtr& f);

class FacePtrRangeIterator {
 public:
  FacePtrRangeIterator(FacePtr startingFace);
  const FacePtrRangeIterator& operator++();
  bool operator==(const FacePtrRangeIterator& other) const;
  bool operator!=(const FacePtrRangeIterator& other) const;
  FacePtr operator*() const;

 private:
  FacePtr currFace;
};
class FacePtrSet {
 public:
  FacePtrSet(FacePtr beginptr_, FacePtr endptr_);
  FacePtrRangeIterator begin();
  FacePtrRangeIterator end();

 private:
  FacePtr beginptr, endptr;
};

// Boundary (currently just a renaming of Face---if we wanted
// stronger type checking we could instead inherit from Face*)

typedef Face Boundary;
typedef FacePtr BoundaryPtr;
typedef FacePtrSet BoundaryPtrSet;
typedef FacePtrRangeIterator BoundaryRangeIterator;

}  // namespace geometrycentral