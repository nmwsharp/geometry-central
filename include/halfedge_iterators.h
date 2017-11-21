#pragma once

#include <iterator>

#include <halfedge_mesh.h>

namespace geometrycentral {

// TODO add const interators across the board

// NOTE: These iterators are not STL compliant (so you can use them with
// <algorithm> and friends.
// This is mainly becuase the STL notion of iterators seems to strongly imply
// each "container"
// has exactly one set of data to be iterated over. Obviously we break this
// here, we don't even
// have containers.

// The ugly inlining throughout these iterators was chosen after some
// halfhearted performance
// testing. When first implemented, these functions seemed to be a factor of 4
// slower than
// the equivalent do{} while() loops. Now, they're are 1.0x-1.5x the cost.

// ==========================================================
// ================    Vertex Iterators    ==================
// ==========================================================

// Iterate around all incoming halfedges (both on the interior and the boundary
// of the domain)
class VertexIncomingHalfedgeIterator {
 public:
  VertexIncomingHalfedgeIterator(HalfedgePtr startingEdge, bool justStarted);
  const VertexIncomingHalfedgeIterator& operator++();
  bool operator==(const VertexIncomingHalfedgeIterator& other) const;
  bool operator!=(const VertexIncomingHalfedgeIterator& other) const;
  HalfedgePtr operator*() const;
  // HalfedgePtr operator-> () const;

 private:
  HalfedgePtr currHe;
  bool justStarted;
};
class VertexIncomingHalfedgeSet {
 public:
  VertexIncomingHalfedgeSet(HalfedgePtr he);
  VertexIncomingHalfedgeIterator begin();
  VertexIncomingHalfedgeIterator end();

 private:
  HalfedgePtr firstHe;
};

// Iterate around all incoming halfedges that are strictly on the interior of
// the domain
class VertexIncomingInteriorHalfedgeIterator {
 public:
  VertexIncomingInteriorHalfedgeIterator(HalfedgePtr startingEdge,
                                         bool justStarted);
  const VertexIncomingInteriorHalfedgeIterator& operator++();
  bool operator==(const VertexIncomingInteriorHalfedgeIterator& other) const;
  bool operator!=(const VertexIncomingInteriorHalfedgeIterator& other) const;
  HalfedgePtr operator*() const;
  // HalfedgePtr operator-> () const;

 private:
  HalfedgePtr currHe;
  bool justStarted;
};
class VertexIncomingInteriorHalfedgeSet {
 public:
  VertexIncomingInteriorHalfedgeSet(HalfedgePtr he);
  VertexIncomingInteriorHalfedgeIterator begin();
  VertexIncomingInteriorHalfedgeIterator end();

 private:
  HalfedgePtr firstHe;
};

// Iterate around all outgoing halfedges (both on the interior and the boundary
// of the domain)
class VertexOutgoingHalfedgeIterator {
 public:
  VertexOutgoingHalfedgeIterator(HalfedgePtr startingEdge, bool justStarted);
  const VertexOutgoingHalfedgeIterator& operator++();
  bool operator==(const VertexOutgoingHalfedgeIterator& other) const;
  bool operator!=(const VertexOutgoingHalfedgeIterator& other) const;
  HalfedgePtr operator*() const;
  // HalfedgePtr operator-> () const;

 private:
  HalfedgePtr currHe;
  bool justStarted;
};
class VertexOutgoingHalfedgeSet {
 public:
  VertexOutgoingHalfedgeSet(HalfedgePtr he);
  VertexOutgoingHalfedgeIterator begin();
  VertexOutgoingHalfedgeIterator end();

 private:
  HalfedgePtr firstHe;
};

// Iterate around all outgoing halfedges that are strictly on the interior of
// the domain
class VertexOutgoingInteriorHalfedgeIterator {
 public:
  VertexOutgoingInteriorHalfedgeIterator(HalfedgePtr startingEdge,
                                         bool justStarted);
  const VertexOutgoingInteriorHalfedgeIterator& operator++();
  bool operator==(const VertexOutgoingInteriorHalfedgeIterator& other) const;
  bool operator!=(const VertexOutgoingInteriorHalfedgeIterator& other) const;
  HalfedgePtr operator*() const;
  // HalfedgePtr operator-> () const;

 private:
  HalfedgePtr currHe;
  bool justStarted;
};
class VertexOutgoingInteriorHalfedgeSet {
 public:
  VertexOutgoingInteriorHalfedgeSet(HalfedgePtr he);
  VertexOutgoingInteriorHalfedgeIterator begin();
  VertexOutgoingInteriorHalfedgeIterator end();

 private:
  HalfedgePtr firstHe;
};

// Iterate around adjacent vertices
class VertexAdjacentVertexIterator {
 public:
  VertexAdjacentVertexIterator(HalfedgePtr startingEdge, bool justStarted);
  const VertexAdjacentVertexIterator& operator++();
  bool operator==(const VertexAdjacentVertexIterator& other) const;
  bool operator!=(const VertexAdjacentVertexIterator& other) const;
  VertexPtr operator*() const;
  // VertexPtr operator-> () const;

 private:
  HalfedgePtr currHe;
  bool justStarted;
};
class VertexAdjacentVertexSet {
 public:
  VertexAdjacentVertexSet(HalfedgePtr he);
  VertexAdjacentVertexIterator begin();
  VertexAdjacentVertexIterator end();

 private:
  HalfedgePtr firstHe;
};

// Iterate around adjacent (real) faces
class VertexAdjacentFaceIterator {
 public:
  VertexAdjacentFaceIterator(HalfedgePtr startingEdge, bool justStarted);
  const VertexAdjacentFaceIterator& operator++();
  bool operator==(const VertexAdjacentFaceIterator& other) const;
  bool operator!=(const VertexAdjacentFaceIterator& other) const;
  FacePtr operator*() const;
  // FacePtr operator-> () const;

 private:
  HalfedgePtr currHe;
  bool justStarted;
};
class VertexAdjacentFaceSet {
 public:
  VertexAdjacentFaceSet(HalfedgePtr he);
  VertexAdjacentFaceIterator begin();
  VertexAdjacentFaceIterator end();

 private:
  HalfedgePtr firstHe;
};

// Iterate around adjacent edges
class VertexAdjacentEdgeIterator {
 public:
  VertexAdjacentEdgeIterator(HalfedgePtr startingEdge, bool justStarted);
  const VertexAdjacentEdgeIterator& operator++();
  bool operator==(const VertexAdjacentEdgeIterator& other) const;
  bool operator!=(const VertexAdjacentEdgeIterator& other) const;
  EdgePtr operator*() const;
  // EdgePtr operator-> () const;

 private:
  HalfedgePtr currHe;
  bool justStarted;
};
class VertexAdjacentEdgeSet {
 public:
  VertexAdjacentEdgeSet(HalfedgePtr he);
  VertexAdjacentEdgeIterator begin();
  VertexAdjacentEdgeIterator end();

 private:
  HalfedgePtr firstHe;
};

// Iterate around all adjacent corners
class VertexAdjacentCornerIterator {
 public:
  VertexAdjacentCornerIterator(HalfedgePtr startingEdge, bool justStarted);
  const VertexAdjacentCornerIterator& operator++();
  bool operator==(const VertexAdjacentCornerIterator& other) const;
  bool operator!=(const VertexAdjacentCornerIterator& other) const;
  CornerPtr operator*() const;
  // HalfedgePtr operator-> () const;

 private:
  HalfedgePtr currHe;
  bool justStarted;
};
class VertexAdjacentCornerSet {
 public:
  VertexAdjacentCornerSet(HalfedgePtr he);
  VertexAdjacentCornerIterator begin();
  VertexAdjacentCornerIterator end();

 private:
  HalfedgePtr firstHe;
};

// ==========================================================
// ================     Face Iterators     ==================
// ==========================================================

// Iterate around adjacent halfedges
class FaceAdjacentHalfedgeIterator {
 public:
  FaceAdjacentHalfedgeIterator(HalfedgePtr startingEdge, bool justStarted);
  const FaceAdjacentHalfedgeIterator& operator++();
  bool operator==(const FaceAdjacentHalfedgeIterator& other) const;
  bool operator!=(const FaceAdjacentHalfedgeIterator& other) const;
  HalfedgePtr operator*() const;

 private:
  HalfedgePtr currHe;
  bool justStarted;
};
class FaceAdjacentHalfedgeSet {
 public:
  FaceAdjacentHalfedgeSet(HalfedgePtr he);
  FaceAdjacentHalfedgeIterator begin();
  FaceAdjacentHalfedgeIterator end();

 private:
  HalfedgePtr firstHe;
};

// Iterate around adjacent vertices
class FaceAdjacentVertexIterator {
 public:
  FaceAdjacentVertexIterator(HalfedgePtr startingEdge, bool justStarted);
  const FaceAdjacentVertexIterator& operator++();
  bool operator==(const FaceAdjacentVertexIterator& other) const;
  bool operator!=(const FaceAdjacentVertexIterator& other) const;
  VertexPtr operator*() const;

 private:
  HalfedgePtr currHe;
  bool justStarted;
};
class FaceAdjacentVertexSet {
 public:
  FaceAdjacentVertexSet(HalfedgePtr he);
  FaceAdjacentVertexIterator begin();
  FaceAdjacentVertexIterator end();

 private:
  HalfedgePtr firstHe;
};

// Iterate around adjacent edges
class FaceAdjacentEdgeIterator {
 public:
  FaceAdjacentEdgeIterator(HalfedgePtr startingEdge, bool justStarted);
  const FaceAdjacentEdgeIterator& operator++();
  bool operator==(const FaceAdjacentEdgeIterator& other) const;
  bool operator!=(const FaceAdjacentEdgeIterator& other) const;
  EdgePtr operator*() const;

 private:
  HalfedgePtr currHe;
  bool justStarted;
};
class FaceAdjacentEdgeSet {
 public:
  FaceAdjacentEdgeSet(HalfedgePtr he);
  FaceAdjacentEdgeIterator begin();
  FaceAdjacentEdgeIterator end();

 private:
  HalfedgePtr firstHe;
};

// Iterate around adjacent (real) faces
class FaceAdjacentFaceIterator {
 public:
  FaceAdjacentFaceIterator(HalfedgePtr startingEdge, bool justStarted);
  const FaceAdjacentFaceIterator& operator++();
  bool operator==(const FaceAdjacentFaceIterator& other) const;
  bool operator!=(const FaceAdjacentFaceIterator& other) const;
  FacePtr operator*() const;

 private:
  HalfedgePtr currHe;
  bool justStarted;
};
class FaceAdjacentFaceSet {
 public:
  FaceAdjacentFaceSet(HalfedgePtr he);
  FaceAdjacentFaceIterator begin();
  FaceAdjacentFaceIterator end();

 private:
  HalfedgePtr firstHe;
};

// Iterate around all adjacent corners
class FaceAdjacentCornerIterator {
 public:
  FaceAdjacentCornerIterator(HalfedgePtr startingEdge, bool justStarted_);
  const FaceAdjacentCornerIterator& operator++();
  bool operator==(const FaceAdjacentCornerIterator& other) const;
  bool operator!=(const FaceAdjacentCornerIterator& other) const;
  CornerPtr operator*() const;
  // HalfedgePtr operator-> () const;

 private:
  HalfedgePtr currHe;
  bool justStarted;
};
class FaceAdjacentCornerSet {
 public:
  FaceAdjacentCornerSet(HalfedgePtr he);
  FaceAdjacentCornerIterator begin();
  FaceAdjacentCornerIterator end();

 private:
  HalfedgePtr firstHe;
};

} // namespace geometrycentral