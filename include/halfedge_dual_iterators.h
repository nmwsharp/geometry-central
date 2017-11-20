#pragma once

#include <iterator>

#include <halfedge_mesh.h>

// ===============================================================
// ================    Dual Vertex Iterators    ==================
// ===============================================================

// TODO "Interior" iterators have not yet been implemented; some careful thought
// about what this means for dual meshes is needed.

// Iterate around all incoming halfedges (both on the interior and the boundary
// of the domain)
class DualVertexIncomingDualHalfedgeIterator {
 public:
  DualVertexIncomingDualHalfedgeIterator(DualHalfedgePtr startingDualEdge,
                                         bool justStarted);
  const DualVertexIncomingDualHalfedgeIterator& operator++();
  bool operator==(const DualVertexIncomingDualHalfedgeIterator& other) const;
  bool operator!=(const DualVertexIncomingDualHalfedgeIterator& other) const;
  DualHalfedgePtr operator*() const;

 private:
  DualHalfedgePtr currHe;
  bool justStarted;
};
class DualVertexIncomingDualHalfedgeSet {
 public:
  DualVertexIncomingDualHalfedgeSet(DualHalfedgePtr he);
  DualVertexIncomingDualHalfedgeIterator begin();
  DualVertexIncomingDualHalfedgeIterator end();

 private:
  DualHalfedgePtr firstHe;
};

// Iterate around all outgoing halfedges (both on the interior and the boundary
// of the domain)
class DualVertexOutgoingDualHalfedgeIterator {
 public:
  DualVertexOutgoingDualHalfedgeIterator(DualHalfedgePtr startingDualEdge,
                                         bool justStarted);
  const DualVertexOutgoingDualHalfedgeIterator& operator++();
  bool operator==(const DualVertexOutgoingDualHalfedgeIterator& other) const;
  bool operator!=(const DualVertexOutgoingDualHalfedgeIterator& other) const;
  DualHalfedgePtr operator*() const;

 private:
  DualHalfedgePtr currHe;
  bool justStarted;
};
class DualVertexOutgoingDualHalfedgeSet {
 public:
  DualVertexOutgoingDualHalfedgeSet(DualHalfedgePtr he);
  DualVertexOutgoingDualHalfedgeIterator begin();
  DualVertexOutgoingDualHalfedgeIterator end();

 private:
  DualHalfedgePtr firstHe;
};

// Iterate around adjacent vertices
class DualVertexAdjacentDualVertexIterator {
 public:
  DualVertexAdjacentDualVertexIterator(DualHalfedgePtr startingDualEdge,
                                       bool justStarted);
  const DualVertexAdjacentDualVertexIterator& operator++();
  bool operator==(const DualVertexAdjacentDualVertexIterator& other) const;
  bool operator!=(const DualVertexAdjacentDualVertexIterator& other) const;
  DualVertexPtr operator*() const;

 private:
  DualHalfedgePtr currHe;
  bool justStarted;
};
class DualVertexAdjacentDualVertexSet {
 public:
  DualVertexAdjacentDualVertexSet(DualHalfedgePtr he);
  DualVertexAdjacentDualVertexIterator begin();
  DualVertexAdjacentDualVertexIterator end();

 private:
  DualHalfedgePtr firstHe;
};

// Iterate around adjacent (real) faces
class DualVertexAdjacentDualFaceIterator {
 public:
  DualVertexAdjacentDualFaceIterator(DualHalfedgePtr startingDualEdge,
                                     bool justStarted);
  const DualVertexAdjacentDualFaceIterator& operator++();
  bool operator==(const DualVertexAdjacentDualFaceIterator& other) const;
  bool operator!=(const DualVertexAdjacentDualFaceIterator& other) const;
  DualFacePtr operator*() const;

 private:
  DualHalfedgePtr currHe;
  bool justStarted;
};
class DualVertexAdjacentDualFaceSet {
 public:
  DualVertexAdjacentDualFaceSet(DualHalfedgePtr he);
  DualVertexAdjacentDualFaceIterator begin();
  DualVertexAdjacentDualFaceIterator end();

 private:
  DualHalfedgePtr firstHe;
};

// Iterate around adjacent edges
class DualVertexAdjacentDualEdgeIterator {
 public:
  DualVertexAdjacentDualEdgeIterator(DualHalfedgePtr startingDualEdge,
                                     bool justStarted);
  const DualVertexAdjacentDualEdgeIterator& operator++();
  bool operator==(const DualVertexAdjacentDualEdgeIterator& other) const;
  bool operator!=(const DualVertexAdjacentDualEdgeIterator& other) const;
  DualEdgePtr operator*() const;

 private:
  DualHalfedgePtr currHe;
  bool justStarted;
};
class DualVertexAdjacentDualEdgeSet {
 public:
  DualVertexAdjacentDualEdgeSet(DualHalfedgePtr he);
  DualVertexAdjacentDualEdgeIterator begin();
  DualVertexAdjacentDualEdgeIterator end();

 private:
  DualHalfedgePtr firstHe;
};

// ===============================================================
// ================     Dual Face Iterators     ==================
// ===============================================================

// Iterate around adjacent halfedges
class DualFaceAdjacentDualHalfedgeIterator {
 public:
  DualFaceAdjacentDualHalfedgeIterator(DualHalfedgePtr startingDualEdge,
                                       bool justStarted);
  const DualFaceAdjacentDualHalfedgeIterator& operator++();
  bool operator==(const DualFaceAdjacentDualHalfedgeIterator& other) const;
  bool operator!=(const DualFaceAdjacentDualHalfedgeIterator& other) const;
  DualHalfedgePtr operator*() const;

 private:
  DualHalfedgePtr currHe;
  bool justStarted;
};
class DualFaceAdjacentDualHalfedgeSet {
 public:
  DualFaceAdjacentDualHalfedgeSet(DualHalfedgePtr he);
  DualFaceAdjacentDualHalfedgeIterator begin();
  DualFaceAdjacentDualHalfedgeIterator end();

 private:
  DualHalfedgePtr firstHe;
};

// Iterate around adjacent vertices
class DualFaceAdjacentDualVertexIterator {
 public:
  DualFaceAdjacentDualVertexIterator(DualHalfedgePtr startingDualEdge,
                                     bool justStarted);
  const DualFaceAdjacentDualVertexIterator& operator++();
  bool operator==(const DualFaceAdjacentDualVertexIterator& other) const;
  bool operator!=(const DualFaceAdjacentDualVertexIterator& other) const;
  DualVertexPtr operator*() const;

 private:
  DualHalfedgePtr currHe;
  bool justStarted;
};
class DualFaceAdjacentDualVertexSet {
 public:
  DualFaceAdjacentDualVertexSet(DualHalfedgePtr he);
  DualFaceAdjacentDualVertexIterator begin();
  DualFaceAdjacentDualVertexIterator end();

 private:
  DualHalfedgePtr firstHe;
};

// Iterate around adjacent edges
class DualFaceAdjacentDualEdgeIterator {
 public:
  DualFaceAdjacentDualEdgeIterator(DualHalfedgePtr startingDualEdge,
                                   bool justStarted);
  const DualFaceAdjacentDualEdgeIterator& operator++();
  bool operator==(const DualFaceAdjacentDualEdgeIterator& other) const;
  bool operator!=(const DualFaceAdjacentDualEdgeIterator& other) const;
  DualEdgePtr operator*() const;

 private:
  DualHalfedgePtr currHe;
  bool justStarted;
};
class DualFaceAdjacentDualEdgeSet {
 public:
  DualFaceAdjacentDualEdgeSet(DualHalfedgePtr he);
  DualFaceAdjacentDualEdgeIterator begin();
  DualFaceAdjacentDualEdgeIterator end();

 private:
  DualHalfedgePtr firstHe;
};

// Iterate around adjacent (real) faces
class DualFaceAdjacentDualFaceIterator {
 public:
  DualFaceAdjacentDualFaceIterator(DualHalfedgePtr startingDualEdge,
                                   bool justStarted);
  const DualFaceAdjacentDualFaceIterator& operator++();
  bool operator==(const DualFaceAdjacentDualFaceIterator& other) const;
  bool operator!=(const DualFaceAdjacentDualFaceIterator& other) const;
  DualFacePtr operator*() const;

 private:
  DualHalfedgePtr currHe;
  bool justStarted;
};
class DualFaceAdjacentDualFaceSet {
 public:
  DualFaceAdjacentDualFaceSet(DualHalfedgePtr he);
  DualFaceAdjacentDualFaceIterator begin();
  DualFaceAdjacentDualFaceIterator end();

 private:
  DualHalfedgePtr firstHe;
};
