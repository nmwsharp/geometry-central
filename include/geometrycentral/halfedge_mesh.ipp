#pragma once

namespace geometrycentral {

// Methods for getting number of mesh elements

// Primal

inline size_t HalfedgeMesh::nHalfedges(void) const { return nRealHalfedges; }
inline size_t HalfedgeMesh::nCorners(void) const { return rawHalfedges.size(); }
inline size_t HalfedgeMesh::nVertices(void) const { return rawVertices.size(); }
inline size_t HalfedgeMesh::nEdges(void) const { return rawEdges.size(); }
inline size_t HalfedgeMesh::nFaces(void) const { return rawFaces.size(); }
inline size_t HalfedgeMesh::nBoundaryLoops(void) const { return rawBoundaryLoops.size(); }
inline size_t HalfedgeMesh::nImaginaryHalfedges(void) const { return rawHalfedges.size() - nRealHalfedges; }

// Methods for iterating over mesh elements w/ range-based for loops ===========

// Primal

inline HalfedgePtrSet HalfedgeMesh::realHalfedges(void) {
  HalfedgePtr beginptr(&rawHalfedges[0]);
  HalfedgePtr endptr(&rawHalfedges[nRealHalfedges]); // note that we really want nH and not
                                                     // nH-1, since we want the address one
                                                     // past the final element

  return HalfedgePtrSet(beginptr, endptr);
}

inline HalfedgePtrSet HalfedgeMesh::imaginaryHalfedges(void) {
  size_t nH = rawHalfedges.size();
  HalfedgePtr beginptr(&rawHalfedges[nRealHalfedges]);
  HalfedgePtr endptr(&rawHalfedges[nH]);

  return HalfedgePtrSet(beginptr, endptr);
}

inline HalfedgePtrSet HalfedgeMesh::allHalfedges(void) {
  size_t nH = rawHalfedges.size();
  HalfedgePtr beginptr(&rawHalfedges[0]);
  HalfedgePtr endptr(&rawHalfedges[nH]);

  return HalfedgePtrSet(beginptr, endptr);
}

inline CornerPtrSet HalfedgeMesh::corners(void) {
  size_t nC = rawHalfedges.size();
  CornerPtr beginptr(&rawHalfedges[0]);
  CornerPtr endptr(&rawHalfedges[nC]);

  return CornerPtrSet(beginptr, endptr);
}

inline VertexPtrSet HalfedgeMesh::vertices(void) {
  size_t nV = rawVertices.size();
  VertexPtr beginptr{&rawVertices[0]};
  VertexPtr endptr{&rawVertices[nV]};

  return VertexPtrSet(beginptr, endptr);
}

inline EdgePtrSet HalfedgeMesh::edges(void) {
  size_t nE = rawEdges.size();
  EdgePtr beginptr{&rawEdges[0]};
  EdgePtr endptr{&rawEdges[nE]};

  return EdgePtrSet(beginptr, endptr);
}

inline FacePtrSet HalfedgeMesh::faces(void) {
  size_t nF = rawFaces.size();
  FacePtr beginptr{&rawFaces[0]};
  FacePtr endptr{&rawFaces[nF]};

  return FacePtrSet(beginptr, endptr);
}

inline BoundaryPtrSet HalfedgeMesh::boundaryLoops(void) {
  size_t nBL = rawBoundaryLoops.size();
  BoundaryPtr beginptr{&rawBoundaryLoops[0]};
  BoundaryPtr endptr{&rawBoundaryLoops[nBL]};

  return BoundaryPtrSet(beginptr, endptr);
}

// Methods for accessing elements by index =====================================

// Primal

inline HalfedgePtr HalfedgeMesh::realHalfedge(size_t index) { return HalfedgePtr{&rawHalfedges[index]}; }

inline HalfedgePtr HalfedgeMesh::imaginaryHalfedge(size_t index) {
  return HalfedgePtr{&rawHalfedges[nRealHalfedges + index]};
}

inline HalfedgePtr HalfedgeMesh::allHalfedge(size_t index) { return HalfedgePtr{&rawHalfedges[index]}; }

inline CornerPtr HalfedgeMesh::corner(size_t index) { return CornerPtr{&rawHalfedges[index]}; }

inline VertexPtr HalfedgeMesh::vertex(size_t index) { return VertexPtr{&rawVertices[index]}; }

inline EdgePtr HalfedgeMesh::edge(size_t index) { return EdgePtr{&rawEdges[index]}; }

inline FacePtr HalfedgeMesh::face(size_t index) { return FacePtr{&rawFaces[index]}; }

inline BoundaryPtr HalfedgeMesh::boundaryLoop(size_t index) { return BoundaryPtr{&rawBoundaryLoops[index]}; }

} // namespace geometrycentral
