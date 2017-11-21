#include <halfedge_mesh.h>

// Implementation of various methods to support iterators in the halfedge mesh.
// Note that many functions are inlined, and just defined in
// halfedge_iterators.ipp

namespace geometrycentral {

// ==========================================================
// ================    Vertex Iterators    ==================
// ==========================================================

// === Incoming halfedges (with imaginary)
VertexIncomingHalfedgeSet::VertexIncomingHalfedgeSet(HalfedgePtr he)
    : firstHe(he) {}

VertexIncomingHalfedgeIterator::VertexIncomingHalfedgeIterator(
    HalfedgePtr startingEdge, bool justStarted_)
    : currHe(startingEdge), justStarted(justStarted_) {}

// === Incoming halfedges (without imaginary)
VertexIncomingInteriorHalfedgeSet::VertexIncomingInteriorHalfedgeSet(
    HalfedgePtr he)
    : firstHe(he) {}

VertexIncomingInteriorHalfedgeIterator::VertexIncomingInteriorHalfedgeIterator(
    HalfedgePtr startingEdge, bool justStarted_)
    : currHe(startingEdge), justStarted(justStarted_) {
  // If startingEdge is a boundary edge, the first state of this iterator would
  // be a halfedge which should not be returned,
  // so advance the iterator until we find valid one.
  // TODO initialize the mesh so that the starting edge satisfies this property
  // by default

  while (!currHe.isReal()) {
    currHe = currHe.next().twin();
  }
}

// === Outgoing halfedges (with imaginary)
VertexOutgoingHalfedgeSet::VertexOutgoingHalfedgeSet(HalfedgePtr he)
    : firstHe(he) {}

VertexOutgoingHalfedgeIterator::VertexOutgoingHalfedgeIterator(
    HalfedgePtr startingEdge, bool justStarted_)
    : currHe(startingEdge), justStarted(justStarted_) {}

// === Outgoing halfedges (without imaginary)
VertexOutgoingInteriorHalfedgeSet::VertexOutgoingInteriorHalfedgeSet(
    HalfedgePtr he)
    : firstHe(he) {}

VertexOutgoingInteriorHalfedgeIterator::VertexOutgoingInteriorHalfedgeIterator(
    HalfedgePtr startingEdge, bool justStarted_)
    : currHe(startingEdge), justStarted(justStarted_) {
  // If startingEdge is a boundary edge, the first state of this iterator would
  // be a halfedge which should not be returned,
  // so advance the iterator until we find valid one.
  // TODO initialize the mesh so that the starting edge satisfies this property
  // by default

  while (!currHe.isReal()) {
    currHe = currHe.twin().next();
  }
}

// === Adjacent vertices
VertexAdjacentVertexSet::VertexAdjacentVertexSet(HalfedgePtr he)
    : firstHe(he) {}

VertexAdjacentVertexIterator::VertexAdjacentVertexIterator(
    HalfedgePtr startingEdge, bool justStarted_)
    : currHe(startingEdge), justStarted(justStarted_) {}

// === Adjacent faces
VertexAdjacentFaceSet::VertexAdjacentFaceSet(HalfedgePtr he) : firstHe(he) {}

VertexAdjacentFaceIterator::VertexAdjacentFaceIterator(HalfedgePtr startingEdge,
                                                       bool justStarted_)
    : currHe(startingEdge), justStarted(justStarted_) {}

// === Adjacent edges
VertexAdjacentEdgeSet::VertexAdjacentEdgeSet(HalfedgePtr he) : firstHe(he) {}

VertexAdjacentEdgeIterator::VertexAdjacentEdgeIterator(HalfedgePtr startingEdge,
                                                       bool justStarted_)
    : currHe(startingEdge), justStarted(justStarted_) {}

// === Adjacent corners
VertexAdjacentCornerSet::VertexAdjacentCornerSet(HalfedgePtr he) : firstHe(he) {
  // Set firstHe to a real halfedge
  if (!firstHe.isReal()) firstHe = firstHe.twin().next();
}

VertexAdjacentCornerIterator::VertexAdjacentCornerIterator(
    HalfedgePtr startingEdge, bool justStarted_)
    : currHe(startingEdge), justStarted(justStarted_) {}

// ==========================================================
// ================     Face Iterators     ==================
// ==========================================================

// === Adjacent halfedges
FaceAdjacentHalfedgeSet::FaceAdjacentHalfedgeSet(HalfedgePtr he)
    : firstHe(he) {}

FaceAdjacentHalfedgeIterator::FaceAdjacentHalfedgeIterator(
    HalfedgePtr startingEdge, bool justStarted_)
    : currHe(startingEdge), justStarted(justStarted_) {}

// === Adjacent vertices
FaceAdjacentVertexSet::FaceAdjacentVertexSet(HalfedgePtr he) : firstHe(he) {}

FaceAdjacentVertexIterator::FaceAdjacentVertexIterator(HalfedgePtr startingEdge,
                                                       bool justStarted_)
    : currHe(startingEdge), justStarted(justStarted_) {}

// === Adjacent edges
FaceAdjacentEdgeSet::FaceAdjacentEdgeSet(HalfedgePtr he) : firstHe(he) {}

FaceAdjacentEdgeIterator::FaceAdjacentEdgeIterator(HalfedgePtr startingEdge,
                                                   bool justStarted_)
    : currHe(startingEdge), justStarted(justStarted_) {}

// === Adjacent faces
FaceAdjacentFaceSet::FaceAdjacentFaceSet(HalfedgePtr he) : firstHe(he) {}

FaceAdjacentFaceIterator::FaceAdjacentFaceIterator(HalfedgePtr startingEdge,
                                                   bool justStarted_)
    : currHe(startingEdge), justStarted(justStarted_) {
  // face->halfedge could point to a halfedge on the boundary, whose
  // corresponding face should
  // not be returned by this iterator. As such, advance the iterator during
  // construction until
  // it points to a non-border halfedge (or until it loops around).
  // An alternative would be to require that face->halfedge points to a
  // non-boundary halfedge,
  // but this feels like an unnatural requirement (for instance, it disallows a
  // single lonely
  // triangle as a mesh element).
  // As an aside, in the "lonely triangle" case, this code is correct but does
  // an extra
  // circuit of the face.
  // TODO can we initialize the mesh such that this loop is not needed (i.e.,
  // such that it's satisfied by default?)
  if (!currHe.twin().isReal()) {
    do {
      currHe = currHe.next();
    } while (!currHe.twin().isReal() && currHe != startingEdge);
  }
}

// === Adjacent corners
FaceAdjacentCornerSet::FaceAdjacentCornerSet(HalfedgePtr he) : firstHe(he) {}

FaceAdjacentCornerIterator::FaceAdjacentCornerIterator(HalfedgePtr startingEdge,
                                                       bool justStarted_)
    : currHe(startingEdge), justStarted(justStarted_) {}

} // namespace geometrycentral