#include <halfedge_mesh.h>

// Implementation of various methods to support dual iterators in the halfedge
// mesh.
// Note that many functions are inlined, and defined in halfedge_iterators.ipp

// ==============================================================
// ================    DualVertex Iterators    ==================
// ==============================================================

// === Incoming halfedges (with imaginary)
DualVertexIncomingDualHalfedgeSet::DualVertexIncomingDualHalfedgeSet(
    DualHalfedgePtr he)
    : firstHe(he) {}

DualVertexIncomingDualHalfedgeIterator::DualVertexIncomingDualHalfedgeIterator(
    DualHalfedgePtr startingDualEdge, bool justStarted_)
    : currHe(startingDualEdge), justStarted(justStarted_) {}

// === Adjacent vertices
DualVertexAdjacentDualVertexSet::DualVertexAdjacentDualVertexSet(
    DualHalfedgePtr he)
    : firstHe(he) {}

DualVertexAdjacentDualVertexIterator::DualVertexAdjacentDualVertexIterator(
    DualHalfedgePtr startingDualEdge, bool justStarted_)
    : currHe(startingDualEdge), justStarted(justStarted_) {}

// === Adjacent faces
DualVertexAdjacentDualFaceSet::DualVertexAdjacentDualFaceSet(DualHalfedgePtr he)
    : firstHe(he) {}

DualVertexAdjacentDualFaceIterator::DualVertexAdjacentDualFaceIterator(
    DualHalfedgePtr startingDualEdge, bool justStarted_)
    : currHe(startingDualEdge), justStarted(justStarted_) {}

// === Adjacent edges
DualVertexAdjacentDualEdgeSet::DualVertexAdjacentDualEdgeSet(DualHalfedgePtr he)
    : firstHe(he) {}

DualVertexAdjacentDualEdgeIterator::DualVertexAdjacentDualEdgeIterator(
    DualHalfedgePtr startingDualEdge, bool justStarted_)
    : currHe(startingDualEdge), justStarted(justStarted_) {}

// ==============================================================
// ================     DualFace Iterators     ==================
// ==============================================================

// === Adjacent halfedges
DualFaceAdjacentDualHalfedgeSet::DualFaceAdjacentDualHalfedgeSet(
    DualHalfedgePtr he)
    : firstHe(he) {}

DualFaceAdjacentDualHalfedgeIterator::DualFaceAdjacentDualHalfedgeIterator(
    DualHalfedgePtr startingDualEdge, bool justStarted_)
    : currHe(startingDualEdge), justStarted(justStarted_) {}

// === Adjacent vertices
DualFaceAdjacentDualVertexSet::DualFaceAdjacentDualVertexSet(DualHalfedgePtr he)
    : firstHe(he) {}

DualFaceAdjacentDualVertexIterator::DualFaceAdjacentDualVertexIterator(
    DualHalfedgePtr startingDualEdge, bool justStarted_)
    : currHe(startingDualEdge), justStarted(justStarted_) {}

// === Adjacent edges
DualFaceAdjacentDualEdgeSet::DualFaceAdjacentDualEdgeSet(DualHalfedgePtr he)
    : firstHe(he) {}

DualFaceAdjacentDualEdgeIterator::DualFaceAdjacentDualEdgeIterator(
    DualHalfedgePtr startingDualEdge, bool justStarted_)
    : currHe(startingDualEdge), justStarted(justStarted_) {}

// === Adjacent faces
DualFaceAdjacentDualFaceSet::DualFaceAdjacentDualFaceSet(DualHalfedgePtr he)
    : firstHe(he) {}

DualFaceAdjacentDualFaceIterator::DualFaceAdjacentDualFaceIterator(
    DualHalfedgePtr startingDualEdge, bool justStarted_)
    : currHe(startingDualEdge), justStarted(justStarted_) {}
