#pragma once

// ==============================================================
// ================    DualVertex Iterators    ==================
// ==============================================================

// === Incoming dual halfedges, including those on the domain boundary
inline DualVertexIncomingDualHalfedgeIterator
DualVertexIncomingDualHalfedgeSet::begin() {
  return DualVertexIncomingDualHalfedgeIterator(firstHe, true);
}
inline DualVertexIncomingDualHalfedgeIterator
DualVertexIncomingDualHalfedgeSet::end() {
  return DualVertexIncomingDualHalfedgeIterator(firstHe, false);
}
inline const DualVertexIncomingDualHalfedgeIterator&
    DualVertexIncomingDualHalfedgeIterator::
    operator++() {
  justStarted = false;
  currHe = currHe.next().twin();
  return *this;
}
inline bool DualVertexIncomingDualHalfedgeIterator::operator==(
    const DualVertexIncomingDualHalfedgeIterator& other) const {
  return currHe == other.currHe && justStarted == other.justStarted;
}
inline bool DualVertexIncomingDualHalfedgeIterator::operator!=(
    const DualVertexIncomingDualHalfedgeIterator& other) const {
  return !(*this == other);
}
inline DualHalfedgePtr DualVertexIncomingDualHalfedgeIterator::operator*()
    const {
  return currHe;
}

// === Outgoing dual halfedges, including those on the domain boundary
inline DualVertexOutgoingDualHalfedgeIterator
DualVertexOutgoingDualHalfedgeSet::begin() {
  return DualVertexOutgoingDualHalfedgeIterator(firstHe, true);
}
inline DualVertexOutgoingDualHalfedgeIterator
DualVertexOutgoingDualHalfedgeSet::end() {
  return DualVertexOutgoingDualHalfedgeIterator(firstHe, false);
}
inline const DualVertexOutgoingDualHalfedgeIterator&
    DualVertexOutgoingDualHalfedgeIterator::
    operator++() {
  justStarted = false;
  currHe = currHe.twin().next();
  return *this;
}
inline bool DualVertexOutgoingDualHalfedgeIterator::operator==(
    const DualVertexOutgoingDualHalfedgeIterator& other) const {
  return currHe == other.currHe && justStarted == other.justStarted;
}
inline bool DualVertexOutgoingDualHalfedgeIterator::operator!=(
    const DualVertexOutgoingDualHalfedgeIterator& other) const {
  return !(*this == other);
}
inline DualHalfedgePtr DualVertexOutgoingDualHalfedgeIterator::operator*()
    const {
  return currHe;
}

// === Adjacent dual vertices
inline DualVertexAdjacentDualVertexIterator
DualVertexAdjacentDualVertexSet::begin() {
  return DualVertexAdjacentDualVertexIterator(firstHe, true);
}
inline DualVertexAdjacentDualVertexIterator
DualVertexAdjacentDualVertexSet::end() {
  return DualVertexAdjacentDualVertexIterator(firstHe, false);
}
inline const DualVertexAdjacentDualVertexIterator&
    DualVertexAdjacentDualVertexIterator::
    operator++() {
  justStarted = false;
  currHe = currHe.next().twin();
  return *this;
}
inline bool DualVertexAdjacentDualVertexIterator::operator==(
    const DualVertexAdjacentDualVertexIterator& other) const {
  return currHe == other.currHe && justStarted == other.justStarted;
}
inline bool DualVertexAdjacentDualVertexIterator::operator!=(
    const DualVertexAdjacentDualVertexIterator& other) const {
  return !(*this == other);
}
inline DualVertexPtr DualVertexAdjacentDualVertexIterator::operator*() const {
  return currHe.vertex();
}

// === Adjacent dual faces
inline DualVertexAdjacentDualFaceIterator
DualVertexAdjacentDualFaceSet::begin() {
  return DualVertexAdjacentDualFaceIterator(firstHe, true);
}
inline DualVertexAdjacentDualFaceIterator DualVertexAdjacentDualFaceSet::end() {
  return DualVertexAdjacentDualFaceIterator(firstHe, false);
}
inline const DualVertexAdjacentDualFaceIterator&
    DualVertexAdjacentDualFaceIterator::
    operator++() {
  justStarted = false;
  do {
    currHe = currHe.twin().next();
  } while (!currHe.isReal());
  return *this;
}
inline bool DualVertexAdjacentDualFaceIterator::operator==(
    const DualVertexAdjacentDualFaceIterator& other) const {
  return currHe == other.currHe && justStarted == other.justStarted;
}
inline bool DualVertexAdjacentDualFaceIterator::operator!=(
    const DualVertexAdjacentDualFaceIterator& other) const {
  return !(*this == other);
}
inline DualFacePtr DualVertexAdjacentDualFaceIterator::operator*() const {
  return currHe.face();
}

// === Adjacent dual edges
inline DualVertexAdjacentDualEdgeIterator
DualVertexAdjacentDualEdgeSet::begin() {
  return DualVertexAdjacentDualEdgeIterator(firstHe, true);
}
inline DualVertexAdjacentDualEdgeIterator DualVertexAdjacentDualEdgeSet::end() {
  return DualVertexAdjacentDualEdgeIterator(firstHe, false);
}
inline const DualVertexAdjacentDualEdgeIterator&
    DualVertexAdjacentDualEdgeIterator::
    operator++() {
  justStarted = false;
  currHe = currHe.twin().next();
  return *this;
}
inline bool DualVertexAdjacentDualEdgeIterator::operator==(
    const DualVertexAdjacentDualEdgeIterator& other) const {
  return currHe == other.currHe && justStarted == other.justStarted;
}
inline bool DualVertexAdjacentDualEdgeIterator::operator!=(
    const DualVertexAdjacentDualEdgeIterator& other) const {
  return !(*this == other);
}
inline DualEdgePtr DualVertexAdjacentDualEdgeIterator::operator*() const {
  return currHe.edge();
}

// ===============================================================
// ================     Dual Face Iterators     ==================
// ===============================================================

// === Adjacent dual halfedges
inline DualFaceAdjacentDualHalfedgeIterator
DualFaceAdjacentDualHalfedgeSet::begin() {
  return DualFaceAdjacentDualHalfedgeIterator(firstHe, true);
}
inline DualFaceAdjacentDualHalfedgeIterator
DualFaceAdjacentDualHalfedgeSet::end() {
  return DualFaceAdjacentDualHalfedgeIterator(firstHe, false);
}
inline const DualFaceAdjacentDualHalfedgeIterator&
    DualFaceAdjacentDualHalfedgeIterator::
    operator++() {
  justStarted = false;
  currHe = currHe.next();
  return *this;
}
inline bool DualFaceAdjacentDualHalfedgeIterator::operator==(
    const DualFaceAdjacentDualHalfedgeIterator& other) const {
  return currHe == other.currHe && justStarted == other.justStarted;
}
inline bool DualFaceAdjacentDualHalfedgeIterator::operator!=(
    const DualFaceAdjacentDualHalfedgeIterator& other) const {
  return !(*this == other);
}
inline DualHalfedgePtr DualFaceAdjacentDualHalfedgeIterator::operator*() const {
  return currHe;
}

// === Adjacent dual vertices
inline DualFaceAdjacentDualVertexIterator
DualFaceAdjacentDualVertexSet::begin() {
  return DualFaceAdjacentDualVertexIterator(firstHe, true);
}
inline DualFaceAdjacentDualVertexIterator DualFaceAdjacentDualVertexSet::end() {
  return DualFaceAdjacentDualVertexIterator(firstHe, false);
}
inline const DualFaceAdjacentDualVertexIterator&
    DualFaceAdjacentDualVertexIterator::
    operator++() {
  justStarted = false;
  currHe = currHe.next();
  return *this;
}
inline bool DualFaceAdjacentDualVertexIterator::operator==(
    const DualFaceAdjacentDualVertexIterator& other) const {
  return currHe == other.currHe && justStarted == other.justStarted;
}
inline bool DualFaceAdjacentDualVertexIterator::operator!=(
    const DualFaceAdjacentDualVertexIterator& other) const {
  return !(*this == other);
}
inline DualVertexPtr DualFaceAdjacentDualVertexIterator::operator*() const {
  return currHe.vertex();
}

// === Adjacent dual edges
inline DualFaceAdjacentDualEdgeIterator DualFaceAdjacentDualEdgeSet::begin() {
  return DualFaceAdjacentDualEdgeIterator(firstHe, true);
}
inline DualFaceAdjacentDualEdgeIterator DualFaceAdjacentDualEdgeSet::end() {
  return DualFaceAdjacentDualEdgeIterator(firstHe, false);
}
inline const DualFaceAdjacentDualEdgeIterator&
    DualFaceAdjacentDualEdgeIterator::
    operator++() {
  justStarted = false;
  currHe = currHe.next();
  return *this;
}
inline bool DualFaceAdjacentDualEdgeIterator::operator==(
    const DualFaceAdjacentDualEdgeIterator& other) const {
  return currHe == other.currHe && justStarted == other.justStarted;
}
inline bool DualFaceAdjacentDualEdgeIterator::operator!=(
    const DualFaceAdjacentDualEdgeIterator& other) const {
  return !(*this == other);
}
inline DualEdgePtr DualFaceAdjacentDualEdgeIterator::operator*() const {
  return currHe.edge();
}

// === Adjacent dual faces
inline DualFaceAdjacentDualFaceIterator DualFaceAdjacentDualFaceSet::begin() {
  return DualFaceAdjacentDualFaceIterator(firstHe, true);
}
inline DualFaceAdjacentDualFaceIterator DualFaceAdjacentDualFaceSet::end() {
  return DualFaceAdjacentDualFaceIterator(firstHe, false);
}
inline const DualFaceAdjacentDualFaceIterator&
    DualFaceAdjacentDualFaceIterator::
    operator++() {
  justStarted = false;
  do {
    currHe = currHe.next();
  } while (!currHe.twin().isReal());
  return *this;
}
inline bool DualFaceAdjacentDualFaceIterator::operator==(
    const DualFaceAdjacentDualFaceIterator& other) const {
  return currHe == other.currHe && justStarted == other.justStarted;
}
inline bool DualFaceAdjacentDualFaceIterator::operator!=(
    const DualFaceAdjacentDualFaceIterator& other) const {
  return !(*this == other);
}
inline DualFacePtr DualFaceAdjacentDualFaceIterator::operator*() const {
  return currHe.twin().face();
}
