#pragma once

namespace geometrycentral {
namespace surface {

// ==========================================================
// ================     Base  Iterator     ==================
// ==========================================================


template <typename N>
inline NavigationIteratorBase<N>::NavigationIteratorBase(typename N::Etype e, bool justStarted_)
    : state{e}, justStarted(justStarted_) {

  // checking startingE is a weird hack, but ensures the iterators don't infinite loop if there are no valid elemetnts
  // to return
  typename N::Etype startingE;

  // Advance to first valid element
  while (!state.isValid()) {
    state.advance();
    if (state.currE == startingE) break;
  }
}

template <typename N>
inline const NavigationIteratorBase<N>& NavigationIteratorBase<N>::operator++() {
  state.advance();
  while (!state.isValid()) {
    state.advance();
  }
  justStarted = false;
  return *this;
}

template <typename N>
inline bool NavigationIteratorBase<N>::operator==(const NavigationIteratorBase<N>& other) const {
  return justStarted == other.justStarted && state.currE == other.state.currE;
}

template <typename N>
inline bool NavigationIteratorBase<N>::operator!=(const NavigationIteratorBase<N>& other) const {
  return !(*this == other);
}

template <typename N>
inline typename N::Rtype NavigationIteratorBase<N>::operator*() const {
  return state.getCurrent();
}

template <typename N>
NavigationSetBase<N>::NavigationSetBase(typename N::Etype firstE_)
    : firstE(firstE_), cachedEndIter{NavigationIteratorBase<N>(firstE, false)} {}

template <typename N>
inline NavigationIteratorBase<N> NavigationSetBase<N>::begin() const {
  return NavigationIteratorBase<N>(firstE, true);
}

template <typename N>
inline NavigationIteratorBase<N> NavigationSetBase<N>::end() const {
  return cachedEndIter;
}


// ==========================================================
// ================    Vertex Iterators    ==================
// ==========================================================

inline void VertexAdjacentVertexNavigator::advance() { currE = currE.twin().next(); }
inline bool VertexAdjacentVertexNavigator::isValid() const { return true; }
inline Vertex VertexAdjacentVertexNavigator::getCurrent() const { return currE.twin().vertex(); }

inline void VertexIncomingHalfedgeNavigator::advance() { currE = currE.twin().next(); }
inline bool VertexIncomingHalfedgeNavigator::isValid() const { return true; }
inline Halfedge VertexIncomingHalfedgeNavigator::getCurrent() const { return currE.twin(); }

inline void VertexOutgoingHalfedgeNavigator::advance() { currE = currE.twin().next(); }
inline bool VertexOutgoingHalfedgeNavigator::isValid() const { return true; }
inline Halfedge VertexOutgoingHalfedgeNavigator::getCurrent() const { return currE; }

inline void VertexAdjacentCornerNavigator::advance() { currE = currE.twin().next(); }
inline bool VertexAdjacentCornerNavigator::isValid() const { return currE.isInterior(); }
inline Corner VertexAdjacentCornerNavigator::getCurrent() const { return currE.corner(); }

inline void VertexAdjacentEdgeNavigator::advance() { currE = currE.twin().next(); }
inline bool VertexAdjacentEdgeNavigator::isValid() const { return true; }
inline Edge VertexAdjacentEdgeNavigator::getCurrent() const { return currE.edge(); }

inline void VertexAdjacentFaceNavigator::advance() { currE = currE.twin().next(); }
inline bool VertexAdjacentFaceNavigator::isValid() const { return currE.isInterior(); }
inline Face VertexAdjacentFaceNavigator::getCurrent() const { return currE.face(); }


// ==========================================================
// ================     Face Iterators     ==================
// ==========================================================

inline void FaceAdjacentVertexNavigator::advance() { currE = currE.next(); }
inline bool FaceAdjacentVertexNavigator::isValid() const { return true; }
inline Vertex FaceAdjacentVertexNavigator::getCurrent() const { return currE.vertex(); }

inline void FaceAdjacentHalfedgeNavigator::advance() { currE = currE.next(); }
inline bool FaceAdjacentHalfedgeNavigator::isValid() const { return true; }
inline Halfedge FaceAdjacentHalfedgeNavigator::getCurrent() const { return currE; }

inline void FaceAdjacentCornerNavigator::advance() { currE = currE.next(); }
inline bool FaceAdjacentCornerNavigator::isValid() const { return true; }
inline Corner FaceAdjacentCornerNavigator::getCurrent() const { return currE.corner(); }

inline void FaceAdjacentEdgeNavigator::advance() { currE = currE.next(); }
inline bool FaceAdjacentEdgeNavigator::isValid() const { return true; }
inline Edge FaceAdjacentEdgeNavigator::getCurrent() const { return currE.edge(); }

inline void FaceAdjacentFaceNavigator::advance() { currE = currE.next(); }
inline bool FaceAdjacentFaceNavigator::isValid() const { return currE.twin().isInterior(); }
inline Face FaceAdjacentFaceNavigator::getCurrent() const { return currE.twin().face(); }

// ==========================================================
// ==============   Boundary Loop Iterators   ===============
// ==========================================================

inline void BoundaryLoopAdjacentVertexNavigator::advance() { currE = currE.next(); }
inline bool BoundaryLoopAdjacentVertexNavigator::isValid() const { return true; }
inline Vertex BoundaryLoopAdjacentVertexNavigator::getCurrent() const { return currE.vertex(); }

inline void BoundaryLoopAdjacentHalfedgeNavigator::advance() { currE = currE.next(); }
inline bool BoundaryLoopAdjacentHalfedgeNavigator::isValid() const { return true; }
inline Halfedge BoundaryLoopAdjacentHalfedgeNavigator::getCurrent() const { return currE; }

inline void BoundaryLoopAdjacentEdgeNavigator::advance() { currE = currE.next(); }
inline bool BoundaryLoopAdjacentEdgeNavigator::isValid() const { return true; }
inline Edge BoundaryLoopAdjacentEdgeNavigator::getCurrent() const { return currE.edge(); }

} // namespace surface
} // namespace geometrycentral
