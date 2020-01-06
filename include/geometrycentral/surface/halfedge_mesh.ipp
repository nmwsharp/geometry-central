#pragma once

namespace geometrycentral {
namespace surface {

// clang-format off

// Methods for getting number of mesh elements
inline size_t HalfedgeMesh::nHalfedges()         const { return nHalfedgesCount; }
inline size_t HalfedgeMesh::nInteriorHalfedges() const { return nInteriorHalfedgesCount; }
inline size_t HalfedgeMesh::nExteriorHalfedges() const { return nHalfedgesCount - nInteriorHalfedgesCount; }
inline size_t HalfedgeMesh::nCorners()           const { return nInteriorHalfedgesCount; }
inline size_t HalfedgeMesh::nVertices()          const { return nVerticesCount; }
inline size_t HalfedgeMesh::nEdges()             const { return nHalfedgesCount / 2; }
inline size_t HalfedgeMesh::nFaces()             const { return nFacesCount; }
inline size_t HalfedgeMesh::nBoundaryLoops()     const { return nBoundaryLoopsCount; }

// Capacities
inline size_t HalfedgeMesh::nHalfedgesCapacity()        const { return nHalfedgesCapacityCount; }
inline size_t HalfedgeMesh::nVerticesCapacity()         const { return nVerticesCapacityCount; }
inline size_t HalfedgeMesh::nEdgesCapacity()            const { return nHalfedgesCapacityCount / 2; }
inline size_t HalfedgeMesh::nFacesCapacity()            const { return nFacesCapacityCount - nBoundaryLoopsFillCount; }
inline size_t HalfedgeMesh::nBoundaryLoopsCapacity()    const { return nFacesCapacityCount - nFacesFillCount; }

// Implicit relationships
inline size_t HalfedgeMesh::heTwin(size_t iHe)   { return iHe ^ 1; }
inline size_t HalfedgeMesh::heEdge(size_t iHe)   { return iHe / 2; }
inline size_t HalfedgeMesh::eHalfedge(size_t iE) { return 2 * iE; }
inline size_t HalfedgeMesh::nEdgesFillCount() const { return nHalfedgesFillCount/2; }

// Other getters
inline bool HalfedgeMesh::heIsInterior(size_t iHe) const { return !faceIsBoundaryLoop(heFace[iHe]); }
inline bool HalfedgeMesh::faceIsBoundaryLoop(size_t iF) const { return iF >= nFacesFillCount; }
inline size_t HalfedgeMesh::faceIndToBoundaryLoopInd(size_t iF) const { return nFacesCapacityCount - 1 - iF;}
inline size_t HalfedgeMesh::boundaryLoopIndToFaceInd(size_t iB) const { return nFacesCapacityCount - 1 - iB;}

// Detect dead elements
inline bool HalfedgeMesh::vertexIsDead(size_t iV)      const { return vHalfedge[iV] == INVALID_IND; }
inline bool HalfedgeMesh::halfedgeIsDead(size_t iHe)   const { return heNext[iHe] == INVALID_IND; }
inline bool HalfedgeMesh::edgeIsDead(size_t iE)        const { return heNext[eHalfedge(iE)] == INVALID_IND; }
inline bool HalfedgeMesh::faceIsDead(size_t iF)        const { return fHalfedge[iF] == INVALID_IND;}

// Methods for iterating over mesh elements w/ range-based for loops ===========

inline VertexSet HalfedgeMesh::vertices()                        { return VertexSet(this, 0, nVerticesFillCount); }
inline HalfedgeSet HalfedgeMesh::halfedges()                     { return HalfedgeSet(this, 0, nHalfedgesFillCount); }
inline HalfedgeInteriorSet HalfedgeMesh::interiorHalfedges()     { return HalfedgeInteriorSet(this, 0, nHalfedgesFillCount); }
inline HalfedgeExteriorSet HalfedgeMesh::exteriorHalfedges()     { return HalfedgeExteriorSet(this, 0, nHalfedgesFillCount); }
inline CornerSet HalfedgeMesh::corners()                         { return CornerSet(this, 0, nHalfedgesFillCount); }
inline EdgeSet HalfedgeMesh::edges()                             { return EdgeSet(this, 0, nEdgesFillCount()); }
inline FaceSet HalfedgeMesh::faces()                             { return FaceSet(this, 0, nFacesFillCount); }
inline BoundaryLoopSet HalfedgeMesh::boundaryLoops()             { return BoundaryLoopSet(this, 0, nBoundaryLoopsFillCount); }

// Methods for accessing elements by index =====================================
// Note that these are only valid when the mesh is compressed.

inline Vertex HalfedgeMesh::vertex(size_t index)             { return Vertex(this, index); }
inline Halfedge HalfedgeMesh::halfedge(size_t index)         { return Halfedge(this, index); } 
inline Corner HalfedgeMesh::corner(size_t index)             { return Corner(this,index); }
inline Edge HalfedgeMesh::edge(size_t index)                 { return Edge(this, index); }
inline Face HalfedgeMesh::face(size_t index)                 { return Face(this, index); }
inline BoundaryLoop HalfedgeMesh::boundaryLoop(size_t index) { return BoundaryLoop(this, index); }

// Misc utility methods =====================================

inline bool HalfedgeMesh::isCompressed() const { return isCompressedFlag; }
inline bool HalfedgeMesh::hasBoundary() const { return nBoundaryLoopsCount > 0; }

// clang-format on


} // namespace surface
} // namespace geometrycentral
