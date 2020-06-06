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
inline size_t HalfedgeMesh::nEdges()             const { return nEdgesCount; }
inline size_t HalfedgeMesh::nFaces()             const { return nFacesCount; }
inline size_t HalfedgeMesh::nBoundaryLoops()     const { return nBoundaryLoopsCount; }

// Capacities
inline size_t HalfedgeMesh::nHalfedgesCapacity()        const { return nHalfedgesCapacityCount; }
inline size_t HalfedgeMesh::nVerticesCapacity()         const { return nVerticesCapacityCount; }
inline size_t HalfedgeMesh::nEdgesCapacity()            const { return nEdgesCapacityCount; }
inline size_t HalfedgeMesh::nFacesCapacity()            const { return nFacesCapacityCount - nBoundaryLoopsFillCount; }
inline size_t HalfedgeMesh::nBoundaryLoopsCapacity()    const { return nFacesCapacityCount - nFacesFillCount; }

// Connectivity
inline size_t HalfedgeMesh::heNext(size_t iHe)      const { return heNextArr[iHe]; }
inline size_t HalfedgeMesh::heTwin(size_t iHe)      const { if(usesImplictTwin()) return heTwinImplicit(iHe); 
                                                            //throw std::runtime_error("called he.twin() on not-necessarily-manifold mesh. Try he.sibling() instead"); 
                                                            return heSiblingArr[iHe]; }
inline size_t HalfedgeMesh::heSibling(size_t iHe)   const { return usesImplictTwin() ? heTwinImplicit(iHe) : heSiblingArr[iHe]; }
inline size_t HalfedgeMesh::heEdge(size_t iHe)      const { return usesImplictTwin() ? heEdgeImplicit(iHe) : heEdgeArr[iHe]; }
inline size_t HalfedgeMesh::heVertex(size_t iHe)    const { return heVertexArr[iHe]; }
inline size_t HalfedgeMesh::heFace(size_t iHe)      const { return heFaceArr[iHe]; }
inline size_t HalfedgeMesh::eHalfedge(size_t iE)    const { return usesImplictTwin() ? eHalfedgeImplicit(iE) : eHalfedgeArr[iE]; }
inline size_t HalfedgeMesh::vHalfedge(size_t iV)    const { return vHalfedgeArr[iV]; }
inline size_t HalfedgeMesh::fHalfedge(size_t iF)    const { return fHalfedgeArr[iF]; }

// Implicit relationships
inline bool HalfedgeMesh::usesImplictTwin() const           { return useImplicitTwinFlag; }
inline size_t HalfedgeMesh::heTwinImplicit(size_t iHe)      { return iHe ^ 1; }     // static
inline size_t HalfedgeMesh::heEdgeImplicit(size_t iHe)      { return iHe / 2; }     // static
inline size_t HalfedgeMesh::eHalfedgeImplicit(size_t iE)    { return 2 * iE; }      // static

// Other getters
inline bool HalfedgeMesh::heIsInterior(size_t iHe) const { return !faceIsBoundaryLoop(heFaceArr[iHe]); }
inline bool HalfedgeMesh::faceIsBoundaryLoop(size_t iF) const { return iF >= nFacesFillCount; }
inline size_t HalfedgeMesh::faceIndToBoundaryLoopInd(size_t iF) const { return nFacesCapacityCount - 1 - iF;}
inline size_t HalfedgeMesh::boundaryLoopIndToFaceInd(size_t iB) const { return nFacesCapacityCount - 1 - iB;}

// Detect dead elements
inline bool HalfedgeMesh::vertexIsDead(size_t iV)      const { return vHalfedgeArr[iV] == INVALID_IND; }
inline bool HalfedgeMesh::halfedgeIsDead(size_t iHe)   const { return heNextArr[iHe] == INVALID_IND; }
inline bool HalfedgeMesh::edgeIsDead(size_t iE)        const { return usesImplictTwin() ? heNextArr[eHalfedgeImplicit(iE)] == INVALID_IND : eHalfedgeArr[iE] == INVALID_IND; }
inline bool HalfedgeMesh::faceIsDead(size_t iF)        const { return fHalfedgeArr[iF] == INVALID_IND;}

// Methods for iterating over mesh elements w/ range-based for loops ===========

inline VertexSet HalfedgeMesh::vertices()                        { return VertexSet(this, 0, nVerticesFillCount); }
inline HalfedgeSet HalfedgeMesh::halfedges()                     { return HalfedgeSet(this, 0, nHalfedgesFillCount); }
inline HalfedgeInteriorSet HalfedgeMesh::interiorHalfedges()     { return HalfedgeInteriorSet(this, 0, nHalfedgesFillCount); }
inline HalfedgeExteriorSet HalfedgeMesh::exteriorHalfedges()     { return HalfedgeExteriorSet(this, 0, nHalfedgesFillCount); }
inline CornerSet HalfedgeMesh::corners()                         { return CornerSet(this, 0, nHalfedgesFillCount); }
inline EdgeSet HalfedgeMesh::edges()                             { return EdgeSet(this, 0, nEdgesFillCount); }
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
inline void HalfedgeMesh::ensureVertexIterationCachePopulated() { if(vertexIterationCacheTick != modificationTick) populateVertexIterationCache(); }

// clang-format on


} // namespace surface
} // namespace geometrycentral
