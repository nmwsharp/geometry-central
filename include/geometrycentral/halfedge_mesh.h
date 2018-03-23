#pragma once

#include <vector>

#include "geometrycentral/vector3.h"
#include <geometrycentral/utilities.h>
// NOTE: More includes at bottom of file

namespace geometrycentral {

// A HalfedgeMesh encodes the connectivity---but not the geometry---of a
// manifold surface, possibly with boundary.

// Forward declare classes for primary mesh objects
class HalfedgeMesh;
class Halfedge;
class Vertex;
class Edge;
class Face;

// Forward declare classes used for conversion from other mesh types
class PolygonSoupMesh;
template <class T>
class Geometry;

} // namespace geometrycentral

// Order MATTERS for these includes
// 1
#include "geometrycentral/halfedge_pointer_types.h"
// 3
#include "geometrycentral/halfedge_data_types.h"
#include "geometrycentral/halfedge_iterators.h"
// 4
#include "geometrycentral/halfedge_mesh_data_transfer.h"

namespace geometrycentral {


class HalfedgeMesh {

public:
  HalfedgeMesh();
  HalfedgeMesh(const PolygonSoupMesh& soup, Geometry<Vector3>*& geometry);


  // Number of mesh elements of each type
  size_t nHalfedges() const;
  size_t nCorners() const;
  size_t nVertices() const;
  size_t nInteriorVertices() const;
  size_t nEdges() const;
  size_t nFaces() const;
  size_t nBoundaryLoops() const;
  size_t nImaginaryHalfedges() const;

  // Methods for range-based for loops
  // Example: for(VertexPtr v : mesh.vertices()) { ... }
  HalfedgePtrSet realHalfedges();
  HalfedgePtrSet imaginaryHalfedges();
  HalfedgePtrSet allHalfedges();
  CornerPtrSet corners();
  VertexPtrSet vertices();
  EdgePtrSet edges();
  FacePtrSet faces();
  BoundaryPtrSet boundaryLoops();

  // Methods for accessing elements by index
  // Example: VertexPtr v = mesh.vertex(123);
  HalfedgePtr realHalfedge(size_t index);
  HalfedgePtr imaginaryHalfedge(size_t index);
  HalfedgePtr allHalfedge(size_t index);
  CornerPtr corner(size_t index);
  VertexPtr vertex(size_t index);
  EdgePtr edge(size_t index);
  FacePtr face(size_t index);
  BoundaryPtr boundaryLoop(size_t index);

  // Methods for obtaining canonical indices for mesh elements
  // (Note that in some situations, custom indices might instead be needed)
  VertexData<size_t> getVertexIndices();
  VertexData<size_t> getInteriorVertexIndices();
  FaceData<size_t> getFaceIndices();
  EdgeData<size_t> getEdgeIndices();
  HalfedgeData<size_t> getHalfedgeIndices();
  CornerData<size_t> getCornerIndices();

  // Utility functions
  bool isSimplicial() const;          // returns true if and only if all faces are triangles
  size_t nFacesTriangulation() const; // returns the number of triangles in the
                                      // triangulation determined by
                                      // Face::triangulate()
  size_t longestBoundaryLoop() const;
  int eulerCharacteristic() const;
  size_t nConnectedComponents() const;
  HalfedgeMesh* copy();                            // returns a deep copy
  HalfedgeMesh* copy(HalfedgeMeshDataTransfer& t); // returns a deep copy

private:
  // The contiguous chunks of memory which hold the actual structs.
  // Don't modify them after construction.
  size_t nRealHalfedges;
  std::vector<Halfedge> rawHalfedges; // first real, then imaginary
  std::vector<Vertex> rawVertices;
  std::vector<Edge> rawEdges;
  std::vector<Face> rawFaces;
  std::vector<Face> rawBoundaryLoops;

  // Hide copy and move constructors, we don't wanna mess with that
  HalfedgeMesh(const HalfedgeMesh& other) = delete;
  HalfedgeMesh& operator=(const HalfedgeMesh& other) = delete;
  HalfedgeMesh(HalfedgeMesh&& other) = delete;
  HalfedgeMesh& operator=(HalfedgeMesh&& other) = delete;

  // Cache some basic information that may be queried many
  // times, but require O(n) computation to determine.
  void cacheInfo();
  void cache_isSimplicial();
  void cache_nFacesTriangulation();
  void cache_longestBoundaryLoop();
  void cache_nInteriorVertices();
  void cache_nConnectedComponents();
  bool _isSimplicial;
  size_t _nFacesTriangulation;
  size_t _longestBoundaryLoop;
  size_t _nInteriorVertices;
  size_t _nConnectedComponents;
};

class Halfedge {
  friend class Edge;
  friend class HalfedgeMesh;
  friend class HalfedgePtr;

protected:
  Halfedge* twin;
  Halfedge* next;
  Vertex* vertex;
  Edge* edge;
  Face* face;

  bool isReal = true;

#ifndef NDEBUG
public:
  // The mesh that this is a part of. Should only be used for debugging, so
  // exclude it unless debug is enabled.
  HalfedgeMesh* parentMesh;
#endif
};

class Vertex {
  friend class Edge;
  friend class HalfedgeMesh;
  friend class VertexPtr;

protected:
  // Data structure
  Halfedge* halfedge; // some halfedge that emanates from this vertex
                      // (guaranteed to be real)
  bool isBoundary = false;

#ifndef NDEBUG
public:
  // The mesh that this is a part of. Should only be used for debugging, so
  // exclude it unless debug is enabled.
  HalfedgeMesh* parentMesh;
#endif

  Vector3 projectToTangentSpace(Vector3 inVec);
};

class Edge {
  friend class HalfedgeMesh;
  friend class EdgePtr;

protected:
  Halfedge* halfedge;

  bool flip(void);

  bool isBoundary = false;

#ifndef NDEBUG
public:
  // The mesh that this is a part of. Should only be used for debugging, so
  // exclude it unless debug is enabled.
  HalfedgeMesh* parentMesh;
#endif
};

class Face {
  friend class Edge;
  friend class HalfedgeMesh;
  friend class FacePtr;

protected:
  Halfedge* halfedge;

  bool isBoundary = false;
  bool isReal = false;

#ifndef NDEBUG
public:
  // The mesh that this is a part of. Should only be used for debugging, so
  // exclude it unless debug is enabled.
  HalfedgeMesh* parentMesh;
#endif
};

} // namespace geometrycentral

#include "geometrycentral/halfedge_data_types.ipp"
#include "geometrycentral/halfedge_iterators.ipp"
#include "geometrycentral/halfedge_mesh.ipp"
#include "geometrycentral/halfedge_mesh_data_transfer.ipp"
#include "geometrycentral/halfedge_pointer_types.ipp"
