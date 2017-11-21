#pragma once

#include <vector>

#include <geometrycentral/utilities.h>
#include "geometrycentral/vector3.h"
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

class HalfedgeDual;
class DualHalfedge;
class DualVertex;
class DualEdge;
class DualFace;

// Forward declare classes used for conversion from other mesh types
class PolygonSoupMesh;
template <class T>
class Geometry;

}  // namespace geometrycentral

// Order MATTERS for these includes
// 1
#include "geometrycentral/halfedge_pointer_types.h"
// 2
#include "geometrycentral/halfedge_dual_pointer_types.h"
// rest
#include "geometrycentral/halfedge_data_types.h"
#include "geometrycentral/halfedge_dual_data_types.h"
#include "geometrycentral/halfedge_dual_iterators.h"
#include "geometrycentral/halfedge_iterators.h"

namespace geometrycentral {

class HalfedgeDual;

class HalfedgeMesh {
  friend class HalfedgeDual;

 public:
  HalfedgeMesh();
  HalfedgeMesh(const PolygonSoupMesh& soup, Geometry<Vector3>*& geometry);

  HalfedgeDual dual();

  // Number of mesh elements of each type
  size_t nHalfedges() const;
  size_t nCorners() const;
  size_t nVertices() const;
  size_t
  nInteriorVertices();  // WARNING: O(n) at the moment TODO make this const
  size_t nEdges() const;
  size_t nFaces() const;
  size_t nBoundaryLoops() const;
  size_t nImaginaryHalfedges() const;

  // Methods for range-based for loops
  // Example: for(VertexPtr v : mesh.vertices()) { ... }
  HalfedgePtrSet halfedges();
  CornerPtrSet corners();
  VertexPtrSet vertices();
  EdgePtrSet edges();
  FacePtrSet faces();
  BoundaryPtrSet boundaryLoops();
  HalfedgePtrSet imaginaryHalfedges();
  CutPtrSet cutBoundary(
      int loop = -1);  // NOTE: Designed for 1 boundary loop. With more
                       // than 1 loop, only halfedges on the longest loop
                       // and its associated cuts are visited if loop is not
                       // specified

  // Methods for accessing elements by index
  // Example: VertexPtr v = mesh.vertex(123);
  HalfedgePtr halfedge(size_t index);
  CornerPtr corner(size_t index);
  VertexPtr vertex(size_t index);
  EdgePtr edge(size_t index);
  FacePtr face(size_t index);
  HalfedgePtr imaginaryHalfedge(size_t index);
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
  bool isSimplicial();  // returns true if and only if all faces are triangles
  bool hasCut();        // returns true if any edge is marked as a cut edge
  size_t nFacesTriangulation();  // returns the number of triangles in the
                                 // triangulation determined by
                                 // Face::triangulate()
  size_t longestBoundaryLoop();
  HalfedgeMesh*
  copy();  // returns a deep copy (note: might not preserve all indices)

 private:
  // The contiguous chunks of memory which hold the actual structs.
  // Don't modify them after construction.
  std::vector<Halfedge> rawHalfedges;
  std::vector<Vertex> rawVertices;
  std::vector<Edge> rawEdges;
  std::vector<Face> rawFaces;
  std::vector<Halfedge> rawImaginaryHalfedges;
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
  bool _isSimplicial;
  size_t _nFacesTriangulation;
  size_t _longestBoundaryLoop;
  size_t _nInteriorVertices;
};

class HalfedgeDual {
 public:
  HalfedgeDual(HalfedgeMesh& mesh);

  // Number of mesh elements of each type
  size_t nHalfedges() const;
  size_t nVertices() const;
  size_t nEdges() const;
  size_t nFaces() const;

  // Methods for range-based for loops
  // Example: for(DualVertexPtr v : mesh.vertices()) { ... }
  DualHalfedgePtrSet halfedges();
  DualVertexPtrSet vertices();
  DualEdgePtrSet edges();
  DualFacePtrSet faces();

  // Methods for accessing elements by index
  // Example: VertexPtr v = mesh.vertex[123];
  DualHalfedgePtr halfedge(size_t index);
  DualVertexPtr vertex(size_t index);
  DualEdgePtr edge(size_t index);
  DualFacePtr face(size_t index);

 private:
#ifndef NDEBUG
 public:
#endif
  // References to primal data
  HalfedgeMesh* mesh = nullptr;
  std::vector<Halfedge>& rawHalfedges;
  std::vector<Face>& rawVertices;
  std::vector<Edge>& rawEdges;
  std::vector<Vertex>& rawFaces;
  // Notice that vertices and faces are swapped---the hope is that
  // by doing this swap once (here), it avoids subsequent errors where
  // methods forget to interchange vertices and faces.
};

class Halfedge {
  friend class Edge;
  friend class HalfedgeMesh;
  friend class HalfedgePtr;
  friend class DualHalfedgePtr;

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
  friend class DualFacePtr;

 protected:
  // Data structure
  Halfedge* halfedge;  // some halfedge that emanates from this vertex
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
  friend class DualEdgePtr;

 protected:
  Halfedge* halfedge;

  bool flip(void);

  bool isBoundary = false;
  bool isCut = false;

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
  friend class DualVertexPtr;

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

}  // namespace geometrycentral

#include "geometrycentral/halfedge_data_types.ipp"
#include "geometrycentral/halfedge_dual_data_types.ipp"
#include "geometrycentral/halfedge_dual_iterators.ipp"
#include "geometrycentral/halfedge_dual_pointer_types.ipp"
#include "geometrycentral/halfedge_iterators.ipp"
#include "geometrycentral/halfedge_mesh.ipp"
#include "geometrycentral/halfedge_pointer_types.ipp"
