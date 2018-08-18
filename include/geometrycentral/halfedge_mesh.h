#pragma once

#include <list>
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
  HalfedgePtr halfedge(size_t index);
  CornerPtr corner(size_t index);
  VertexPtr vertex(size_t index);
  EdgePtr edge(size_t index);
  FacePtr face(size_t index);
  BoundaryPtr boundaryLoop(size_t index);

  // Methods that mutate the mesh. Note that these occasionally trigger a resize, which invaliates
  // any outstanding VertexPtr or MeshData<> objects.
  // TODOs: support removing elements, support adding boundary

  // Adds a vertex along an edge, increasing degree of faces. Returns ptr along the new edge, with he.vertex() as new
  // vertex
  HalfedgePtr insertVertexAlongEdge(EdgePtr e);

  // Split an edge, also splitting adjacent faces. Returns new vertex.
  VertexPtr splitEdge(EdgePtr e);

  // Add vertex inside face and triangulate. Returns new vertex.
  VertexPtr insertVertex(FacePtr f);

  // Add an edge connecting two vertices inside the same face. Returns new halfedge with vA at tail. he.twin().face() is the new face.
  HalfedgePtr connectVertices(VertexPtr vA, VertexPtr vB);

  // Same as above. Faster if you know the face.
  HalfedgePtr connectVertices(FacePtr face, VertexPtr vA, VertexPtr vB);

  // Triangulate in a face, returns all subfaces
  std::vector<FacePtr> triangulate(FacePtr face);

  void permuteToCanonical(); // permute to the same indexing convention as after construction

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
  std::vector<std::vector<size_t>> getPolygonSoupFaces();
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

  size_t nextElemID = 77; // used to assign unique ID to elements

  // Hide copy and move constructors, we don't wanna mess with that
  HalfedgeMesh(const HalfedgeMesh& other) = delete;
  HalfedgeMesh& operator=(const HalfedgeMesh& other) = delete;
  HalfedgeMesh(HalfedgeMesh&& other) = delete;
  HalfedgeMesh& operator=(HalfedgeMesh&& other) = delete;

  // Used to resize the halfedge mesh. Expands and shifts vectors as necessary.
  Halfedge* getNewRealHalfedge();
  Halfedge* getNewImaginaryHalfedge();
  Vertex* getNewVertex();
  Edge* getNewEdge();
  Face* getNewFace();

  // Performs a sanity checks on halfedge structure; throws on fail
  void validateConnectivity();

  // Cache some basic information that may be queried many
  // times, but require O(n) computation to determine.
  // FIXME this is now broken with respect to modifications
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

  size_t indexOf(Halfedge* ptr);
  size_t indexOf(Vertex* ptr);
  size_t indexOf(Edge* ptr);
  size_t indexOf(Face* ptr);

  friend class DynamicHalfedgePtr;
  friend class DynamicVertexPtr;
  friend class DynamicEdgePtr;
  friend class DynamicFacePtr;
};

class Halfedge {
  friend class Edge;
  friend class HalfedgeMesh;
  friend class HalfedgePtr;
  friend struct std::hash<HalfedgePtr>;
  friend struct std::hash<DynamicHalfedgePtr>;

protected:
  Halfedge* twin;
  Halfedge* next;
  Vertex* vertex;
  Edge* edge;
  Face* face;

  bool isReal = true;
  size_t ID; // a unique value useful for hashing (etc). NOT an index

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
  friend struct std::hash<VertexPtr>;
  friend struct std::hash<DynamicVertexPtr>;

protected:
  // Data structure
  Halfedge* halfedge; // some halfedge that emanates from this vertex
                      // (guaranteed to be real)
  bool isBoundary = false;
  size_t ID; // a unique value useful for hashing (etc). NOT an index

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
  friend struct std::hash<EdgePtr>;
  friend struct std::hash<DynamicEdgePtr>;

protected:
  Halfedge* halfedge;

  bool flip(void);

  bool isBoundary = false;
  size_t ID; // a unique value useful for hashing (etc). NOT an index

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
  friend struct std::hash<FacePtr>;
  friend struct std::hash<DynamicFacePtr>;

protected:
  Halfedge* halfedge;

  bool isBoundary = false;
  bool isReal = false;
  size_t ID; // a unique value useful for hashing (etc). NOT an index

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
