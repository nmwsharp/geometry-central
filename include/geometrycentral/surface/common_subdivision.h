#pragma once

#include "geometrycentral/surface/barycentric_coordinate_helpers.h"
#include "geometrycentral/surface/edge_length_geometry.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/normal_coordinates.h"
#include "geometrycentral/surface/simple_polygon_mesh.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/utilities/elementary_geometry.h"

#include <deque>

namespace geometrycentral {
namespace surface {

// EDGE_PARALLEL intersections are kind of a hack. In the case of shared edges,
// we need to store which edge of mesh A corresponds to which edge of mesh B,
// since edges may not be determined by their endpoints. To do so, we store the
// two edges in an EDGE_PARALLEL intersection
// TODO: add in VERTEX_FACE, FACE_VERTEX, VERTEX_EDGE, etc
enum class CSIntersectionType {
  VERTEX_VERTEX,
  EDGE_TRANSVERSE,
  EDGE_PARALLEL,
  FACE_VERTEX, // Face of mesh A, Vertex of mesh B
  EDGE_VERTEX  // Edge of mesh A, Vertex of mesh B
};

// on the common subdivision between meshes "A" and "B"
struct CommonSubdivisionPoint {
  CSIntersectionType intersectionType;
  SurfacePoint posA;
  SurfacePoint posB;
  bool orientation; // WARNING: some constructors (e.g. signposts) do not properly populate this
};
std::ostream& operator<<(std::ostream& out, const CSIntersectionType& type);
std::ostream& operator<<(std::ostream& out, const CommonSubdivisionPoint& pt);

// HACK: Represent shared edge paths using [start vertex, shared halfedge,
// end vertex] where the middle subdivision point is of type EDGE_PARALLEL
// This is necessary because there could be multiple halfedges in either mesh
// from start vertex to end vertex, and we need to disambiguate somehow
class CommonSubdivision {
public:
  CommonSubdivision(ManifoldSurfaceMesh& meshA_, ManifoldSurfaceMesh& meshB_);

  ManifoldSurfaceMesh& meshA;
  ManifoldSurfaceMesh& meshB;

  // Use std::deque instead of std::vector because std::deque will not
  // reallocate when items are added, meaning that pointers to elements of
  // std::deque stay valid
  // If we used a vector, pointers into it would get invalidated when the
  // vector changed size
  std::deque<CommonSubdivisionPoint> subdivisionPoints;

  // INCLUDES start and end point
  EdgeData<std::vector<CommonSubdivisionPoint*>> pointsAlongA;
  EdgeData<std::vector<CommonSubdivisionPoint*>> pointsAlongB;

  // === Accessors

  // For any halfedge of A, return its sequence of points on B (and vice versa)
  // Includes endpoints.
  std::vector<SurfacePoint> getHalfedgePathAonB(Halfedge heA);
  std::vector<SurfacePoint> getHalfedgePathBonA(Halfedge heB);
  // Number of elements in common subdivision
  // Warning: not constant time: requires mesh traversal to compute
  size_t nVertices() const;
  size_t nEdges() const;
  size_t nFaces() const;
  // compute (nVertices, nEdges, nFaces) all at once
  std::tuple<size_t, size_t, size_t> elementCounts() const;

  size_t intersectionsA(Edge eA) const; // count intersections along edge eA of meshA
  size_t intersectionsB(Edge eB) const; // count intersections along edge eB of meshB

  // Construct and return a SimplePolygonMesh of the common subdivision.
  // In cases where there are some errors in intersection data defining the common subvidision, it may not be possible
  // to construct the `ManifoldSurfaceMesh` with `constructMesh()`, and this is the only option.
  std::unique_ptr<SimplePolygonMesh> buildSimpleMesh();

  // ===============================================
  // ==== Common subdivision mesh connectivity
  // ===============================================
  //
  // Optionally, the raw intersections stored in the common subdivision can be used to explicitly construct a
  // SurfaceMesh (with all the usual halfedge connectivity) of the common subdivision; this mesh can then be used for
  // many higher-level geometric operations.
  //
  // **The routines and members in this section all require that constructMesh() has been called first.**
  //
  // Note that in cases where there are some errors in intersection data defining the common subdivision, it may not be
  // possible to construct a manifold mesh of the common subdivision. `constructMesh()` will fail with an exception in
  // this case. Note that `buildSimpleMesh()` above can be used to build a plain old vertex-face adjacency list
  // representation of the mesh, although it cannot be used for the various routines in this section.

  // === Members

  // The mesh connectivity for the common subdivision.
  std::unique_ptr<ManifoldSurfaceMesh> mesh;

  // For each vertex in mesh, what source point did it come from
  VertexData<CommonSubdivisionPoint*> sourcePoints;

  // For each face in the common subdivision, which face in meshA is it a sub-face of?
  FaceData<Face> sourceFaceA;

  // For each face in the common subdivision, which face in meshB is it a sub-face of?
  FaceData<Face> sourceFaceB;


  // === Methods

  // Construct `mesh` and auxiliary data. Throws on failure (see note above)
  void constructMesh(bool triangulate = true, bool skipIfAlreadyConstructed = true);
  void triangulateMesh();

  // Interpolate data at vertices from one of the meshes to the common subdivision. The return value is defined
  // per-vertex of the commonm subdivision mesh.
  template <typename T>
  VertexData<T> interpolateAcrossA(const VertexData<T>& dataA) const;
  template <typename T>
  VertexData<T> interpolateAcrossB(const VertexData<T>& dataB) const;

  // Copy data at faces from one of the meshes to the common subdivision. Each face of the common subdivision gets
  // the value from the face which contains it. The return value is defined per-face of the common subdivision mesh.
  template <typename T>
  FaceData<T> copyFromA(const FaceData<T>& dataA) const;
  template <typename T>
  FaceData<T> copyFromB(const FaceData<T>& dataB) const;

  // Interpolation matrices
  // Gives a |V_c| x |V_a| sparse matrix (where |V_c| is number of vertice subdivision) such that
  // multplying by a vector of values at the vertices of A gives interpolated values at the vertices of the subdivision.
  SparseMatrix<double> interpolationMatrixA();
  SparseMatrix<double> interpolationMatrixB(); // And respectively for B.

  // Use edge lengths on either of the source triangulations to get edge lengths
  // for the common subdivision.
  // Note that in the standard case of an intrinsic triangulation with Euclidean metric-preserving edge flips, calling
  // either of these methods with the edge lengths from the respective triangulation will produce identical outputs (up
  // to floating-point error).
  EdgeData<double> interpolateEdgeLengthsA(const EdgeData<double>& lengthA);
  EdgeData<double> interpolateEdgeLengthsB(const EdgeData<double>& lengthB);

  // Warning: triangulates common subdivision mesh
  SparseMatrix<double> vertexGalerkinMassMatrixFromPositionsA(const VertexData<Vector3>& positionsA);
  SparseMatrix<double> vertexGalerkinMassMatrixFromPositionsB(const VertexData<Vector3>& positionsB);
  SparseMatrix<double> vertexGalerkinMassMatrixFromLengthsA(const EdgeData<double>& lengthsA);
  SparseMatrix<double> vertexGalerkinMassMatrixFromLengthsB(const EdgeData<double>& lengthsB);


  // Write mesh A and common subdivision to obj files
  // Vertex positions should be for mesh A
  void writeToFile(std::string filename, const VertexData<Vector3>& vertexPositions, int kColors = 7);


  // ===============================================
  // ==== Debugging routines, etc
  // ===============================================

  // Global index
  // Exhaustive search, very expensive!
  int getIndex(const CommonSubdivisionPoint& p);
  int getIndex(const CommonSubdivisionPoint* p);

  // If p is the i'th point along edge e, returns i
  // Uses an simple search along the edge; generally this function is not needed.
  // precondition: must lie along e, must not be a face point
  int getOrderAlongEdgeA(const CommonSubdivisionPoint& p, Edge eA);
  int getOrderAlongEdgeB(const CommonSubdivisionPoint& p, Edge eB);

  // Helper routine abstracting common logic for mesh creation.
  void constructMeshData(std::vector<std::vector<size_t>>& faces_out, std::vector<CommonSubdivisionPoint*>& parents_out,
                         std::vector<Face>& sourceFaceA_out, std::vector<Face>& sourceFaceB_out);

  // Throws an error if the mesh has not yet been constructed.
  void checkMeshConstructed() const;
};

std::vector<std::vector<size_t>> sliceFace(const std::vector<size_t>& pij, const std::vector<size_t>& pjk,
                                           const std::vector<size_t>& pki);

// Precondition: pij.size() >= pjk.size(), pki.size()
std::vector<std::vector<size_t>> sliceNicelyOrderedFace(const std::vector<size_t>& pij, const std::vector<size_t>& pjk,
                                                        const std::vector<size_t>& pki);

// Take evenly spaced values on [0,1] and treat them as categorical labels. These labels are assigned to mesh faces in
// the sense of a graph coloring, with further heuristics to try to avoid neighbors-of-neighbors at vertices.
FaceData<double> niceColors(ManifoldSurfaceMesh& mesh, int kColors = 7);

} // namespace surface
} // namespace geometrycentral

#include "geometrycentral/surface/common_subdivision.ipp"
