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
  bool orientation;
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

  std::unique_ptr<ManifoldSurfaceMesh> mesh;        // optional: mesh connectivity
  std::unique_ptr<SimplePolygonMesh> simpleMesh;    // fallback if manifold mesh construction fails
  VertexData<CommonSubdivisionPoint*> sourcePoints; // for each vertex in mesh, what source point did it come from

  FaceData<Face> sourceFaceA;
  FaceData<Face> sourceFaceB;
  std::vector<Face> sourceFaceB_vec;

  // === Accessors
  // precondition: must lie along e
  // must not be a face point
  int getOrderAlongEdgeA(const CommonSubdivisionPoint& p, Edge eA);
  int getOrderAlongEdgeB(const CommonSubdivisionPoint& p, Edge eB);

  std::vector<SurfacePoint> getHalfedgePathAonB(Halfedge heA);
  std::vector<SurfacePoint> getHalfedgePathBonA(Halfedge heB);

  // Global index
  int getIndex(const CommonSubdivisionPoint& p);
  int getIndex(const CommonSubdivisionPoint* p);

  // Interpolate data from one of the meshes to the common subdivision
  template <typename T>
  VertexData<T> interpolateAcrossA(const VertexData<T>& dataA) const;
  template <typename T>
  VertexData<T> interpolateAcrossB(const VertexData<T>& dataB) const;

  // Interpolation matrices
  SparseMatrix<double> interpolationMatrixA();
  SparseMatrix<double> interpolationMatrixB();

  EdgeData<double> interpolateEdgeLengthsA(const EdgeData<double>& lengthA);
  EdgeData<double> interpolateEdgeLengthsB(const EdgeData<double>& lengthB);

  // Warning: triangulates common subdivision mesh
  SparseMatrix<double> vertexGalerkinMassMatrixFromPositionsA(const VertexData<Vector3>& positionsA);
  SparseMatrix<double> vertexGalerkinMassMatrixFromPositionsB(const VertexData<Vector3>& positionsB);
  SparseMatrix<double> vertexGalerkinMassMatrixFromLengthsA(const EdgeData<double>& lengthsA);
  SparseMatrix<double> vertexGalerkinMassMatrixFromLengthsB(const EdgeData<double>& lengthsB);

  // Number of elements in common subdivision
  // Warning: not constant time: requires mesh traversal to compute
  // TODO: assumes mesh B is finer than mesh A
  size_t nVertices() const;
  size_t nEdges() const;
  size_t nFaces() const;
  // compute (nVertices, nEdges, nFaces) all at once
  std::tuple<size_t, size_t, size_t> elementCounts() const;

  size_t intersectionsA(Edge eA) const; // count intersections along edge eA of meshA
  size_t intersectionsB(Edge eB) const; // count intersections along edge eB of meshB

  // === Mutators

  // TODO: right now, this assumes that mesh B is finer than mesh A. In the
  // future, that might not be the case and we'll have to think harder
  void constructMesh();

  void triangulateMesh();

  // Write mesh A and common subdivision to obj files
  // Vertex positions should be for mesh A
  void writeToFile(std::string filename, const VertexData<Vector3>& vertexPositions, int kColors = 7);
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
