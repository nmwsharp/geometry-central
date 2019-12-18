#pragma once

#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/utilities/elementary_geometry.h"


// An implmementation of the Signpost datastructure from
//   > "Navigating Intrinsic Triangulations". Sharp, Soliman, and Crane. SIGGRAPH 2019


namespace geometrycentral {
namespace surface {


class SignpostIntrinsicTriangulation : public IntrinsicGeometryInterface {

public:
  // Construct an intrinsic triangulation which sits atop this input mesh. Initially, the input triangulation will
  // just be a copy of the input mesh.
  SignpostIntrinsicTriangulation(IntrinsicGeometryInterface& inputGeom);

  // ======================================================
  // ======== Core Members
  // ======================================================

  // The underlying surface on which the intrinsic triangulation has been constructed
  HalfedgeMesh& inputMesh;
  IntrinsicGeometryInterface& inputGeom;

  // The connectivity of the intrinsic triangulation
  // note that somewhat confusingly, there is a .mesh reference which points to this same mesh,
  // inherited from the geometry interface
  std::unique_ptr<HalfedgeMesh> intrinsicMesh;

  // The geometry of the intrinsic triangulation
  EdgeData<double> intrinsicEdgeLengths;            // length of each edge
  HalfedgeData<double> intrinsicHalfedgeDirections; // direction of each halfedge, in radians from [0, angleSum]
  VertexData<double> intrinsicVertexAngleSums;      // vertex cone angle sum
  EdgeData<char> edgeIsOriginal; // did this edge come from the original triangulation? used mainly for optimizations.

  // NOTE: To enable use to make efficient use of the surface tracers, this class always automatically updates the
  // halfedgeVectorsInVertex and halfedgeVectorsInFace geometry members. Could remove this requirement if we change the
  // way the tracer works.

  // Vertex locations for the intrinsic triangulation
  VertexData<SurfacePoint> vertexLocations;

  // Minor parameters
  double delaunayEPS = 1e-6; // epsilon to use when testing if an edge is Delaunay

  // ======================================================
  // ======== Queries & Accessors
  // ======================================================

  // Given a point on the input triangulation, returns the corresponding point on the intrinsic triangulation
  SurfacePoint equivalentPointOnIntrinsic(const SurfacePoint& pointOnInput);

  // Given a point on the intrinsic triangulation, returns the corresponding point on the input triangulation
  SurfacePoint equivalentPointOnInput(const SurfacePoint& pointOnIntrinsic);

  // Trace out the edges of the intrinsic triangulation along the surface of the input mesh.
  // Each path is ordered along edge.halfedge(), and includes both the start and end points
  EdgeData<std::vector<SurfacePoint>> traceEdges();
  std::vector<SurfacePoint> traceHalfedge(Halfedge he); // trace a single intrinsic halfedge

  // Given data defined on the intrinsic triangulation, samples it at the vertices of the input triangulation
  template <typename T>
  VertexData<T> sampleAtInput(const VertexData<T>& dataOnIntrinsic);

  // Returns true if the intrinsic triangulation (or edge) satisifies the intrinsic Delaunay criterion
  bool isDelaunay();
  bool isDelaunay(Edge e);

  // Returns the smallest angle in the intrinsic triangulation, in degrees
  double minAngleDegrees();

  // ======================================================
  // ======== High-Level Mutators
  // ======================================================
  //
  // Call once to build a useful triangulation

  // Flips edge in the intrinsic triangulation until is satisfies teh intrinsic Delaunay criterion
  void flipToDelaunay();

  // Perform intrinsic Delaunay refinement the intrinsic triangulation until it simultaneously:
  //   - satisfies the intrinsic Delaunay criterion
  //   - has no angles smaller than `angleThreshDegrees` (values > 30 degrees may not terminate)
  //   - has no triangles larger than `circumradiusThresh`
  // Terminates no matter what after maxInsertions insertions (infinite by default)
  void delaunayRefine(double angleThreshDegrees = 25.,
                      double circumradiusThresh = std::numeric_limits<double>::infinity(),
                      size_t maxInsertions = INVALID_IND);


  // General version of intrinsic Delaunay refinement, taking a function which will be called
  // to determine if a triangle should be refined.
  // Will return only when all triangles pass this function, or maxInsertions is exceeded, so
  // be sure to chose arguments such that the function terminates.
  void delaunayRefine(const std::function<bool(Face)>& shouldRefine, size_t maxInsertions = INVALID_IND);


  // ======================================================
  // ======== Low-Level Mutators
  // ======================================================
  //
  // Basic operations to modify the intrinsic triangulation

  // If the edge is not Delaunay, flip it. Returns true if flipped.
  bool flipEdgeIfNotDelaunay(Edge e);

  // If the edge can be flipped, flip it (must be combinatorially flippable and inside a convex quad). Returns true if
  // flipped.
  bool flipEdgeIfPossible(Edge e, double possibleEPS = 1e-6);

  // Insert a new vertex in to the intrinsic triangulation
  Vertex insertVertex(SurfacePoint newPositionOnIntrinsic);

  // Insert the circumcenter of a face in to the triangulation. Returns the newly created intrinsic vertex.
  Vertex insertCircumcenter(Face f);

  // Insert the barycenter of a face in to the triangulation. Returns the newly created intrinsic vertex.
  Vertex insertBarycenter(Face f);

  // ======================================================
  // ======== Geometry Immediates
  // ======================================================

  // Computed on the intrinsic triangulation
  double cornerAngle(Corner c) const;
  double halfedgeCotanWeight(Halfedge he) const;
  double edgeCotanWeight(Edge e) const;
  double area(Face f) const;
  double shortestEdge(Face f) const;
  double circumradius(Face f) const;


private:
  // ======================================================
  // ======== Geometry Interface
  // ======================================================
  //
  // Satisfy the requirements of the IntrinsicGeometryInterface

  // Override the compute edge lengths method from intrinsic geometry.
  virtual void computeEdgeLengths() override;

  // Override the halfedge vectors method from intrinsic geometry
  virtual void computeHalfedgeVectorsInVertex() override;


  // ======================================================
  // ======== Helpers
  // ======================================================

  // Insertion helpers
  Vertex insertVertex_face(SurfacePoint newPositionOnIntrinsic);
  Vertex insertVertex_edge(SurfacePoint newPositionOnIntrinsic);
  void resolveNewVertex(Vertex newV);

  // Update a signpost angle from the clockwise neighboring angle
  void updateAngleFromCWNeighor(Halfedge he);

  // Map angle to range [0, angleSum)
  double standardizeAngle(Vertex vert, double angle) const;

  // Get vector in rescaled vertex coordinates
  Vector2 halfedgeVector(Halfedge he) const;

  // Scale factor to take Euclidean data to cone data
  double vertexAngleScaling(Vertex v) const;

  // Repopulate the member halfedgeVectorInFace
  void updateFaceBasis(Face f);

  // Isometrically lay out the vertices around a halfedge in 2D coordinates
  // he points from vertex 2 to 0; others are numbered CCW
  std::array<Vector2, 4> layoutDiamond(Halfedge he);
  std::array<Vector2, 3> vertexCoordinatesInTriangle(Face face);

  // Helper for layoutDiamond()
  static Vector2 layoutTriangleVertex(const Vector2& pA, const Vector2& pB, const double& lBC, const double& lCA);
};


} // namespace surface
} // namespace geometrycentral


#include "geometrycentral/surface/signpost_intrinsic_triangulation.ipp"

