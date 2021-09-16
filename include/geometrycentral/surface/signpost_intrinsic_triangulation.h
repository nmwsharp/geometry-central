#pragma once

#include "geometrycentral/surface/embedded_geometry_interface.h"
#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/intrinsic_triangulation.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/utilities/elementary_geometry.h"


// An implmementation of the Signpost datastructure from
//   > "Navigating Intrinsic Triangulations". Sharp, Soliman, and Crane. SIGGRAPH 2019


namespace geometrycentral {
namespace surface {


class SignpostIntrinsicTriangulation : public IntrinsicTriangulation {

public:
  // Construct an intrinsic triangulation which sits atop this input mesh. Initially, the input triangulation will
  // just be a copy of the input mesh.
  SignpostIntrinsicTriangulation(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& inputGeom);

  // ======================================================
  // ======== Core Members
  // ======================================================

  // Direction of each halfedge, in radians from [0, angleSum]. Importantly, this is different from
  // IntrinsicGeometryInterface's `halfedgeVectorsInVertex`, because it is an _angle only_ measured in radians, and it
  // is _not_ rescaled according to the vertex curvature.
  HalfedgeData<double> signpostAngle;

  // did this edge come from the original triangulation? used mainly for optimizations.
  EdgeData<bool> edgeIsOriginal;

  // ======================================================
  // ======== Queries & Accessors
  // ======================================================

  SurfacePoint equivalentPointOnIntrinsic(const SurfacePoint& pointOnInput) override;

  SurfacePoint equivalentPointOnInput(const SurfacePoint& pointOnIntrinsic) override;

  std::vector<SurfacePoint> traceIntrinsicHalfedgeAlongInput(Halfedge intrinsicHe) override;
  std::vector<SurfacePoint> traceIntrinsicHalfedgeAlongInput(Halfedge intrinsicHe, bool trimEnd);

  std::vector<SurfacePoint> traceInputHalfedgeAlongIntrinsic(Halfedge inputHe) override;
  std::vector<SurfacePoint> traceInputHalfedgeAlongIntrinsic(Halfedge inputHe, bool trimEnd);

  // ======================================================
  // ======== Low-Level Mutators
  // ======================================================

  // If the edge is not Delaunay, flip it. Returns true if flipped.
  bool flipEdgeIfNotDelaunay(Edge e) override;

  // If the edge can be flipped, flip it (must be combinatorially flippable and inside a convex quad). Returns true if
  // flipped.
  bool flipEdgeIfPossible(Edge e) override;

  // Flip an edge, where the caller specifies geometric data for the updated edge, rather than it being computed. Must
  // be flippable. Experts only.
  void flipEdgeManual(Edge e, double newLength, double forwardAngle, double reverseAngle, bool isOrig,
                      bool reverseFlip = false);

  // Insert a new vertex in to the intrinsic triangulation
  Vertex insertVertex(SurfacePoint newPositionOnIntrinsic) override;

  // Remove an (inserted) vertex from the triangulation.
  // Note: if something goes terribly (numerically?) wrong, will exit without removing the vertex.
  Face removeInsertedVertex(Vertex v) override;

  // Split an edge
  Halfedge splitEdge(Halfedge he, double tSplit) override;

protected:
private:
  // ======================================================
  // ======== Geometry Interface
  // ======================================================
  //
  // Satisfy the requirements of the IntrinsicGeometryInterface

  // Override the halfedge vectors method from intrinsic geometry
  virtual void computeHalfedgeVectorsInVertex() override;


  // ======================================================
  // ======== Helpers
  // ======================================================

  // Insertion helpers
  Vertex insertVertex_face(SurfacePoint newPositionOnIntrinsic);
  Halfedge insertVertex_edge(SurfacePoint newPositionOnIntrinsic);
  void resolveNewVertex(Vertex newV, SurfacePoint intrinsicPoint);

  // Update a signpost angle from the (counter-)clockwise neighboring angle
  void updateAngleFromCWNeighor(Halfedge he);

  // Map angle to range [0, angleSum)
  double standardizeAngle(Vertex vert, double angle) const;

  // Get vector in rescaled vertex coordinates (either halfedge or arbitrary vector)
  Vector2 halfedgeVector(Halfedge he) const;
  Vector2 rescaledVertexVector(Vertex v, double angle, double len) const;

  // Scale factor to take Euclidean data to cone data
  double vertexAngleScaling(Vertex v) const;

  // Recover t-values after tracing
  // Note that really we ought to just report these back from the tracing routine itself, which computes them internally. We don't have a nice API for passing that data around, so this lazily recovers it
  std::vector<double> recoverTraceTValues(const std::vector<SurfacePoint>& edgeTrace);

  // Construct the common subdivision for the current triangulation.
  // WARNING: this signpost implementation does not properly populate some fields
  // of the common subdivision, in particular the `orientation` field for intersections.
  void constructCommonSubdivision() override;
};


} // namespace surface
} // namespace geometrycentral


#include "geometrycentral/surface/signpost_intrinsic_triangulation.ipp"
