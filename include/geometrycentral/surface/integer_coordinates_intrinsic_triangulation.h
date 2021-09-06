#pragma once

#include "geometrycentral/surface/barycentric_coordinate_helpers.h"
#include "geometrycentral/surface/common_subdivision.h"
#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/intrinsic_mollification.h"
#include "geometrycentral/surface/intrinsic_triangulation.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/normal_coordinates.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/surface/trace_geodesic.h"
#include "geometrycentral/utilities/elementary_geometry.h"

#include <deque>

namespace geometrycentral {
namespace surface {

class IntegerCoordinatesIntrinsicTriangulation : public IntrinsicTriangulation {

public:
  // Construct an intrinsic triangulation which sits atop this input mesh.
  // Initially, the input triangulation will just be a copy of the input mesh.
  IntegerCoordinatesIntrinsicTriangulation(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& inputGeom,
                                           double mollifyEPS = 1e-5);

  // ======================================================
  //                   Core Members
  // ======================================================

  // The actual normal coordinates (and roundabouts) encoding the triangulation. These normal coordinates are defined on
  // top of the _intrinsic_ mesh---for each intrinsic edge, the encode how many original edges cross it.
  NormalCoordinates normalCoordinates;

  // ======================================================
  // ======== Queries & Accessors
  // ======================================================

  // See intrinsic_triangulation.h for docs.

  SurfacePoint equivalentPointOnIntrinsic(const SurfacePoint& pointOnInput) override;

  SurfacePoint equivalentPointOnInput(const SurfacePoint& pointOnIntrinsic) override;

  EdgeData<std::vector<SurfacePoint>> traceAllIntrinsicEdgesAlongInput() override;
  std::vector<SurfacePoint> traceIntrinsicHalfedgeAlongInput(Halfedge intrinsicHe) override;

  EdgeData<std::vector<SurfacePoint>> traceAllInputEdgesAlongIntrinsic() override;
  std::vector<SurfacePoint> traceInputHalfedgeAlongIntrinsic(Halfedge inputHe) override;

  // ======================================================
  // ======== Low-Level Mutators
  // ======================================================

  // See intrinsic_triangulation.h for docs.

  bool flipEdgeIfNotDelaunay(Edge e) override;

  bool flipEdgeIfPossible(Edge e) override;

  Vertex insertVertex(SurfacePoint newPositionOnIntrinsic) override;

  Face removeInsertedVertex(Vertex v) override;

  Halfedge splitEdge(Halfedge he, double tSplit) override;

  // Check if an edge can be flipped geometrically, as defined by the (relative) signed areas of the resulting triangles; positive values mean flippable.
  double checkFlip(Edge e);

  // Insert circumcenter or split segment
  Vertex insertCircumcenterOrSplitSegment(Face f, bool verbose = false);

  Vertex splitFace(Face f, Vector3 bary, bool verbose = false);
  Vertex splitEdge(Edge e, double bary, bool verbose = false);
  Vertex splitInteriorEdge(Edge e, double bary, bool verbose = false);
  Vertex splitBoundaryEdge(Edge e, double bary, bool verbose = false);

  // Move a vertex `v` in direction `vec`, represented as a vector in the
  // vertex's tangent space.
  Vertex moveVertex(Vertex v, Vector2 vec);

  // Assumes intrinsicEdgeLengths is up to date
  void updateCornerAngle(Corner c);

  // Assumes cornerAngles, vertexAngleSums exist and are up to date
  void updateHalfedgeVectorsInVertex(Vertex v);

  // ======================================================
  //                Low-Level Queries
  // ======================================================

  // Takes in a halfedge of the intrinsic mesh whose edge's normal coordinate
  // is negative (meaning that it lies along an edge of the input mesh) and
  // returns the halfedge in the input mesh pointing in the same direction
  // e.vertex() must live in both meshes
  Halfedge getSharedInputEdge(Halfedge e) const;

  // Takes in an intrinsic point, represented as an intrinsic face and barycentric coordinate,
  // and computes the corresponding point on the input mesh, as well as the normal coordinates of
  // the edges connecting this new point to f's vertices.
  std::pair<SurfacePoint, std::array<int, 3>> computeFaceSplitData(Face f, Vector3 bary, bool verbose = false);

  // Compute the number of vertices in the common subdivision
  // i.e. intrinsicMesh->nVertices() plus the sum of the normal coordinates
  size_t nSubdividedVertices() const;

  // HACK: represents arcs parallel to a mesh edge with a single pair {-n,
  // he} where n is the number of arcs parallel to he.edge() Trace an edge
  // of the input mesh over the intrinsic triangulation
  NormalCoordinatesCompoundCurve traceInputEdge(Edge e, bool verbose = false) const;
  NormalCoordinatesCompoundCurve traceInputHalfedge(Halfedge inputHe, bool verbose = false) const;

  std::pair<bool, NormalCoordinatesCurve> traceNextCurve(const NormalCoordinatesCurve& oldCurve,
                                                         bool verbose = false) const;

  // Inverse to traceInputEdge
  Halfedge identifyInputEdge(const NormalCoordinatesCurve& path, bool verbose = false) const;

  // Identify shared halfedge, throw exception if halfedge is not shared
  // (i.e. edgeCoords[he.edge()] must be negative)
  Halfedge identifyInputEdge(Halfedge he) const;

  std::array<Vector2, 3> vertexCoordinatesInFace(Face face) const;

  void setFixedEdges(const EdgeData<bool>& fixedEdges);

  // If f is entirely contained in some face of the input mesh, return that
  // face Otherwise return Face()
  Face getParentFace(Face f) const;

private:
  // Implementation details

  // Construct the common subdivision for the current triangulation.
  void constructCommonSubdivision() override;
};

// Compute the cotan weight of halfedge ij in terms of the lengths of its
// neighbors
double halfedgeCotanWeight(double lij, double ljk, double lki);

FaceData<Vector2> interpolateTangentVectorsB(const IntegerCoordinatesIntrinsicTriangulation& tri,
                                             const CommonSubdivision& cs, const FaceData<Vector2>& dataB);


} // namespace surface
} // namespace geometrycentral

#include "integer_coordinates_intrinsic_triangulation.ipp"
