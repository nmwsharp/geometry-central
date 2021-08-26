#pragma once

#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/utilities/elementary_geometry.h"


namespace geometrycentral {
namespace surface {

// Halfedges always point "left"
struct NormalCoordinatesCurve {
  std::vector<std::pair<int, Halfedge>> crossings;
};
struct NormalCoordinatesCompoundCurve {
  std::vector<NormalCoordinatesCurve> components;
};


// A set of normal coordinates atop a triangulation. Represents general
// non-crossing curves which may terminate at vertices.
class NormalCoordinates {
public:
  // == Constructors

  // Initializes an empty set of normal coords
  NormalCoordinates(ManifoldSurfaceMesh& mesh);

  // == Accessors
  int operator[](Edge e) const;   // gives normal coodinates across the edge
  int& operator[](Edge e);        // set the normal coordinates across the edge
  int operator[](Corner c) const; // gives corner coordinates implied by edge coordinates

  // == Members

  // The surface atop which the coordinates are defined
  ManifoldSurfaceMesh& mesh;

  // Normal coordinates per edge
  EdgeData<int> edgeCoords;

  // In counterclockwise order. If a vertex is on the boundary, the first
  // halfedge must be the one pointing along the boundary in the interior of
  // the mesh, and the last one must be the one pointing along the boundary in
  // the exterior of the mesh
  // As always, assumes vertices are manifold
  HalfedgeData<int> roundabouts;

  // Count of edges emanating from each corner (nonnegative)
  VertexData<int> roundaboutDegrees;
  
  // === Initialization
  void setCurvesFromEdges(ManifoldSurfaceMesh& mesh);

  // === General routines

  // Call immediately before mesh->flip(e);
  std::tuple<int, size_t, size_t> computeFlippedData(Edge e);

  // Call after mesh->flip(e) to update normal coordinates
  void applyFlippedData(Edge e, const std::tuple<int, size_t, size_t>& update);


  // === Mutation

  // Vertex Insertion
  // Call immediately before mesh->insertVertex(f);
  // Input are the desired normal coordinates of the edges connecting face f's
  // vertices to the new vertex, in the order given by f.adjacentVertices()
  // TODO: validate that these inputs are possible
  std::array<int, 3> computeVertexInsertionData(Face f, const std::array<int, 3>& newCrossingCounts);
  std::array<int, 3> computeVertexInsertionDataGeodesic(IntrinsicGeometryInterface& geo, Face f,
                                                        const Vector3 location);

  // Call after mesh->insertVertex(f) to update normal coordinates;
  void applyVertexInsertionData(Vertex newVertex, const std::array<int, 3>& update);

  std::array<int, 4> computeInteriorEdgeSplitDataGeodesic(IntrinsicGeometryInterface& geo, Edge e, double location);

  std::array<int, 3> computeBoundaryEdgeSplitDataGeodesic(IntrinsicGeometryInterface& geo, Edge e, double location);


  // Check that these are valid normal coordinates, throw if not
  void validate() const;

  // *** CURVES ARE ALWAYS ORIENTED SUCH THAT WHEN YOU WALK FORWARD ALONG THE
  // CURVE, HALFEDGE POINT TO THE LEFT ***

  // Given a curve identified by point where it starts, topological trace
  // it out. Produces a list of halfedges which cross the curve from right to
  // left
  // As usual, index is 0-indexed
  NormalCoordinatesCurve topologicalTrace(Halfedge he, int index) const;
  // As usual, index is 0-indexed
  NormalCoordinatesCurve topologicalTrace(Corner c, int index) const;

  // Traces forward and back
  // return value is a tuple of the list of all crossings, and an index in to
  // the crossing list indicating which was the crossing you gave as input
  // Returns a curve that intersects he w/ positive orientation (i.e. he is in
  // the list of halfedges crossed by the curve, rather than he.twin()) Recall
  // also that curves always intersect left-pointing halfedges
  std::tuple<NormalCoordinatesCurve, int> topologicalTraceBidirectional(Halfedge he, int index) const;

  // HACK: represents arcs parallel to a mesh edge with a single pair {-n, he}
  // where n is the number of arcs parallel to he.edge()
  // TODO: Doesn't handle boundary properly
  std::vector<NormalCoordinatesCurve> topologicalTraceAllCurves() const;

  // Gives ANY noncrossing geometry which is topologically correct
  // (equally spaced)
  std::vector<std::vector<SurfacePoint>> generateAnyGeometry() const;


  // === Geodesic routines
  std::vector<std::vector<SurfacePoint>> generateGeodesicGeometry(IntrinsicGeometryInterface& geo) const;

  // Returns the t-value for the crossing of the ind'th curve along this edge
  // (under the assumption that it is geodesic)
  double generateGeodesicCrossingLocation(IntrinsicGeometryInterface& geo, Halfedge he, int ind) const;

  // Returns the ordered t-values for the crossing of all curves along this
  // edge (under the assumption that they are geodesic)
  std::vector<double> generateGeodesicCrossingLocations(IntrinsicGeometryInterface& geo, Halfedge he) const;

  // === Helpers
  // Compute a corner's coordinate from the edge coordinates in its triangle
  // TODO: this counts twice the number of arcs that cross this corner, since
  // this computation also gives an integer answer when arcs emanate from the
  // corner. Is this the right tradeoff?
  int cornerCoord(Corner c) const;

  // Counts how many arcs come out of this corner and exit at the opposite
  // edge (i.e. ignores arcs parallel to triangle edges)
  size_t strictDegree(Corner c) const;

  // Counts how many arcs clip off this corner, ignoring arcs that emanate
  size_t strictCornerCoord(Corner c) const;

  // He is the twin of the halfedge in a face that the curve starts from,
  // index is the index of the curve, with 0 being closest to the source
  // vertex of he.
  // If the curve terminates at a vertex, returns true
  // Otherwise returns false and udpates he and index to be the next crossing
  // As usual, index is 0-indexed
  bool stepTopologicalCurve(Halfedge& he, int& index) const;

  // Checks for a violation of the normal coordinate triangle inequality in
  // face f
  // If there is one, returns true, and sets violatingHe to the long halfedge
  // Otherwise returns false
  bool triangleInequalityViolation(Face f, Halfedge& violatingHe) const;

  // at least one loop wraps directly around the vertex
  bool isEncircledByLoopCurve(Vertex v) const;

  // loop in a "rubber band" configuration is hooked around the vertex,
  // entering and exiting the 1-ring through the same edge
  bool isHookedByCurve(Vertex v) const;

  // Set roundabout relative to roundabout for previous (in counterclockwise
  // order) halfedge
  void setRoundaboutFromPrevRoundabout(Halfedge he);
};

// == Helper function
// They're easier to test if they're outside of the class
int flipNormalCoordinates(int nij, int njk, int nki, int nil, int nlj);

// See NormalCoordinates::flipEdgeUpdate for an explanation of all these
// variable names;
std::pair<size_t, size_t> flipRoundabouts(int nij, int njk, int nki, int nil, int nlj, int nkl, size_t rki, size_t rlj,
                                          size_t dk, size_t dl);

// Count the number of arcs strictly leaving corner k (i.e. not parallel to jk
// or ki)
size_t strictDegree(int nij, int njk, int nki);

// Counts the number of arcs clipping off corner k. If arcs leave corner k, then
// they get counted as -1/2 each
int cornerCoord(int nij, int njk, int nki);

// Counts the number of arcs clipping off corner k, returns 0 if arcs leave
// corner k
size_t strictCornerCoord(int nij, int njk, int nki);

// == Geodesic Helpers
// Given trace info (like the output from topotrace), generate geodesic geometry.
// Barycentric coordinates are 0 at the src of and edge and 1 at the dst
std::vector<std::vector<SurfacePoint>> generateGeodesicGeometry(ManifoldSurfaceMesh& mesh,
                                                                IntrinsicGeometryInterface& geo,
                                                                const std::vector<NormalCoordinatesCurve>& traceCounts);


std::vector<SurfacePoint> generateSingleGeodesicGeometry(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geo,
                                                         const NormalCoordinatesCurve& curve);

// Record barycentric coordinate of each SurfacePoint along the curve
// Barycentric coordinates are 0 at the src of and edge and 1 at the dst
std::vector<std::pair<SurfacePoint, double>> generateFullSingleGeodesicGeometry(ManifoldSurfaceMesh& mesh,
                                                                                IntrinsicGeometryInterface& geo,
                                                                                const NormalCoordinatesCurve& curve);

// Compute the new normal coordinates for an inserted vertex v in face f given
// that v's barycentric coordinates and the barycentric coordinates of the
// crossings along f's edges
std::array<int, 3> computeVertexInsertionCrossingCounts(Vector3 bary,
                                                        const std::array<std::vector<double>, 3>& boundaryCrossings);

} // namespace surface
} // namespace geometrycentral
