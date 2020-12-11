#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/signpost_intrinsic_triangulation.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include <list>
#include <memory>
#include <queue>
#include <unordered_map>
#include <vector>


// Implements the edge-flip based algorithm for finding geodesics on triangular meshes, from
//      You Can Find Geodesic Paths in Triangle Meshes by Just Flipping Edges
//      Nicholas Sharp and Keenan Crane
//      ACM Trans. on Graph. (SIGGRAPH Asia 2020)
//
// The FlipEdgeNetwork class contains the high-level functionality.

namespace geometrycentral {
namespace surface {

// Forward declare
class FlipEdgeNetwork;

using HalfedgeListIter = std::list<Halfedge>::iterator;
using SegmentID = size_t;
enum class SegmentAngleType { Shortest = 0, LeftTurn, RightTurn }; // what direction does the min wedge point in

// Despite the name, this class represents
//
// - open paths
// - closed loops
// - path networks

// Simple mini class representing a path.
// Tightly coupled with FlipEdgeNetwork, where most of the functionality is implemented
//
// Each entry is
class FlipEdgePath {

public:
  FlipEdgePath(FlipEdgeNetwork& network, std::vector<Halfedge> halfedges, bool isClosed);

  // Members
  FlipEdgeNetwork& network; // the parent network of which this path is a part
  bool isClosed;

  // Data about the segments in path. For a path segment id, contains the corresponding halfedge, as well as the ID of
  // the previous/next segment (or INVALID_IND if they are path endpoints)
  std::unordered_map<SegmentID, std::tuple<Halfedge, SegmentID, SegmentID>>
      pathHeInfo; // ID --> <Halfedge, prevID, nextID>

  // Take path segment formed by the two halfedges prev(heNextIter) and heNextIter and replace it with the list of
  // halfedges in newHalfedges. Does all relevant bookkeeping, both here and in the FlipEdgeNetwork
  void replacePathSegment(SegmentID halfedgeNextID, SegmentAngleType angleType,
                          const std::vector<Halfedge>& newHalfedges);


  // Helpers to access the path in a more convenient format
  size_t size();
  std::vector<Halfedge> getHalfedgeList();
};


// A reference to a point along a path, which which tells us both the path of interest and the ID of the segment in that
// path.
struct FlipPathSegment {
  FlipEdgePath* path;
  SegmentID id;

  bool isEndpoint();
  FlipPathSegment next();
  FlipPathSegment prev();
  Halfedge halfedge();

  // This segment gets re-used as rear, this->next() is a new segment
  void splitEdge(double tSplit);

  // Inherit operators from tuple
  bool operator==(const FlipPathSegment& other) const;
  bool operator!=(const FlipPathSegment& other) const;
  bool operator>(const FlipPathSegment& other) const;
  bool operator>=(const FlipPathSegment& other) const;
  bool operator<(const FlipPathSegment& other) const;
  bool operator<=(const FlipPathSegment& other) const;
};

// A collection of paths. Contains all the high-level routines for interacting with path and applying algorithms.
class FlipEdgeNetwork {

public:
  // === Constructors

  // Construct a network from a collection of paths
  FlipEdgeNetwork(ManifoldSurfaceMesh& mesh_, IntrinsicGeometryInterface& inputGeom,
                  std::vector<std::vector<Halfedge>> paths, VertexData<bool> extraMarkedVerts = VertexData<bool>());

  // === Static initializers

  // Run Dijkstra between endpoints to initialize path
  static std::unique_ptr<FlipEdgeNetwork> constructFromDijkstraPath(ManifoldSurfaceMesh& mesh,
                                                                    IntrinsicGeometryInterface& geom, Vertex startVert,
                                                                    Vertex endVert);
  // Run Dijkstra between i'th and (i+1)'th point to initialize path
  static std::unique_ptr<FlipEdgeNetwork>
  constructFromPiecewiseDijkstraPath(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom,
                                     std::vector<Vertex> points, bool closed = false, bool markInterior = false);

  // Consturct path(s) from marked edges, heuristically inferring endpoints, loopiness, etc
  static std::unique_ptr<FlipEdgeNetwork> constructFromEdgeSet(ManifoldSurfaceMesh& mesh,
                                                               IntrinsicGeometryInterface& geom,
                                                               const EdgeData<bool>& inPath,
                                                               const VertexData<bool>& extraMarkedVertices);


  // add a path to an existing network
  // input should be a vector of INTRINSIC halfedges
  // TODO fornow, does not contain any logic to "layer" path on top of other existing paths; must be unobstructed
  void addPath(const std::vector<Halfedge>& path);

  // === Properties

  // The triangulation that the path is defined on
  std::unique_ptr<SignpostIntrinsicTriangulation> tri;
  ManifoldSurfaceMesh& mesh; // the intrinsic mesh

  // The collection of paths which make up the network
  std::vector<std::unique_ptr<FlipEdgePath>> paths;

  // The paths at each edge (and helpers to access them)
  EdgeData<std::deque<FlipPathSegment>> pathsAtEdge; // front of deque is paths at e.halfedge, back is e.halfedge.twin
  FlipPathSegment getOutsideSegment(Halfedge he);
  void popOutsideSegment(Halfedge he);
  void pushOutsideSegment(Halfedge he, FlipPathSegment p);
  FlipPathSegment getFirst(); // returns any if multiple
  FlipPathSegment getLast();  // returns any if multiple
  // HalfedgeData<int> pathCountAtHalfedge;

  // Vertices at which paths in the network terminate
  VertexData<bool> isMarkedVertex;

  // Queue of angles to be straightened, sorted by smallest angle. The list ref always refers to the wedge at
  // he.vertex().
  using WeightedAngle = std::tuple<double, SegmentAngleType, FlipPathSegment>;
  std::priority_queue<WeightedAngle, std::vector<WeightedAngle>, std::greater<WeightedAngle>> wedgeAngleQueue;
  void addToWedgeAngleQueue(const FlipPathSegment& pathSegmentNext);
  void addAllWedgesToAngleQueue();

  // Manage a unique set of IDs assigned to paths segments. The ordering of these IDs means nothing, and IDs are
  // never re-used. Paths call getNextUniquePathSegmentInd() to get IDs as needed.
  // (technically, each path could have its own pool of IDs, but using a global pool makes debugging easier)
  SegmentID nextUniquePathSegmentInd = 0; // monotonically increasing
  SegmentID getNextUniquePathSegmentInd();


  // === Options
  bool straightenAroundMarkedVertices = true;
  double EPS_ANGLE = 1e-5;

  // === Queries and accessors
  bool edgeInPath(Edge e);
  bool halfedgeInPath(Halfedge he);

  // Measure the length
  double length();

  // Counters
  long long int nFlips = 0;
  long long int nShortenIters = 0;

  // = Test straightness

  // Checks if there any other part of the path inside the wedge region
  bool wedgeIsClear(const FlipPathSegment& pathSegmentNext, SegmentAngleType type);
  bool wedgeIsClearEndpointsOnly(const FlipPathSegment& pathSegmentNext, SegmentAngleType type);

  // workhorse routine: determine the angle of the min wedge as well as the direction it points in
  std::tuple<double, double> measureSideAngles(Halfedge hePrev, Halfedge heNext);
  std::tuple<SegmentAngleType, double> locallyShortestTestWithType(Halfedge hePrev, Halfedge heNext);

  // A bunch of wrappers for locallyShortestTestWithType
  SegmentAngleType locallyShortestTest(Halfedge hePrev, Halfedge heNext);
  double minWedgeAngle(Halfedge hePrev, Halfedge heNext);
  double minWedgeAngle(const FlipPathSegment& segment);
  bool isStraight(double angleThresh = 1e-4);
  double minAngle();        // minimum over all angles
  double minAngleIsotopy(); // minimum over all angles, excluding those blocked by an path endpoint
  struct ShortestReturnBoth {
    SegmentAngleType minType;
    double minAngle;
    SegmentAngleType maxType;
    double maxAngle;
  };
  ShortestReturnBoth locallyShortestTestWithBoth(Halfedge hePrev, Halfedge heNext); // classifies both side angles

  // Get a path as a sequence of surface points along the mesh
  std::vector<std::vector<SurfacePoint>> getPathPolyline();
  std::vector<std::vector<SurfacePoint>> getAllEdgePolyline();
  std::vector<std::vector<SurfacePoint>> getPathPolyline(bool& wasPerfectOut);

  // Get a path as positions in 3D
  std::vector<std::vector<Vector3>> pathTo3D(const std::vector<std::vector<SurfacePoint>>& pathPoints); // helper
  std::vector<std::vector<Vector3>> getPathPolyline3D();
  std::vector<std::vector<Vector3>> getAllEdgePolyline3D();


  // Perform a whole bunch of sanity checks
  void validate();
  void validateHalfedgesOnly();
  bool intrinsicTriIsOriginal();

  // === Mutators

  // Shorten about a single wedge of a path (should probably be the one which makes the smallest angle)
  void locallyShortenAt(FlipPathSegment& pathSegmentNext, SegmentAngleType angleType);

  // Iteratively straighten the smallest angle in any path
  void iterativeShorten(size_t maxIterations = INVALID_IND, double maxRelativeLengthDecrease = 0.);

  // Flip non-path edges to make the triangulation Delaunay away from the path
  void makeDelaunay();

  // Perform intrinsic delaunay refinement while preserving curve
  void delaunayRefine(double areaThresh = std::numeric_limits<double>::infinity(), size_t maxInsertions = INVALID_IND,
                      double angleBound = 25.);

  // Split bent edges in the underlying triangulation. Useful for viz.
  void splitBentEdges(double angleDeg, size_t maxInsertions = INVALID_IND);

  // Perform one round of De Casteljau Bezier subdivision
  // Network must be a connected sequence of paths forming a curve
  void bezierSubdivide(size_t nRounds);
  void bezierSubdivideRecursive(size_t nRoundsRemaining, const Vertex firstControlCall, const Vertex lastControlCall);

  // === Visualization
  VertexPositionGeometry* posGeom = nullptr; // for visualization only!
  void savePath(std::string filenamePLY);
  void savePathOBJLine(std::string filenamePrefix, bool withAll = false);

  // === Helpers
  void updatePathAfterEdgeSplit(Halfedge origHe, Halfedge newHeFront);
  void processSingleEdgeLoop(FlipPathSegment& pathSegment, SegmentAngleType angleType);
  void purgeStaleQueueEntries(); // stop too many stale entries from accumulating
};

} // namespace surface
} // namespace geometrycentral
