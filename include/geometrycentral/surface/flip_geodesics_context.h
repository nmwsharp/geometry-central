#pragma once

#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/surface_point.h"

#include <cstdint>
#include <deque>
#include <memory>
#include <unordered_map>
#include <vector>

// Tuned FlipOut (Sharp & Crane 2020) for point-to-point queries, many per mesh.
//
// The shipped implementation (flip_geodesics.h, FlipEdgeNetwork) builds a fresh
// SignpostIntrinsicTriangulation and a fresh network for every query. FlipOutContext keeps one
// compact intrinsic triangulation per mesh (connectivity, edge lengths, signpost angles: the
// only quantities FlipOut reads), finds the initial edge path with A* over the input mesh,
// runs the same shortening loop as FlipEdgeNetwork::iterativeShorten on dense structures, and
// undoes the query's flips in O(1) each, so the next query starts from the original
// triangulation. The decisions (queue order, flip tests, angle updates) are the shipped ones
// operation for operation, so every query returns the same geodesic as the shipped path.

namespace geometrycentral {
namespace surface {

class VertexPositionGeometry;

class FlipOutContext {
public:
  FlipOutContext(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom);
  ~FlipOutContext();

  // FlipOut geodesic between two vertices; returns its length (0 if start == end or not
  // connected). The triangulation is back in its original state on return.
  double geodesicLength(Vertex start, Vertex end);
  // Same, and the path as points on the input surface (original edges reported directly,
  // flipped edges traced along the input like SignpostIntrinsicTriangulation does).
  double geodesic(Vertex start, Vertex end, std::vector<SurfacePoint>& pathOut);

  // The last query's path as intrinsic halfedge indices (the indices of the shipped
  // intrinsic mesh are the same), for identity checks.
  const std::vector<uint32_t>& lastPathHalfedges() const { return lastPath; }

  // Options (defaults as in FlipEdgeNetwork / IntrinsicTriangulation)
  bool useAStar = true; // Euclidean A* heuristic when the geometry has vertex positions
  double EPS_ANGLE = 1e-5;
  double triangleTestEPS = 1e-6;

  // Counters and timing of the last query (ms), for profiling
  long long nFlips = 0, nShortenIters = 0;
  double lastInitMs = 0, lastShortenMs = 0, lastRewindMs = 0;
  size_t lastPops = 0, lastPushes = 0;

private:
  static constexpr uint32_t NONE = 0xFFFFFFFFu;
  static constexpr uint32_t MULTI = 0xFFFFFFFEu;
  enum AngleType : uint8_t { Shortest = 0, LeftTurn = 1, RightTurn = 2 };

  ManifoldSurfaceMesh& inputMesh;
  IntrinsicGeometryInterface& geom;
  VertexPositionGeometry* positions = nullptr;
  size_t nV, nHe, nE;

  // Intrinsic triangulation. Vertices are the input vertices; halfedge/edge indices start as
  // the input's. Twin of he is he^1, edge is he/2, e.halfedge() is 2e (implicit twins).
  std::vector<uint32_t> heNext, heVertex; // mutated by flips
  std::vector<uint8_t> heInterior, vBoundary; // static
  std::vector<double> edgeLen, angle, angleSum; // angle: signpost angle per halfedge
  std::vector<uint8_t> edgeOrig;               // edge still an input edge (tracing shortcut)
  inline uint32_t twin(uint32_t he) const { return he ^ 1u; }
  inline uint32_t edge(uint32_t he) const { return he >> 1; }
  double cornerAngle(uint32_t he) const;
  void updateAngleFromCWNeighbor(uint32_t he);
  bool flipEdgeIfPossible(uint32_t e);

  // Undo log of the current query's flips
  struct Undo {
    uint32_t e;
    double len, angA, angB;
    uint8_t orig;
  };
  std::vector<Undo> undo;
  void rewind();

  // The path: segments in creation order (ids are their indices), doubly linked
  struct Seg {
    uint32_t he, prev, next;
    uint8_t alive;
  };
  std::vector<Seg> segs;
  // Segments layered on each edge, the shipped deque semantics: front is the side of
  // e.halfedge(), back the twin's. One id per edge in place; more than one spill to a map.
  std::vector<uint32_t> edgeSeg;
  std::unordered_map<uint32_t, std::deque<uint32_t>> multi;
  uint32_t outsideSegment(uint32_t he) const;
  void pushOutside(uint32_t he, uint32_t id);
  void popOutside(uint32_t he);
  inline bool edgeInPath(uint32_t e) const { return edgeSeg[e] != NONE; }

  // Wedge angle queue: min-heap on (angle, type, id)
  struct QEntry {
    double angle;
    uint32_t id;
    uint8_t type;
  };
  std::vector<QEntry> queue;
  void queuePush(double a, uint8_t type, uint32_t id);
  void addToWedgeAngleQueue(uint32_t id);

  // The shortening loop, as FlipEdgeNetwork
  void measureSideAngles(uint32_t hePrev, uint32_t heNext, double& left, double& right) const;
  double minWedgeAngle(uint32_t id) const;
  bool wedgeIsClear(uint32_t id, uint8_t type) const;
  void locallyShortenAt(uint32_t id, uint8_t type);
  void replacePathSegment(uint32_t nextID, uint8_t type, const std::vector<uint32_t>& newHalfedges);
  void iterativeShorten();
  std::vector<uint32_t> newPathBuf;

  // Query driver
  bool runQuery(Vertex start, Vertex end); // path left in place, not rewound
  void finishQuery();
  double pathLength() const;
  std::vector<uint32_t> lastPath;

  // Initial path: A*/Dijkstra over the input mesh with the shipped shortestEdgePath ordering
  struct Adj {
    uint32_t he, to;
    double len;
  };
  struct HeapEntry {
    double key, dist;
    uint32_t he;
  };
  std::vector<size_t> adjStart;
  std::vector<Adj> adj;
  struct SearchState { // per vertex: untouched / queued (best dist, halfedge) / settled (incoming halfedge)
    double dist;
    uint32_t stamp, he;
  };
  std::vector<SearchState> search;
  uint32_t curStamp = 0;
  std::vector<HeapEntry> heap;
  bool initialPath(Vertex start, Vertex end, std::vector<uint32_t>& pathOut);
};

} // namespace surface
} // namespace geometrycentral
