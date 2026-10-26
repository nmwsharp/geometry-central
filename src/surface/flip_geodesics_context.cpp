// Tuned FlipOut on a compact intrinsic triangulation. Mirrors the shipped implementation in
// flip_geodesics.cpp / signpost_intrinsic_triangulation.cpp step for step (same queue order,
// same flip test, same angle updates, same length check) on dense structures, for a single
// open path, with an O(1) undo per flip. The shipped code is untouched.

#include "geometrycentral/surface/flip_geodesics_context.h"

#include "geometrycentral/surface/trace_geodesic.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/utilities/elementary_geometry.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <stdexcept>

namespace geometrycentral {
namespace surface {

namespace {
using Clock = std::chrono::high_resolution_clock;
inline double msSince(Clock::time_point t0) {
  return std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
}
} // namespace

// ======================================================
// Construction: the intrinsic triangulation's initial state
// ======================================================

FlipOutContext::FlipOutContext(ManifoldSurfaceMesh& mesh_, IntrinsicGeometryInterface& geom_)
    : inputMesh(mesh_), geom(geom_) {

  if (!inputMesh.isCompressed()) {
    throw std::runtime_error("FlipOutContext needs a compressed mesh. Call mesh.compress().");
  }
  nV = inputMesh.nVertices();
  nHe = inputMesh.nHalfedges();
  nE = inputMesh.nEdges();
  positions = dynamic_cast<VertexPositionGeometry*>(&geom);
  geom.requireEdgeLengths();

  // Connectivity, as the input's
  heNext.resize(nHe);
  heVertex.resize(nHe);
  heInterior.resize(nHe);
  for (size_t i = 0; i < nHe; i++) {
    Halfedge he = inputMesh.halfedge(i);
    heNext[i] = (uint32_t)he.next().getIndex();
    heVertex[i] = (uint32_t)he.vertex().getIndex();
    heInterior[i] = he.isInterior() ? 1 : 0;
  }
  vBoundary.resize(nV);
  for (size_t v = 0; v < nV; v++) vBoundary[v] = inputMesh.vertex(v).isBoundary() ? 1 : 0;

  // Edge lengths, angle sums and signpost angles, computed as SignpostIntrinsicTriangulation
  // does (law of cosines per corner; corners summed in halfedge order; angles accumulated
  // clockwise from v.halfedge())
  edgeLen.resize(nE);
  for (size_t e = 0; e < nE; e++) edgeLen[e] = geom.edgeLengths[inputMesh.edge(e)];
  edgeOrig.assign(nE, 1);
  angleSum.assign(nV, 0.);
  std::vector<double> corner(nHe, 0.);
  for (size_t he = 0; he < nHe; he++) {
    if (!heInterior[he]) continue;
    corner[he] = cornerAngle((uint32_t)he);
    angleSum[heVertex[he]] += corner[he];
  }
  angle.assign(nHe, 0.);
  for (size_t v = 0; v < nV; v++) {
    double runningAngle = 0.;
    uint32_t firstHe = (uint32_t)inputMesh.vertex(v).halfedge().getIndex();
    uint32_t currHe = firstHe;
    do {
      angle[currHe] = runningAngle;
      if (!heInterior[currHe]) break;
      runningAngle += corner[currHe];
      currHe = twin(heNext[heNext[currHe]]);
    } while (currHe != firstHe);
  }

  // Path structures
  edgeSeg.assign(nE, NONE);

  // Search structures: compressed adjacency of the input mesh (per vertex, its outgoing
  // halfedges), one linear pass over the halfedges
  adjStart.assign(nV + 1, 0);
  for (size_t he = 0; he < nHe; he++) adjStart[heVertex[he] + 1]++;
  for (size_t v = 0; v < nV; v++) adjStart[v + 1] += adjStart[v];
  adj.resize(adjStart[nV]);
  std::vector<size_t> fill(adjStart.begin(), adjStart.end() - 1);
  for (size_t he = 0; he < nHe; he++) {
    Adj& a = adj[fill[heVertex[he]]++];
    a.he = (uint32_t)he;
    a.to = heVertex[twin((uint32_t)he)];
    a.len = edgeLen[edge((uint32_t)he)];
  }
  search.assign(nV, SearchState{0., 0u, 0u});
}

FlipOutContext::~FlipOutContext() { geom.unrequireEdgeLengths(); }

// ======================================================
// Intrinsic triangulation operations (SignpostIntrinsicTriangulation)
// ======================================================

// EdgeLengthGeometry::cornerAngle for the corner at the tail of he
double FlipOutContext::cornerAngle(uint32_t he) const {
  uint32_t heOpp = heNext[he];
  uint32_t heB = heNext[heOpp];
  double lOpp = edgeLen[edge(heOpp)];
  double lA = edgeLen[edge(he)];
  double lB = edgeLen[edge(heB)];
  double q = (lA * lA + lB * lB - lOpp * lOpp) / (2. * lA * lB);
  q = clamp(q, -1.0, 1.0);
  return std::acos(q);
}

// SignpostIntrinsicTriangulation::updateAngleFromCWNeighor (without the vector caches FlipOut
// never reads)
void FlipOutContext::updateAngleFromCWNeighbor(uint32_t he) {
  if (!heInterior[he]) {
    angle[he] = angleSum[heVertex[he]];
    return;
  }
  if (!heInterior[twin(he)]) {
    angle[he] = 0.;
    return;
  }
  uint32_t cwHe = heNext[twin(he)];
  double neighAngle = angle[cwHe];
  double cAngle = cornerAngle(cwHe);
  uint32_t v = heVertex[he];
  double a = neighAngle + cAngle;
  angle[he] = vBoundary[v] ? a : std::fmod(a, angleSum[v]); // standardizeAngle
}

// SignpostIntrinsicTriangulation::flipEdgeIfPossible + SurfaceMesh::flip
bool FlipOutContext::flipEdgeIfPossible(uint32_t e) {
  uint32_t ha1 = 2 * e, hb1 = ha1 + 1;
  if (!heInterior[ha1] || !heInterior[hb1]) return false; // boundary edge: the mesh flip refuses

  // layoutDiamond(he)
  uint32_t iHeA1 = heNext[ha1], iHeA2 = heNext[iHeA1];
  uint32_t iHeB1 = heNext[hb1], iHeB2 = heNext[iHeB1];
  double l01 = edgeLen[edge(iHeA1)];
  double l12 = edgeLen[edge(iHeA2)];
  double l23 = edgeLen[edge(iHeB1)];
  double l30 = edgeLen[edge(iHeB2)];
  double l02 = edgeLen[e];
  Vector2 p3{0., 0.};
  Vector2 p0{l30, 0.};
  Vector2 p2 = layoutTriangleVertex(p3, p0, l02, l23);
  Vector2 p1 = layoutTriangleVertex(p2, p0, l01, l12);

  // Geometrically flippable: both new triangles have positive signed area
  double A1 = cross(p1 - p0, p3 - p0);
  double A2 = cross(p3 - p2, p1 - p2);
  double areaEPS = triangleTestEPS * (A1 + A2);
  if (A1 < areaEPS || A2 < areaEPS) return false;

  double newLength = norm(p1 - p3);
  if (!std::isfinite(newLength)) return false;

  // Combinatorial flip
  uint32_t ha2 = iHeA1, ha3 = iHeA2;
  if (heNext[ha3] != ha1) return false;
  uint32_t hb2 = iHeB1, hb3 = iHeB2;
  if (heNext[hb3] != hb1) return false;
  if (ha2 == hb1 || hb2 == ha1) return false; // degree-1 vertex
  uint32_t vc = heVertex[ha3], vd = heVertex[hb3];
  heNext[ha1] = hb3;
  heNext[hb3] = ha2;
  heNext[ha2] = ha1;
  heNext[hb1] = ha3;
  heNext[ha3] = hb2;
  heNext[hb2] = hb1;
  heVertex[ha1] = vc;
  heVertex[hb1] = vd;

  edgeLen[e] = newLength;
  updateAngleFromCWNeighbor(ha1);
  updateAngleFromCWNeighbor(hb1);
  edgeOrig[e] = 0;
  return true;
}

// Undo every flip of the query, last first. After a flip the two faces read
// (ha1, hb3, ha2) and (hb1, ha3, hb2); the tails of ha2 and hb2 are the old endpoints.
void FlipOutContext::rewind() {
  while (!undo.empty()) {
    const Undo& u = undo.back();
    uint32_t ha1 = 2 * u.e, hb1 = ha1 + 1;
    uint32_t hb3 = heNext[ha1], ha2 = heNext[hb3];
    uint32_t ha3 = heNext[hb1], hb2 = heNext[ha3];
    heNext[ha1] = ha2;
    heNext[ha2] = ha3;
    heNext[ha3] = ha1;
    heNext[hb1] = hb2;
    heNext[hb2] = hb3;
    heNext[hb3] = hb1;
    heVertex[ha1] = heVertex[hb2];
    heVertex[hb1] = heVertex[ha2];
    edgeLen[u.e] = u.len;
    angle[ha1] = u.angA;
    angle[hb1] = u.angB;
    edgeOrig[u.e] = u.orig;
    undo.pop_back();
  }
}

// ======================================================
// Path segments on edges (FlipEdgeNetwork::pathsAtEdge deque semantics)
// ======================================================

uint32_t FlipOutContext::outsideSegment(uint32_t he) const {
  uint32_t e = edge(he), s = edgeSeg[e];
  if (s != MULTI) return s; // NONE or the single segment
  const std::deque<uint32_t>& d = multi.find(e)->second;
  return (he & 1u) == 0 ? d.front() : d.back();
}

void FlipOutContext::pushOutside(uint32_t he, uint32_t id) {
  uint32_t e = edge(he), s = edgeSeg[e];
  if (s == NONE) {
    edgeSeg[e] = id;
    return;
  }
  std::deque<uint32_t>& d = multi[e];
  if (s != MULTI) {
    d.clear();
    d.push_back(s);
    edgeSeg[e] = MULTI;
  }
  if ((he & 1u) == 0) {
    d.push_front(id);
  } else {
    d.push_back(id);
  }
}

void FlipOutContext::popOutside(uint32_t he) {
  uint32_t e = edge(he), s = edgeSeg[e];
  if (s == NONE) throw std::runtime_error("FlipOutContext: popped a segment from an empty edge");
  if (s != MULTI) {
    edgeSeg[e] = NONE;
    return;
  }
  std::deque<uint32_t>& d = multi[e];
  if ((he & 1u) == 0) {
    d.pop_front();
  } else {
    d.pop_back();
  }
  if (d.size() <= 1) {
    edgeSeg[e] = d.empty() ? NONE : d.front();
    multi.erase(e);
  }
}

// ======================================================
// Wedge angle queue
// ======================================================

namespace {
// the shipped queue is a std::priority_queue<tuple<angle, type, segment>, ..., greater<>>:
// the smallest (angle, type, id) comes out first
inline bool qLess(const double a, const uint8_t ta, const uint32_t ia, const double b, const uint8_t tb,
                  const uint32_t ib) {
  if (a != b) return a < b;
  if (ta != tb) return ta < tb;
  return ia < ib;
}
} // namespace

#define FOC_QCMP                                                                                                       \
  [](const QEntry& x, const QEntry& y) { return qLess(y.angle, y.type, y.id, x.angle, x.type, x.id); }

void FlipOutContext::queuePush(double a, uint8_t type, uint32_t id) {
  queue.push_back(QEntry{a, id, type});
  std::push_heap(queue.begin(), queue.end(), FOC_QCMP);
}

// FlipEdgeNetwork::measureSideAngles
void FlipOutContext::measureSideAngles(uint32_t hePrev, uint32_t heNext_, double& leftAngle, double& rightAngle) const {
  uint32_t v = heVertex[heNext_];
  double s = angleSum[v];
  double angleIn = angle[twin(hePrev)];
  double angleOut = angle[heNext_];
  bool isBoundary = vBoundary[v] != 0;

  if (angleIn < angleOut) {
    rightAngle = angleOut - angleIn;
  } else {
    rightAngle = isBoundary ? std::numeric_limits<double>::infinity() : (s - angleIn) + angleOut;
  }
  if (angleOut < angleIn) {
    leftAngle = angleIn - angleOut;
  } else {
    leftAngle = isBoundary ? std::numeric_limits<double>::infinity() : (s - angleOut) + angleIn;
  }
}

// FlipEdgeNetwork::addToWedgeAngleQueue, with locallyShortestTestWithBoth inlined
void FlipOutContext::addToWedgeAngleQueue(uint32_t id) {
  if (id == NONE) return;
  const Seg& seg = segs[id];
  if (seg.prev == NONE) return; // first halfedge of the open path
  uint32_t hePrev = segs[seg.prev].he, heNext_ = seg.he;

  double leftAngle, rightAngle;
  measureSideAngles(hePrev, heNext_, leftAngle, rightAngle);
  double minAngle, maxAngle;
  uint8_t minType, maxType;
  if (leftAngle < rightAngle) {
    minAngle = leftAngle;
    maxAngle = rightAngle;
    minType = leftAngle > (M_PI - EPS_ANGLE) ? Shortest : LeftTurn;
    maxType = rightAngle > (M_PI - EPS_ANGLE) ? Shortest : RightTurn;
  } else {
    minAngle = rightAngle;
    maxAngle = leftAngle;
    minType = rightAngle > (M_PI - EPS_ANGLE) ? Shortest : RightTurn;
    maxType = leftAngle > (M_PI - EPS_ANGLE) ? Shortest : LeftTurn;
  }
  if (minType == Shortest) return;
  queuePush(minAngle, minType, id);
  if (maxType == Shortest) return;
  queuePush(maxAngle, maxType, id);
}

// ======================================================
// Shortening (FlipEdgeNetwork)
// ======================================================

double FlipOutContext::minWedgeAngle(uint32_t id) const {
  const Seg& seg = segs[id];
  if (seg.prev == NONE) return M_PI;
  double leftAngle, rightAngle;
  measureSideAngles(segs[seg.prev].he, seg.he, leftAngle, rightAngle);
  return std::fmin(leftAngle, rightAngle);
}

bool FlipOutContext::wedgeIsClear(uint32_t id, uint8_t type) const {
  const Seg& seg = segs[id];
  if (seg.prev == NONE) throw std::runtime_error("called wedgeIsClear() beginning of openPath");
  uint32_t prevID = seg.prev;
  uint32_t hePrev = segs[prevID].he, heNext_ = seg.he;

  if (type == LeftTurn) {
    if (outsideSegment(hePrev) != prevID) return false;
    if (outsideSegment(heNext_) != id) return false;
    uint32_t heCurr = heNext[hePrev];
    while (heCurr != heNext_) {
      if (edgeInPath(edge(heCurr))) return false;
      heCurr = heNext[twin(heCurr)];
    }
  } else if (type == RightTurn) {
    if (outsideSegment(twin(hePrev)) != prevID) return false;
    if (outsideSegment(twin(heNext_)) != id) return false;
    uint32_t heCurr = twin(heNext[heNext[twin(hePrev)]]);
    while (heCurr != heNext_) {
      if (edgeInPath(edge(heCurr))) return false;
      heCurr = twin(heNext[heNext[heCurr]]);
    }
  } else {
    throw std::runtime_error("checked wedgeIsClear() with straight wedge, which doesn't make sense");
  }
  return true;
}

void FlipOutContext::locallyShortenAt(uint32_t nextID, uint8_t angleType) {
  const Seg seg = segs[nextID];
  uint32_t prevID = seg.prev;
  if (prevID == NONE) return; // first halfedge of the open path
  uint32_t heNext_ = seg.he, hePrev = segs[prevID].he;
  if (angleType == Shortest) return;
  nShortenIters++;
  if (prevID == nextID) throw std::runtime_error("FlipOutContext handles open paths only");

  double initPathLength = edgeLen[edge(hePrev)] + edgeLen[edge(heNext_)];

  // The straightening logic walks CW; a right turn is handled as the reversed left turn
  uint32_t sPrev, sNext;
  bool reversed;
  if (angleType == LeftTurn) {
    sPrev = hePrev;
    sNext = heNext_;
    reversed = false;
  } else {
    sPrev = twin(heNext_);
    sNext = twin(hePrev);
    reversed = true;
  }

  { // flip until a shorter path exists
    uint32_t sCurr = heNext[sPrev];
    uint32_t sPrevTwin = twin(sPrev);
    while (sCurr != sNext) {
      if (sCurr == sPrevTwin) { // never flip the first edge of the wedge
        sCurr = heNext[twin(sCurr)];
        continue;
      }
      uint32_t currEdge = edge(sCurr);
      Undo u{currEdge, edgeLen[currEdge], angle[2 * currEdge], angle[2 * currEdge + 1], edgeOrig[currEdge]};
      if (flipEdgeIfPossible(currEdge)) {
        nFlips++;
        undo.push_back(u);
        sCurr = twin(heNext[twin(sCurr)]); // re-check the previous edge
      } else {
        sCurr = heNext[twin(sCurr)];
      }
    }
  }

  // The new path along the boundary of the wedge, and its length
  double newPathLength = 0.;
  std::vector<uint32_t>& newPath = newPathBuf;
  newPath.clear();
  {
    uint32_t sCurr = heNext[sPrev];
    while (true) {
      uint32_t n = heNext[sCurr];
      newPath.push_back(twin(n));
      newPathLength += edgeLen[edge(n)];
      if (sCurr == sNext) break;
      sCurr = heNext[twin(sCurr)];
    }
  }

  // Only accept a path that is actually shorter (numerically unflippable edges can prevent it)
  if (newPathLength > initPathLength) return;

  if (reversed) {
    std::reverse(newPath.begin(), newPath.end());
    for (uint32_t& he : newPath) he = twin(he);
  }

  replacePathSegment(nextID, angleType, newPath);
}

// FlipEdgePath::replacePathSegment
void FlipOutContext::replacePathSegment(uint32_t nextID, uint8_t angleType, const std::vector<uint32_t>& newHalfedges) {
  uint32_t heNext_ = segs[nextID].he, prevID = segs[nextID].prev, nextNextID = segs[nextID].next;
  if (prevID == NONE) throw std::runtime_error("tried to to replace segment with next as first halfedge in list");
  uint32_t hePrev = segs[prevID].he, prevPrevID = segs[prevID].prev;

  // Part 1: remove the old edges from all structures
  if (angleType == LeftTurn) {
    popOutside(hePrev);
    popOutside(heNext_);
  } else {
    popOutside(twin(hePrev));
    popOutside(twin(heNext_));
  }
  segs[prevID].alive = 0;
  segs[nextID].alive = 0;

  bool replacedLength2Loop = false;
  if (prevPrevID == nextID) {
    replacedLength2Loop = true;
    prevPrevID = NONE;
    nextNextID = NONE;
  }

  // Part 2: add the new edges
  uint32_t loopPrevID = prevPrevID;
  uint32_t firstAddedID = NONE;
  for (uint32_t newHe : newHalfedges) {
    uint32_t newHeID = (uint32_t)segs.size();
    segs.push_back(Seg{newHe, loopPrevID, NONE, 1});
    if (angleType == LeftTurn) {
      pushOutside(twin(newHe), newHeID);
    } else {
      pushOutside(newHe, newHeID);
    }
    if (loopPrevID != NONE) segs[loopPrevID].next = newHeID;
    addToWedgeAngleQueue(newHeID);
    if (firstAddedID == NONE) firstAddedID = newHeID;
    loopPrevID = newHeID;
  }

  if (loopPrevID != NONE) segs[loopPrevID].next = nextNextID;
  if (nextNextID != NONE) {
    segs[nextNextID].prev = loopPrevID;
    addToWedgeAngleQueue(nextNextID); // its incoming angle probably changed
  }

  if (replacedLength2Loop) {
    uint32_t lastAddedID = loopPrevID;
    segs[firstAddedID].prev = lastAddedID;
    segs[lastAddedID].next = firstAddedID;
    addToWedgeAngleQueue(firstAddedID);
  }

  // Reconsider segments sharing the old edges: they may have become unblocked
  if (angleType == LeftTurn) {
    addToWedgeAngleQueue(outsideSegment(heNext_));
    addToWedgeAngleQueue(outsideSegment(hePrev));
  } else {
    addToWedgeAngleQueue(outsideSegment(twin(heNext_)));
    addToWedgeAngleQueue(outsideSegment(twin(hePrev)));
  }
}

// FlipEdgeNetwork::iterativeShorten (no iteration or length limit)
void FlipOutContext::iterativeShorten() {
  size_t nIterations = 0;
  while (!queue.empty()) {
    std::pop_heap(queue.begin(), queue.end(), FOC_QCMP);
    QEntry top = queue.back();
    queue.pop_back();

    if (!segs[top.id].alive) continue;                  // segment no longer exists
    if (minWedgeAngle(top.id) != top.angle) continue;   // angle has changed
    if (!wedgeIsClear(top.id, top.type)) continue;

    locallyShortenAt(top.id, top.type);
    nIterations++;

    if (nIterations % 1000 == 0) { // purge stale entries
      queue.clear();
      for (uint32_t id = 0; id < (uint32_t)segs.size(); id++) {
        if (segs[id].alive && segs[id].prev != NONE) addToWedgeAngleQueue(id);
      }
    }
  }
}

// ======================================================
// Queries
// ======================================================

bool FlipOutContext::runQuery(Vertex start, Vertex end) {
  auto t0 = Clock::now();
  std::vector<uint32_t>& path = lastPath;
  bool found = initialPath(start, end, path);
  lastInitMs = msSince(t0);
  if (!found) return false;

  t0 = Clock::now();
  segs.clear();
  queue.clear();
  undo.clear();
  uint32_t loopPrevID = NONE;
  for (uint32_t he : path) {
    uint32_t id = (uint32_t)segs.size();
    segs.push_back(Seg{he, loopPrevID, NONE, 1});
    pushOutside(he, id);
    if (loopPrevID != NONE) segs[loopPrevID].next = id;
    addToWedgeAngleQueue(id);
    loopPrevID = id;
  }
  iterativeShorten();

  // the final path, in order
  path.clear();
  uint32_t first = NONE;
  for (uint32_t id = 0; id < (uint32_t)segs.size() && first == NONE; id++) {
    if (segs[id].alive && segs[id].prev == NONE) first = id;
  }
  for (uint32_t id = first; id != NONE; id = segs[id].next) path.push_back(segs[id].he);
  lastShortenMs = msSince(t0);
  return true;
}

double FlipOutContext::pathLength() const {
  double L = 0.;
  for (uint32_t he : lastPath) L += edgeLen[edge(he)];
  return L;
}

void FlipOutContext::finishQuery() {
  auto t0 = Clock::now();
  for (const Seg& s : segs) {
    if (!s.alive) continue;
    uint32_t e = edge(s.he);
    if (edgeSeg[e] == MULTI) multi.erase(e);
    edgeSeg[e] = NONE;
  }
  rewind();
  segs.clear();
  queue.clear();
  lastRewindMs = msSince(t0);
}

double FlipOutContext::geodesicLength(Vertex start, Vertex end) {
  if (!runQuery(start, end)) return 0.;
  double L = pathLength();
  finishQuery();
  return L;
}

double FlipOutContext::geodesic(Vertex start, Vertex end, std::vector<SurfacePoint>& pathOut) {
  pathOut.clear();
  if (!runQuery(start, end)) return 0.;
  double L = pathLength();

  // trace every path halfedge along the input (SignpostIntrinsicTriangulation::
  // traceIntrinsicHalfedgeAlongInput), joined as FlipEdgeNetwork::getPathPolyline does
  for (uint32_t he : lastPath) {
    Vertex vA = inputMesh.vertex(heVertex[he]);
    Vertex vB = inputMesh.vertex(heVertex[twin(he)]);
    std::vector<SurfacePoint> thisTrace;
    if (edgeOrig[edge(he)]) {
      thisTrace = {SurfacePoint(vA), SurfacePoint(vB)};
    } else {
      uint32_t v = heVertex[he];
      double scaleFac = 1.0 / (angleSum[v] / (vBoundary[v] ? M_PI : 2. * M_PI));
      Vector2 traceVec = Vector2::fromAngle(angle[he] * scaleFac) * edgeLen[edge(he)];
      TraceOptions options;
      options.includePath = true;
      options.maxIters = inputMesh.nFaces() * 10;
      TraceGeodesicResult result = traceGeodesic(geom, SurfacePoint(vA), traceVec, options);
      if (trimTraceResult(result, vB)) {
        result.pathPoints.push_back(SurfacePoint(vB));
      } else {
        result = traceGeodesic(geom, SurfacePoint(vA), traceVec, options);
      }
      thisTrace = result.pathPoints;
    }
    if (!pathOut.empty() && onSameElement(pathOut.back(), thisTrace.front())) pathOut.pop_back();
    pathOut.insert(pathOut.end(), thisTrace.begin(), thisTrace.end());
  }

  finishQuery();
  return L;
}

// ======================================================
// Initial path
// ======================================================

// Dijkstra from start, stopping when end is settled, with the result of the shipped
// shortestEdgePath: that one orders its queue by (dist, halfedge index) and settles a vertex
// with the first entry popped for it. Here an entry for a vertex is only queued when it beats
// the vertex's best queued (dist, halfedge), which leaves the first-popped entry unchanged.
// With vertex positions the key is dist + Euclidean distance to the target (A*): consistent
// for Euclidean edge lengths, so the same entries win except in exact ties, which the
// (key, dist, halfedge) order resolves the Dijkstra way.
bool FlipOutContext::initialPath(Vertex start, Vertex end, std::vector<uint32_t>& path) {
  path.clear();
  if (start == end) return false;
  size_t s = start.getIndex(), t = end.getIndex();
  if (++curStamp >= 0x80000000u) { // the top bit marks queued vertices
    for (SearchState& st : search) st.stamp = 0;
    curStamp = 1;
  }
  const bool astar = useAStar && positions != nullptr;
  Vector3 pt = astar ? positions->inputVertexPositions[end] : Vector3::zero();
  auto h = [&](size_t v) -> double { return astar ? (positions->inputVertexPositions[v] - pt).norm() : 0.; };

  heap.clear();
  auto cmp = [](const HeapEntry& a, const HeapEntry& b) { // min-heap on (key, dist, he)
    if (a.key != b.key) return a.key > b.key;
    if (a.dist != b.dist) return a.dist > b.dist;
    return a.he > b.he;
  };
  const uint32_t QUEUED = 0x80000000u;
  auto relax = [&](size_t v, double dist) {
    for (size_t k = adjStart[v]; k < adjStart[v + 1]; k++) {
      const Adj& a = adj[k];
      SearchState& st = search[a.to];
      if (st.stamp == curStamp) continue; // settled
      double d = dist + a.len;
      if (st.stamp == (curStamp | QUEUED)) {
        if (d > st.dist || (d == st.dist && a.he > st.he)) continue; // cannot be popped first
      }
      st.stamp = curStamp | QUEUED;
      st.dist = d; // best queued (dist, halfedge) so far
      st.he = a.he;
      heap.push_back(HeapEntry{d + h(a.to), d, a.he});
      std::push_heap(heap.begin(), heap.end(), cmp);
      lastPushes++;
    }
  };
  lastPops = lastPushes = 0;
  search[s].stamp = curStamp;
  relax(s, 0.);
  while (!heap.empty()) {
    std::pop_heap(heap.begin(), heap.end(), cmp);
    HeapEntry e = heap.back();
    heap.pop_back();
    size_t v = heVertex[twin(e.he)];
    SearchState& st = search[v];
    if (st.stamp == curStamp) continue;
    st.stamp = curStamp; // settled: st.he is the incoming halfedge
    st.he = e.he;
    lastPops++;
    if (v == t) {
      size_t w = v;
      while (w != s) {
        uint32_t in = search[w].he;
        path.push_back(in);
        w = heVertex[in];
      }
      std::reverse(path.begin(), path.end());
      return true;
    }
    relax(v, e.dist);
  }
  return false; // not connected
}

} // namespace surface
} // namespace geometrycentral
