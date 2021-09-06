#include "geometrycentral/surface/flip_geodesics.h"

#include "geometrycentral/surface/mesh_graph_algorithms.h"

#include "happly.h"

namespace geometrycentral {
namespace surface {


FlipEdgePath::FlipEdgePath(FlipEdgeNetwork& network_, std::vector<Halfedge> halfedges_, bool isClosed_)
    : network(network_), isClosed(isClosed_) {

  if (halfedges_.empty()) {
    throw std::runtime_error("tried to create path from empty halfege list");
  }

  // Populate the list
  SegmentID loopPrevID = INVALID_IND;
  SegmentID firstID = INVALID_IND;
  for (Halfedge newHe : halfedges_) {


    // Get a new ID and entry
    SegmentID newHeID = network.getNextUniquePathSegmentInd();
    std::tuple<Halfedge, SegmentID, SegmentID> newEntry{newHe, loopPrevID, INVALID_IND};
    pathHeInfo[newHeID] = newEntry;

    // Try layering this path along the outside of the edge, on the side of the halfedge. Note that this MIGHT NOT be
    // the right thing to do, but we cannot distinguish more complicated cases without adding more info to the
    // constructor. Use other construction strategies to encode paths with interesting coincidence.
    network.pushOutsideSegment(newHe, {this, newHeID});


    // Track first
    if (firstID == INVALID_IND) {
      firstID = newHeID;
    }

    // Set the new ID as next on the previous entry
    if (loopPrevID != INVALID_IND) {
      std::get<2>(pathHeInfo[loopPrevID]) = newHeID;
    }

    // Add to the queue of wedge angles
    network.addToWedgeAngleQueue({this, newHeID});

    loopPrevID = newHeID;
  }
  SegmentID lastID = loopPrevID;

  // Gather info about start and end of path
  Halfedge firstHe = halfedges_.front();
  Vertex firstVert = firstHe.vertex();
  Halfedge lastHe = halfedges_.back();
  Vertex lastVert = lastHe.twin().vertex();


  if (isClosed) {

    // Sanity check that input actually connects
    if (lastHe.twin().vertex() != firstVert) {
      throw std::runtime_error("tried to construct closed path, but input halfedges do not form a loop");
    }

    // Connect the back to the front
    std::get<1>(pathHeInfo[firstID]) = lastID;
    std::get<2>(pathHeInfo[lastID]) = firstID;

  } else {

    // Set marked vertex flags for endpoints
    network.isMarkedVertex[firstVert] = true;
    network.isMarkedVertex[lastVert] = true;
  }
}


void FlipEdgePath::replacePathSegment(SegmentID nextID, SegmentAngleType angleType,
                                      const std::vector<Halfedge>& newHalfedges) {

  // We use the following halfedge names in order (note that some may be INVALID_IND)
  //      PrevPrev  ->   Prev    ->   Next   ->   NextNext
  //                     ^^^^^^^^^^^^^^^^^
  //                     straightening here

  // Get the prev and nextnext info
  Halfedge heNext, hePrev;
  SegmentID prevID, nextNextID, prevPrevID, UNUSED;
  std::tie(heNext, prevID, nextNextID) = pathHeInfo[nextID];
  size_t sizeBefore = pathHeInfo.size();

  if (prevID == INVALID_IND) {
    throw std::runtime_error("tried to to replace segment with next as first halfedge in list");
  }

  // Get the prevprev info
  std::tie(hePrev, prevPrevID, UNUSED) = pathHeInfo[prevID];

  // == Part 1: remove the old edges from all structures

  // Remove from the segment stacks on the edge (note that they MUST be on the outside or this whole method invalid)
  if (angleType == SegmentAngleType::LeftTurn) {
    network.popOutsideSegment(hePrev);
    network.popOutsideSegment(heNext);
  } else /* RightTurn */ {
    network.popOutsideSegment(hePrev.twin());
    network.popOutsideSegment(heNext.twin());
  }

  // Remove from this list
  pathHeInfo.erase(prevID);
  pathHeInfo.erase(nextID);

  // These handle deletion when the path is a length-2 loop to start
  bool replacedLength2Loop = false;
  if (prevPrevID == nextID) {
    replacedLength2Loop = true;
    prevPrevID = INVALID_IND;
    nextNextID = INVALID_IND;
  }

  // == Part 2: add the new edges to all structures
  SegmentID loopPrevID = prevPrevID;
  SegmentID firstAddedID = INVALID_IND;
  for (Halfedge newHe : newHalfedges) {
    // Get a new ID and entry
    SegmentID newHeID = network.getNextUniquePathSegmentInd();
    std::tuple<Halfedge, SegmentID, SegmentID> newEntry{newHe, loopPrevID, INVALID_IND};
    pathHeInfo[newHeID] = newEntry;

    // Add to edge stack
    if (angleType == SegmentAngleType::LeftTurn) {
      network.pushOutsideSegment(newHe.twin(), {this, newHeID});
    } else /* RightTurn */ {
      network.pushOutsideSegment(newHe, {this, newHeID});
    }

    // Set the new ID as next on the previous entry
    if (loopPrevID != INVALID_IND) {
      auto tupBefore = pathHeInfo[loopPrevID];
      std::get<2>(pathHeInfo[loopPrevID]) = newHeID;
      auto tupAfter = pathHeInfo[loopPrevID];
    }

    // add any new wedge angles to the list
    network.addToWedgeAngleQueue({this, newHeID});

    if (firstAddedID == INVALID_IND) {
      firstAddedID = newHeID;
    }
    loopPrevID = newHeID;
  }

  // Set the pointers to/from the nextnext segment
  if (loopPrevID != INVALID_IND) {
    std::get<2>(pathHeInfo[loopPrevID]) = nextNextID;
  }
  if (nextNextID != INVALID_IND) {
    std::get<1>(pathHeInfo[nextNextID]) = loopPrevID; // need to set prev ptr!

    // probably changed the angle incoming to the next halfedge
    network.addToWedgeAngleQueue({this, nextNextID});
  }

  // Connect the new front to new back
  if (replacedLength2Loop) {
    SegmentID lastAddedID = loopPrevID;
    std::get<1>(pathHeInfo[firstAddedID]) = lastAddedID;
    std::get<2>(pathHeInfo[lastAddedID]) = firstAddedID;

    network.addToWedgeAngleQueue({this, firstAddedID});
  }

  // Reconsider any outside segments that share the edge, as they may have become unblocked by this move
  if (angleType == SegmentAngleType::LeftTurn) {
    network.addToWedgeAngleQueue(network.getOutsideSegment(heNext)); // auto-exits if null
    network.addToWedgeAngleQueue(network.getOutsideSegment(hePrev));
  } else /* RightTurn */ {
    network.addToWedgeAngleQueue(network.getOutsideSegment(heNext.twin()));
    network.addToWedgeAngleQueue(network.getOutsideSegment(hePrev.twin()));
  }
}

size_t FlipEdgePath::size() { return pathHeInfo.size(); }

std::vector<Halfedge> FlipEdgePath::getHalfedgeList() {

  // Find the last halfedge with `next` == INVALID_IND, if one exists
  SegmentID lastID;
  for (auto it : pathHeInfo) {
    // Gather values
    SegmentID currID = it.first;
    SegmentID prevID, nextID;
    Halfedge currHe;
    std::tie(currHe, prevID, nextID) = it.second;

    lastID = currID;
    if (nextID == INVALID_IND) {
      // since we break after setting, will be left with lastID == currID
      break;
    }
  }

  // Walk backwards along path to reconstruct
  std::vector<Halfedge> result;
  SegmentID walkID = lastID;
  SegmentID firstID = walkID; // used to exit if a loop
  while (walkID != INVALID_IND) {
    SegmentID prevID, nextID;
    Halfedge currHe;
    std::tie(currHe, prevID, nextID) = pathHeInfo[walkID];
    result.push_back(currHe);
    walkID = prevID;
    if (walkID == firstID) break; // test at end so we don't insta-out
  }

  std::reverse(result.begin(), result.end()); // flip to have expected direction
  return result;
}

FlipEdgeNetwork::FlipEdgeNetwork(ManifoldSurfaceMesh& mesh_, IntrinsicGeometryInterface& inputGeom,
                                 const std::vector<std::vector<Halfedge>>& hePaths, VertexData<bool> extraMarkedVerts)
    : tri(std::unique_ptr<SignpostIntrinsicTriangulation>(new SignpostIntrinsicTriangulation(mesh_, inputGeom))),
      mesh(*(tri->intrinsicMesh)), pathsAtEdge(mesh), isMarkedVertex(mesh, false) {


  // Build initial paths from the vectors of edges (path constructor updates other structures of this class)
  for (const std::vector<Halfedge>& hePath : hePaths) {
    // Assumes that path is closed if it ends where it starts
    // (this might be a problem as the default one day, but for now it is overwhelmingly likely to be what we want
    Halfedge firstHe = hePath.front();
    Halfedge lastHe = hePath.back();
    bool isClosed = firstHe.vertex() == lastHe.twin().vertex();


    // Convert to the corresponding intrinsic halfedges
    std::vector<Halfedge> intPath(hePath.size());
    for (size_t i = 0; i < hePath.size(); i++) {
      intPath[i] = mesh.halfedge(hePath[i].getIndex());
    }

    paths.emplace_back(new FlipEdgePath(*this, intPath, isClosed));
  }

  // Mark any additional verts
  if (extraMarkedVerts.size() > 0) {
    for (Vertex v : mesh.vertices()) {
      if (extraMarkedVerts[v.getIndex()]) {
        isMarkedVertex[v] = true;
      }
    }
  }

  // Make sure everything is good to go
  validate();
}


void FlipEdgeNetwork::addPath(const std::vector<Halfedge>& hePath) {
  // Assumes that path is closed if it ends where it starts
  // (this might be a problem as the default one day, but for now it is overwhelmingly likely to be what we want
  Halfedge firstHe = hePath.front();
  Halfedge lastHe = hePath.back();
  bool isClosed = firstHe.vertex() == lastHe.twin().vertex();
  paths.emplace_back(new FlipEdgePath(*this, hePath, isClosed));
}

std::unique_ptr<FlipEdgeNetwork> FlipEdgeNetwork::constructFromDijkstraPath(ManifoldSurfaceMesh& mesh_,
                                                                            IntrinsicGeometryInterface& geom,
                                                                            Vertex startVert, Vertex endVert) {

  // Get the Dijkstra path
  std::vector<Halfedge> dijkstraPath = shortestEdgePath(geom, startVert, endVert);
  if (dijkstraPath.empty()) {
    // Not connected, or same vertex
    return std::unique_ptr<FlipEdgeNetwork>();
  }

  return std::unique_ptr<FlipEdgeNetwork>(new FlipEdgeNetwork(mesh_, geom, {dijkstraPath}));
}

std::unique_ptr<FlipEdgeNetwork> FlipEdgeNetwork::constructFromPiecewiseDijkstraPath(ManifoldSurfaceMesh& mesh_,
                                                                                     IntrinsicGeometryInterface& geom,
                                                                                     std::vector<Vertex> points,
                                                                                     bool closed, bool markInterior) {


  std::vector<Halfedge> halfedges;
  VertexData<bool> extraMark(geom.mesh, false);

  size_t end = closed ? points.size() : points.size() - 1;
  for (size_t i = 0; i < end; i++) {
    Vertex vA = points[i];
    Vertex vB = points[(i + 1) % points.size()];
    std::vector<Halfedge> dijkstraPath = shortestEdgePath(geom, vA, vB);

    if (markInterior) {
      extraMark[vA] = true;
      extraMark[vB] = true;
    }

    if (dijkstraPath.empty()) {
      // Not connected, or same vertex
      return std::unique_ptr<FlipEdgeNetwork>();
    }

    halfedges.insert(halfedges.end(), dijkstraPath.begin(), dijkstraPath.end());
  }

  return std::unique_ptr<FlipEdgeNetwork>(new FlipEdgeNetwork(mesh_, geom, {halfedges}, extraMark));
}


std::unique_ptr<FlipEdgeNetwork> FlipEdgeNetwork::constructFromEdgeSet(ManifoldSurfaceMesh& mesh_,
                                                                       IntrinsicGeometryInterface& geom,
                                                                       const EdgeData<bool>& inPath,
                                                                       const VertexData<bool>& extraMarkedVertices) {
  ManifoldSurfaceMesh& mesh = mesh_;

  std::vector<std::vector<Halfedge>> allHalfedges;

  // Endpoint vertices will be those with != 2 incident path edges
  VertexData<int> endpointCount(mesh, 0);
  for (Edge e : mesh.edges()) {
    if (inPath[e]) {
      endpointCount[e.halfedge().vertex()]++;
      endpointCount[e.halfedge().twin().vertex()]++;
    }
  }

  VertexData<bool> isEndpoint(mesh, false);
  for (Vertex v : mesh.vertices()) {
    if (extraMarkedVertices[v] || (endpointCount[v] != 0 && endpointCount[v] != 2)) {
      isEndpoint[v] = true;
    }
  }


  // Walk paths between the endpoints
  EdgeData<bool> walked(mesh, false);
  for (Halfedge heStart : mesh.halfedges()) {

    // Check if we should start a walk
    if (!inPath[heStart.edge()] || walked[heStart.edge()] || !isEndpoint[heStart.vertex()]) continue;

    // Start a new path
    allHalfedges.emplace_back();
    std::vector<Halfedge>& newPath = allHalfedges.back();

    // Walk along the path until we reach an endpoint
    Halfedge heCurr = heStart;
    while (true) {
      walked[heCurr.edge()] = true;
      newPath.push_back(heCurr);

      if (isEndpoint[heCurr.twin().vertex()]) break;

      // find next edge
      for (Halfedge heOther : heCurr.twin().vertex().outgoingHalfedges()) {
        if (heOther.twin() != heCurr && inPath[heOther.edge()]) {
          heCurr = heOther;
          break;
        }
      }
    }
  }

  // Any remaining paths are closed loops
  for (Halfedge heStart : mesh.halfedges()) {

    // Check if we should start a walk
    if (!inPath[heStart.edge()] || walked[heStart.edge()]) continue;

    // Start a new path
    allHalfedges.emplace_back();
    std::vector<Halfedge>& newPath = allHalfedges.back();

    // Walk along the path until we reach an endpoint
    Halfedge heCurr = heStart;
    do {
      walked[heCurr.edge()] = true;
      newPath.push_back(heCurr);

      // find next edge
      for (Halfedge heOther : heCurr.twin().vertex().outgoingHalfedges()) {
        if (heOther.twin() != heCurr && inPath[heOther.edge()]) {
          heCurr = heOther;
          break;
        }
      }
    } while (heCurr != heStart);
  }

  return std::unique_ptr<FlipEdgeNetwork>(new FlipEdgeNetwork(mesh_, geom, allHalfedges));
}


double FlipEdgeNetwork::minWedgeAngle(const FlipPathSegment& pathSegment) {
  FlipEdgePath& edgePath = *pathSegment.path;
  SegmentID nextID = pathSegment.id;
  SegmentID prevID, UNUSED;
  Halfedge heNext;
  std::tie(heNext, prevID, UNUSED) = edgePath.pathHeInfo[nextID];
  if (prevID == INVALID_IND) {
    // This is the first halfedge in a not-closed path
    return M_PI;
  }
  Halfedge hePrev = std::get<0>(edgePath.pathHeInfo[prevID]);

  return minWedgeAngle(hePrev, heNext);
}

double FlipEdgeNetwork::minWedgeAngle(Halfedge hePrev, Halfedge heNext) {
  return std::get<1>(locallyShortestTestWithType(hePrev, heNext));
}


SegmentAngleType FlipEdgeNetwork::locallyShortestTest(Halfedge hePrev, Halfedge heNext) {
  return std::get<0>(locallyShortestTestWithType(hePrev, heNext));
}

double FlipEdgeNetwork::minAngle() {

  double minAngle = std::numeric_limits<double>::infinity();

  for (auto& epPtr : paths) {
    FlipEdgePath& path = *epPtr;

    int prevInvalidCount = 0;
    int nextInvalidCount = 0;
    for (auto it : path.pathHeInfo) {

      // Gather values
      SegmentID currID = it.first;
      SegmentID prevID, nextID;
      Halfedge currHe;
      std::tie(currHe, prevID, nextID) = it.second;
      if (prevID == INVALID_IND) continue;

      Halfedge prevHe = std::get<0>(path.pathHeInfo[prevID]);

      double angle = minWedgeAngle(prevHe, currHe);
      minAngle = std::fmin(minAngle, angle);
    }
  }

  return minAngle;
}

double FlipEdgeNetwork::minAngleIsotopy() {

  double minAngle = std::numeric_limits<double>::infinity();

  for (auto& epPtr : paths) {
    FlipEdgePath& path = *epPtr;

    int prevInvalidCount = 0;
    int nextInvalidCount = 0;
    for (auto it : path.pathHeInfo) {

      // Gather values
      SegmentID currID = it.first;
      SegmentID prevID, nextID;
      Halfedge currHe;
      std::tie(currHe, prevID, nextID) = it.second;
      if (prevID == INVALID_IND) continue;

      Halfedge prevHe = std::get<0>(path.pathHeInfo[prevID]);

      ShortestReturnBoth result = locallyShortestTestWithBoth(prevHe, currHe);

      // Skip if the wedge if blocked
      FlipPathSegment seg{epPtr.get(), currID};
      if (result.minType != SegmentAngleType::Shortest && !wedgeIsClearEndpointsOnly(seg, result.minType)) {
        continue;
      }

      // Skip if the vertex is marked and not straightening
      if (!straightenAroundMarkedVertices && isMarkedVertex[currHe.twin().vertex()]) continue;

      minAngle = std::fmin(minAngle, result.minAngle);
    }
  }

  return minAngle;
}

std::tuple<double, double> FlipEdgeNetwork::measureSideAngles(Halfedge hePrev, Halfedge heNext) {
  Vertex v = heNext.vertex();
  double s = tri->vertexAngleSums[v];

  double angleIn = tri->signpostAngle[hePrev.twin()];
  double angleOut = tri->signpostAngle[heNext];
  bool isBoundary = v.isBoundary();

  // Compute right angle
  double rightAngle;
  if (angleIn < angleOut) {
    rightAngle = angleOut - angleIn;
  } else {
    if (isBoundary) {
      rightAngle = std::numeric_limits<double>::infinity();
    } else {
      rightAngle = (s - angleIn) + angleOut;
    }
  }

  // Compute left angle
  double leftAngle;
  if (angleOut < angleIn) {
    leftAngle = angleIn - angleOut;
  } else {
    if (isBoundary) {
      leftAngle = std::numeric_limits<double>::infinity();
    } else {
      leftAngle = (s - angleOut) + angleIn;
    }
  }

  return std::tuple<double, double>{leftAngle, rightAngle};
}

std::tuple<SegmentAngleType, double> FlipEdgeNetwork::locallyShortestTestWithType(Halfedge hePrev, Halfedge heNext) {

  if (hePrev == Halfedge()) return std::make_tuple(SegmentAngleType::Shortest, std::numeric_limits<double>::infinity());

  double leftAngle, rightAngle;
  std::tie(leftAngle, rightAngle) = measureSideAngles(hePrev, heNext);
  double minAngle = std::fmin(leftAngle, rightAngle);

  // Classify
  if (leftAngle < rightAngle) {
    if (leftAngle > (M_PI - EPS_ANGLE)) {
      return std::make_tuple(SegmentAngleType::Shortest, minAngle);
    }
    return std::make_tuple(SegmentAngleType::LeftTurn, minAngle);
  } else {
    if (rightAngle > (M_PI - EPS_ANGLE)) {
      return std::make_tuple(SegmentAngleType::Shortest, minAngle);
    }
    return std::make_tuple(SegmentAngleType::RightTurn, minAngle);
  }
}

FlipEdgeNetwork::ShortestReturnBoth FlipEdgeNetwork::locallyShortestTestWithBoth(Halfedge hePrev, Halfedge heNext) {
  ShortestReturnBoth result{SegmentAngleType::Shortest, std::numeric_limits<double>::infinity(),
                            SegmentAngleType::Shortest, std::numeric_limits<double>::infinity()};


  if (hePrev == Halfedge()) return result;

  double leftAngle, rightAngle;
  std::tie(leftAngle, rightAngle) = measureSideAngles(hePrev, heNext);
  double minAngle = std::fmin(leftAngle, rightAngle);

  // Classify
  if (leftAngle < rightAngle) {
    result.minAngle = leftAngle;
    result.maxAngle = rightAngle;
    if (leftAngle > (M_PI - EPS_ANGLE)) {
      result.minType = SegmentAngleType::Shortest;
    } else {
      result.minType = SegmentAngleType::LeftTurn;
    }
    if (rightAngle > (M_PI - EPS_ANGLE)) {
      result.maxType = SegmentAngleType::Shortest;
    } else {
      result.maxType = SegmentAngleType::RightTurn;
    }
  } else {
    result.minAngle = rightAngle;
    result.maxAngle = leftAngle;
    if (rightAngle > (M_PI - EPS_ANGLE)) {
      result.minType = SegmentAngleType::Shortest;
    } else {
      result.minType = SegmentAngleType::RightTurn;
    }
    if (leftAngle > (M_PI - EPS_ANGLE)) {
      result.maxType = SegmentAngleType::Shortest;
    } else {
      result.maxType = SegmentAngleType::LeftTurn;
    }
  }

  return result;
}

bool FlipEdgeNetwork::wedgeIsClear(const FlipPathSegment& pathSegmentNext, SegmentAngleType type) {

  // WARNING code duplications with endpoints only version

  // TODO handle checks in case where the path lollipops out and back along a single edge

  // Gather values
  FlipEdgePath& edgePath = *pathSegmentNext.path;
  SegmentID nextID = pathSegmentNext.id;
  SegmentID prevID, UNUSED;
  Halfedge heNext;
  std::tie(heNext, prevID, UNUSED) = edgePath.pathHeInfo[nextID];
  if (prevID == INVALID_IND) {
    // This is the first halfedge in a not-closed path
    throw std::runtime_error("called wedgeIsClear() beginning of openPath");
  }
  Halfedge hePrev = std::get<0>(edgePath.pathHeInfo[prevID]);
  FlipPathSegment pathSegmentPrev{&edgePath, prevID};


  // Gather values
  Vertex prevVert = hePrev.vertex();
  Vertex middleVert = heNext.vertex();
  Vertex nextVert = heNext.twin().vertex();

  // Used to disable straightening around certain vertices, but usually this isn't what you want: the algorithm should
  // still be able to straighten if a path touches and endpoint vertex but is not obstructed.
  if (!straightenAroundMarkedVertices && isMarkedVertex[middleVert]) {
    return false;
  }

  // Split to cases based on which side the wedge faces. Either way, we're iterating around the wedge making sure there
  // are no path edges in the way.
  switch (type) {
  case SegmentAngleType::Shortest: {
    throw std::runtime_error("checked wedgeIsClear() with straight wedge, which doesn't make sense");
    break;
  }
  case SegmentAngleType::LeftTurn: {

    // Check bounding edges
    if (getOutsideSegment(hePrev) != pathSegmentPrev) return false;
    if (getOutsideSegment(heNext) != pathSegmentNext) return false;

    // Orbit incident edges in wedge, each must have no path edges
    Halfedge heCurr = hePrev.next();
    while (heCurr != heNext) {
      if (edgeInPath(heCurr.edge())) return false;
      heCurr = heCurr.twin().next();
    }
    break;
  }
  case SegmentAngleType::RightTurn: {

    // Check bounding edges
    if (getOutsideSegment(hePrev.twin()) != pathSegmentPrev) return false;
    if (getOutsideSegment(heNext.twin()) != pathSegmentNext) return false;

    // Orbit incident edges in wedge, each must have no path edges
    Halfedge heCurr = hePrev.twin().next().next().twin();
    while (heCurr != heNext) {
      if (edgeInPath(heCurr.edge())) return false;
      heCurr = heCurr.next().next().twin();
    }
    break;
  }
  }

  return true;
}

// Like wedgeIsClear(), but only returns true if block by a path endpoint, rather than an interior portion of hte path.
// Useful for testing isotopy classes.
bool FlipEdgeNetwork::wedgeIsClearEndpointsOnly(const FlipPathSegment& pathSegmentNext, SegmentAngleType type) {

  // WARNING code duplications with endpoints only version

  // Gather values
  FlipEdgePath& edgePath = *pathSegmentNext.path;
  SegmentID nextID = pathSegmentNext.id;
  SegmentID prevID, UNUSED;
  Halfedge heNext;
  std::tie(heNext, prevID, UNUSED) = edgePath.pathHeInfo[nextID];
  if (prevID == INVALID_IND) {
    // This is the first halfedge in a not-closed path
    throw std::runtime_error("called wedgeIsClear() beginning of openPath");
  }
  Halfedge hePrev = std::get<0>(edgePath.pathHeInfo[prevID]);
  FlipPathSegment pathSegmentPrev{&edgePath, prevID};


  // Gather values
  Vertex prevVert = hePrev.vertex();
  Vertex middleVert = heNext.vertex();
  Vertex nextVert = heNext.twin().vertex();


  // Split to cases based on which side the wedge faces. Either way, we're iterating around the wedge making sure there
  // are no path edges in the way.
  switch (type) {
  case SegmentAngleType::Shortest: {
    throw std::runtime_error("checked wedgeIsClear() with straight wedge, which doesn't make sense");
    break;
  }
  case SegmentAngleType::LeftTurn: {

    // Check bounding edges
    if (getOutsideSegment(hePrev) != pathSegmentPrev && getOutsideSegment(hePrev).isEndpoint()) return false;
    if (getOutsideSegment(heNext) != pathSegmentNext && getOutsideSegment(heNext).isEndpoint()) return false;

    // Orbit incident edges in wedge, each must have no path edges
    Halfedge heCurr = hePrev.next();
    while (heCurr != heNext) {
      for (FlipPathSegment& p : pathsAtEdge[heCurr.edge()]) {
        if (p.isEndpoint()) return false;
      }
      heCurr = heCurr.twin().next();
    }
    break;
  }
  case SegmentAngleType::RightTurn: {

    // Check bounding edges
    if (getOutsideSegment(hePrev.twin()) != pathSegmentPrev && getOutsideSegment(hePrev.twin()).isEndpoint())
      return false;
    if (getOutsideSegment(heNext.twin()) != pathSegmentNext && getOutsideSegment(heNext.twin()).isEndpoint())
      return false;

    // Orbit incident edges in wedge, each must have no path edges
    Halfedge heCurr = hePrev.twin().next().next().twin();
    while (heCurr != heNext) {
      for (FlipPathSegment& p : pathsAtEdge[heCurr.edge()]) {
        if (p.isEndpoint()) return false;
      }
      heCurr = heCurr.next().next().twin();
    }
    break;
  }
  }

  return true;
}


void FlipEdgeNetwork::locallyShortenAt(FlipPathSegment& pathSegment, SegmentAngleType angleType) {

  // Gather values
  FlipEdgePath& edgePath = *pathSegment.path;
  SegmentID nextID = pathSegment.id;
  SegmentID prevID, UNUSED;
  Halfedge heNext;
  std::tie(heNext, prevID, UNUSED) = edgePath.pathHeInfo[nextID];
  if (prevID == INVALID_IND) {
    // This is the first halfedge in a not-closed path
    return;
  }
  Halfedge hePrev = std::get<0>(edgePath.pathHeInfo[prevID]);


  // Gather values
  Vertex prevVert = hePrev.vertex();
  Vertex middleVert = heNext.vertex();
  Vertex nextVert = heNext.twin().vertex();

  if (angleType == SegmentAngleType::Shortest) {
    // nothing to do here
    return;
  }
  nShortenIters++;

  // Special case for loop consisting of a single self-edge
  if (prevID == nextID) {
    processSingleEdgeLoop(pathSegment, angleType);
    return;
  }


  // Compute the initial path length
  double initPathLength = tri->edgeLengths[hePrev.edge()] + tri->edgeLengths[heNext.edge()];

  // The straightening logic below always walks CW, so flip the ordering if this is a right turn
  Halfedge sPrev, sNext;
  bool reversed;
  if (angleType == SegmentAngleType::LeftTurn) {
    sPrev = hePrev;
    sNext = heNext;
    reversed = false;
  } else {
    sPrev = heNext.twin();
    sNext = hePrev.twin();
    reversed = true;
  }

  { // == Main logic: flip until a shorter path exists
    Halfedge sCurr = sPrev.next();
    Halfedge sPrevTwin = sPrev.twin();
    while (sCurr != sNext) {

      // Don't want to flip the first edge of the wedge
      if (sCurr == sPrevTwin) {
        sCurr = sCurr.twin().next(); // advance to next edge
        continue;
      }

      // Gather values for the edge to be flipped
      Edge currEdge = sCurr.edge();
      double oldLen = tri->edgeLengths[currEdge]; // old values are used for rewinding
      double oldAngleA = tri->signpostAngle[currEdge.halfedge()];
      double oldAngleB = tri->signpostAngle[currEdge.halfedge().twin()];
      bool oldIsOrig = tri->edgeIsOriginal[currEdge];

      // Try to flip the edge. Note that flipping will only be possible iff \beta < \pi as in the formal algorithm
      // statement
      bool flipped = tri->flipEdgeIfPossible(currEdge);

      if (flipped) {
        nFlips++;

        // track data to support rewinding
        if (supportRewinding) {
          rewindRecord.emplace_back(currEdge, oldLen, oldAngleA, oldAngleB, oldIsOrig);
        }

        // Flip happened! Update data and continue processing
        // Re-check previous edge
        sCurr = sCurr.twin().next().twin();
      } else {

        sCurr = sCurr.twin().next(); // advance to next edge
      }
    }
  }

  // Build the list of edges representing the new path
  // measure the length of the new path along the boundary
  double newPathLength = 0.;
  std::vector<Halfedge> newPath;
  {
    Halfedge sCurr = sPrev.next();
    while (true) {
      newPath.push_back(sCurr.next().twin());
      newPathLength += tri->edgeLengths[sCurr.next().edge()];
      if (sCurr == sNext) break;
      sCurr = sCurr.twin().next();
    }
  }

  // Make sure the new path is actually shorter (this would never happen in the Reals, but can rarely happen if an edge
  // is numerically unflippable for floating point reasons)
  if (newPathLength > initPathLength) return;


  // Make sure the new path orientation matches the orientation of the input edges
  if (reversed) {
    std::reverse(newPath.begin(), newPath.end());
    for (Halfedge& he : newPath) {
      he = he.twin();
    }
  }

  // Replace the path segment with the new path
  // (most of the bookkeeping to update data structures happens in here)
  edgePath.replacePathSegment(nextID, angleType, newPath);
}


void FlipEdgeNetwork::processSingleEdgeLoop(FlipPathSegment& pathSegment, SegmentAngleType angleType) {

  // That annoying special case we need to handle for loops consisting of a single edge

  FlipEdgePath& edgePath = *pathSegment.path;
  SegmentID id = pathSegment.id;
  SegmentID UNUSED1, UNUSED2;
  Halfedge he;
  std::tie(he, UNUSED1, UNUSED2) = edgePath.pathHeInfo[id];


  // == Replace the old segment with the two opposite edges of the triangle

  switch (angleType) {
  case SegmentAngleType::Shortest: {
    // nothing to do, probably shouldn't even be here
    return;
    break;
  }
  case SegmentAngleType::LeftTurn: {

    Halfedge heFirst = he.next().next().twin();
    Halfedge heSecond = he.next().twin();

    SegmentID firstID = getNextUniquePathSegmentInd();
    SegmentID secondID = getNextUniquePathSegmentInd();

    edgePath.pathHeInfo.erase(id);
    popOutsideSegment(he);
    edgePath.pathHeInfo[firstID] = std::make_tuple(heFirst, secondID, secondID);
    edgePath.pathHeInfo[secondID] = std::make_tuple(heSecond, firstID, firstID);

    pushOutsideSegment(heFirst.twin(), FlipPathSegment{&edgePath, firstID});
    pushOutsideSegment(heSecond.twin(), FlipPathSegment{&edgePath, secondID});

    addToWedgeAngleQueue(FlipPathSegment{&edgePath, firstID});
    addToWedgeAngleQueue(FlipPathSegment{&edgePath, secondID});

    break;
  }
  case SegmentAngleType::RightTurn: {

    Halfedge heFirst = he.twin().next();
    Halfedge heSecond = he.twin().next().next();

    SegmentID firstID = getNextUniquePathSegmentInd();
    SegmentID secondID = getNextUniquePathSegmentInd();

    edgePath.pathHeInfo.erase(id);
    popOutsideSegment(he.twin());
    edgePath.pathHeInfo[firstID] = std::make_tuple(heFirst, secondID, secondID);
    edgePath.pathHeInfo[secondID] = std::make_tuple(heSecond, firstID, firstID);

    pushOutsideSegment(heFirst, FlipPathSegment{&edgePath, firstID});
    pushOutsideSegment(heSecond, FlipPathSegment{&edgePath, secondID});

    addToWedgeAngleQueue(FlipPathSegment{&edgePath, firstID});
    addToWedgeAngleQueue(FlipPathSegment{&edgePath, secondID});

    break;
  }
  }
}

void FlipEdgeNetwork::iterativeShorten(size_t maxIterations, double maxRelativeLengthDecrease) {

  bool checkLength = maxRelativeLengthDecrease != 0;
  double initLength = -777;
  if (checkLength) {
    initLength = length();
  }

  size_t nIterations = 0;

  while (!wedgeAngleQueue.empty() && (maxIterations == INVALID_IND || nIterations < maxIterations)) {

    // Get the smallest angle
    double minAngle = std::get<0>(wedgeAngleQueue.top());
    SegmentAngleType angleType = std::get<1>(wedgeAngleQueue.top());
    FlipPathSegment pathSegment = std::get<2>(wedgeAngleQueue.top());
    wedgeAngleQueue.pop();

    // Check if its a stale entry
    FlipEdgePath& path = *pathSegment.path;
    if (path.pathHeInfo.find(pathSegment.id) == path.pathHeInfo.end()) {
      continue; // segment no longer exists
    }
    double currAngle = minWedgeAngle(pathSegment);
    if (currAngle != minAngle) {
      continue; // angle has changed
    }

    // Make sure the wedge is clear
    // TODO I think we _might_ be able to argue that this check isn't necessary, and the wedge will always be clear as
    // long as we check it before inserting in to the queue
    if (!wedgeIsClear(pathSegment, angleType)) {
      continue;
    }

    // Shorten about that wedge
    locallyShortenAt(pathSegment, angleType);

    // validate();

    nIterations++;

    // Periodically purge the queue, since we can't remove from it and thus accumulate stale entries
    size_t purgeInterval = 1000; // note: basically does nothing; our queues pretty much never get this big
    if (nIterations % purgeInterval == 0) {
      purgeStaleQueueEntries();
    }

    // Quit if we pass the length threshold
    if (checkLength) {
      double currLength = length();
      if (currLength < maxRelativeLengthDecrease * initLength) {
        break;
      }
    }
  }
}

SegmentID FlipEdgeNetwork::getNextUniquePathSegmentInd() {
  SegmentID newInd = nextUniquePathSegmentInd;
  nextUniquePathSegmentInd++;
  return newInd;
}

void FlipEdgeNetwork::addToWedgeAngleQueue(const FlipPathSegment& pathSegment) {
  if (pathSegment.path == nullptr) return; // null input

  // Gather values
  FlipEdgePath& edgePath = *pathSegment.path;
  SegmentID nextID = pathSegment.id;
  SegmentID prevID, UNUSED;
  Halfedge heNext;
  std::tie(heNext, prevID, UNUSED) = edgePath.pathHeInfo[nextID];
  if (prevID == INVALID_IND) {
    // This is the first halfedge in a not-closed path
    return;
  }
  Halfedge hePrev = std::get<0>(edgePath.pathHeInfo[prevID]);

  // Measure the angles on both sides
  ShortestReturnBoth testResult = locallyShortestTestWithBoth(hePrev, heNext);

  // Check the smaller of the two angles
  if (testResult.minType == SegmentAngleType::Shortest) return; // If smaller is straight, nothing to do
  wedgeAngleQueue.emplace(testResult.minAngle, testResult.minType, pathSegment);

  // Check the larger of the two angles
  if (testResult.maxType == SegmentAngleType::Shortest) return;
  wedgeAngleQueue.emplace(testResult.maxAngle, testResult.maxType, pathSegment);

  // Note that there's an important implicit thing going on here: after straightening we only add possibly-blocked
  // wedges which are geometrically coincident with the just-straightened wedge. This is because we assume any others
  // that just got blocked have a strictly larger angle, and were added to the queue and are still hanging out there.
  // Thus it's important that we _don't_ try to optimize and test the wedge clearness here, since that would break the
  // logic.
}

void FlipEdgeNetwork::addAllWedgesToAngleQueue() {
  for (auto& epPtr : paths) {
    FlipEdgePath& path = *epPtr;

    for (auto it : path.pathHeInfo) {

      // Gather values
      SegmentID currID = it.first;
      SegmentID prevID, nextID;
      Halfedge currHe;
      std::tie(currHe, prevID, nextID) = it.second;

      if (prevID != INVALID_IND) {
        addToWedgeAngleQueue(FlipPathSegment{epPtr.get(), currID});
      }
    }
  }
}

void FlipEdgeNetwork::makeDelaunay() {

  // == Mark path edges as fixed
  EdgeData<bool> fixedEdges(tri->mesh);
  for (Edge e : tri->mesh.edges()) {
    fixedEdges[e] = edgeInPath(e);
  }
  tri->setMarkedEdges(fixedEdges);

  tri->flipToDelaunay();
}

void FlipEdgeNetwork::delaunayRefine(double areaThresh, size_t maxInsertions, double angleBound) {

  // == Mark path edges as fixed
  EdgeData<bool> fixedEdges(tri->mesh);
  for (Edge e : tri->mesh.edges()) {
    fixedEdges[e] = edgeInPath(e);
  }
  tri->setMarkedEdges(fixedEdges);

  // == Register a callback to maintain the path when edges are split
  auto updatePathOnSplit = [&](Edge oldE, Halfedge newHe1, Halfedge newHe2) {
    updatePathAfterEdgeSplit(oldE.halfedge(), newHe1);
  };
  auto callbackRef = tri->edgeSplitCallbackList.insert(std::end(tri->edgeSplitCallbackList), updatePathOnSplit);

  // == Refine!
  tri->delaunayRefine(angleBound, areaThresh, maxInsertions);

  tri->edgeSplitCallbackList.erase(callbackRef); // remove the callback we registered
}

void FlipEdgeNetwork::bezierSubdivide(size_t nRounds) {

  // disable straightening around marked vertices, we will use them
  bool oldStraightenSetting = straightenAroundMarkedVertices;
  straightenAroundMarkedVertices = false;


  // Ensure the curve is straight to start with
  iterativeShorten();

  bezierSubdivideRecursive(nRounds, getFirst().halfedge().vertex(), getLast().halfedge().twin().vertex());

  // restore
  straightenAroundMarkedVertices = oldStraightenSetting;
}

// Each call of this function introduces a new control point at the midpoint of firstControlRegion--lastControlRegion
// which lies exactly on the limit curve.
//
// We bound regions with vertices, rather than path segments, because the path segments are constantly changing as we
// straighten the path.
//
// lastControlRegion is inclusive
void FlipEdgeNetwork::bezierSubdivideRecursive(size_t nRoundsRemaining, const Vertex firstControlCall,
                                               const Vertex lastControlCall) {

  // NOTE: Much of this implementation assumes the curve is simple at all times
  if (nRoundsRemaining == 0) return;

  // Helpers
  auto findSegmentAfter = [&](Vertex v) {
    FlipPathSegment toReturn{nullptr, INVALID_IND};
    for (Halfedge he : v.outgoingHalfedges()) {
      if (!pathsAtEdge[he.edge()].empty()) {
        FlipPathSegment s = pathsAtEdge[he.edge()].front();
        if (s.halfedge() == he) {
          if (toReturn.id != INVALID_IND) throw std::runtime_error("multiple paths at vertex in bezier");
          toReturn = s;
        }
      }
    }
    return toReturn;
  };
  auto findSegmentBefore = [&](Vertex v) {
    FlipPathSegment toReturn{nullptr, INVALID_IND};
    for (Halfedge he : v.outgoingHalfedges()) {
      if (!pathsAtEdge[he.edge()].empty()) {
        FlipPathSegment s = pathsAtEdge[he.edge()].front();
        if (s.halfedge().twin() == he) {
          if (toReturn.id != INVALID_IND) throw std::runtime_error("multiple paths at vertex in bezier");
          toReturn = s;
        }
      }
    }
    return toReturn;
  };

  // Defines the active region
  Vertex firstControlActive = firstControlCall;
  Vertex lastControlActive = lastControlCall;
  Vertex newMidpoint;

  while (true) {
    FlipPathSegment firstSeg = findSegmentAfter(firstControlActive);
    FlipPathSegment lastSeg = findSegmentBefore(lastControlActive);

    // == Find each of the region between control points
    using SegRegion = std::tuple<FlipPathSegment, FlipPathSegment, double>;
    std::vector<SegRegion> regions;
    { // walk the path to find regions

      FlipPathSegment currSeg = firstSeg;
      FlipPathSegment regionStart = currSeg;
      double length = 0;
      while (true) {
        Halfedge currHe = currSeg.halfedge();

        length += tri->edgeLengths[currHe.edge()];

        // Finish the current region and start a new one
        if (isMarkedVertex[currHe.twin().vertex()]) {
          regions.emplace_back(regionStart, currSeg, length);

          if (currSeg == lastSeg) break;
          currSeg = currSeg.next();

          regionStart = currSeg;
          length = 0;
        } else {
          currSeg = currSeg.next();
        }
      }
    }
    bool isLast = regions.size() == 1;


    // == Unmark all old points on the interior of the active region
    // (except on the last iteration)
    if (!isLast) {
      for (auto& segRegion : regions) {
        FlipPathSegment currP = std::get<0>(segRegion);
        FlipPathSegment lastP = std::get<1>(segRegion);

        while (true) {
          if (currP != firstSeg) {
            Vertex vertBefore = currP.halfedge().vertex();
            isMarkedVertex[vertBefore] = false;
            FlipPathSegment s = findSegmentAfter(vertBefore);
            if (s.id != INVALID_IND) addToWedgeAngleQueue(s); // might need to straighten since we unmarked
          }

          if (currP == lastP) break;
          currP = currP.next();
        }
      }
    }


    // == Insert the midpoint of each region between control points, mark the new midpoints
    std::vector<Vertex> newControlPoints;
    for (auto& segRegion : regions) {

      // Find the segment where the midpoint occurs
      FlipPathSegment currP = std::get<0>(segRegion);
      FlipPathSegment lastP = std::get<1>(segRegion);
      double halfLen = std::get<2>(segRegion) * .5;
      double runningLen = 0;
      double useVertexEPS = 1e-4;
      while (true) {
        double nextLen = runningLen + tri->edgeLengths[currP.halfedge().edge()];
        if ((1. + useVertexEPS) * nextLen > halfLen) {
          break;
        }
        if (currP == lastP) throw std::runtime_error("couldn't find split segment");

        runningLen = nextLen;
        currP = currP.next();
      }


      // Split it
      FlipPathSegment splitP = currP;
      double tSplit = (halfLen - runningLen) / tri->edgeLengths[splitP.halfedge().edge()];

      // Case where the point we were going to insert is already present (or extremely close to) a vertex. Tends to
      // happen on regular grids. Use that point instead.
      if (tSplit > 1. - useVertexEPS) {
        Vertex existingPoint = splitP.halfedge().twin().vertex();

        // make sure it's marked
        isMarkedVertex[existingPoint] = true;
        newControlPoints.push_back(existingPoint);
      } else {


        splitP.splitEdge(tSplit);

        // Mark the new midpoint
        Vertex newControlP = splitP.halfedge().twin().vertex();
        isMarkedVertex[newControlP] = true;
        newControlPoints.push_back(newControlP);
      }
    }

    // Shrink  the active region
    firstControlActive = newControlPoints.front();
    lastControlActive = newControlPoints.back();

    // == Straighten to geodesic
    iterativeShorten();

    if (isLast) {
      newMidpoint = newControlPoints.front(); // must be size 1 list
      break;
    }
  }

  // Recurse on to both halves of the curve
  bezierSubdivideRecursive(nRoundsRemaining - 1, firstControlCall, newMidpoint);
  bezierSubdivideRecursive(nRoundsRemaining - 1, newMidpoint, lastControlCall);
}

void FlipEdgeNetwork::rewind() {
  if (!supportRewinding) {
    throw std::runtime_error(
        "Called FlipEdgeNetwork::rewind(), but rewinding is not supported. Set supportRewinding=true on construction.");
  }

  // TODO in theory, we might want to separate out the idea of undoing a sequence of flips, and of clearing the
  // represented paths. Right now this function always does both.

  // == Clear any stored paths
  // Clear out edge stacks
  for (auto& epPtr : paths) {
    // TODO maybe do this in a desctructor of FlipEdgePath instead? Right now this is the only place removals happen.
    FlipEdgePath& path = *epPtr;
    for (auto it : path.pathHeInfo) {
      SegmentID prevID, nextID;
      Halfedge currHe;
      std::tie(currHe, prevID, nextID) = it.second;
      pathsAtEdge[currHe.edge()].clear();
    }
  }
  // Clear out paths themselves, and any stored data
  paths.clear();
  wedgeAngleQueue = std::priority_queue<WeightedAngle, std::vector<WeightedAngle>,
                                        std::greater<WeightedAngle>>(); // no clear method, so create a new one

  // Undo the rewind stack
  while (!rewindRecord.empty()) {

    // Get the top element
    Edge edge;
    double oldLen, oldAngleA, oldAngleB;
    bool isOrig;
    std::tie(edge, oldLen, oldAngleA, oldAngleB, isOrig) = rewindRecord.back();
    rewindRecord.pop_back();

    // Undo the flip
    // bool flipped = tri->flipEdgeIfPossible(edge, 0.);
    tri->flipEdgeManual(edge, oldLen, oldAngleA, oldAngleB, isOrig, true);
  }
}

void FlipEdgeNetwork::reinitializePath(const std::vector<std::vector<Halfedge>>& newPaths) {

  // reset the data structure to its initial state
  rewind();

  // TODO shared code from constructor
  for (const std::vector<Halfedge>& hePath : newPaths) {
    // Assumes that path is closed if it ends where it starts
    // (this might be a problem as the default one day, but for now it is overwhelmingly likely to be what we want
    Halfedge firstHe = hePath.front();
    Halfedge lastHe = hePath.back();
    bool isClosed = firstHe.vertex() == lastHe.twin().vertex();

    // Convert to the corresponding intrinsic halfedges
    std::vector<Halfedge> intPath(hePath.size());
    for (size_t i = 0; i < hePath.size(); i++) {
      intPath[i] = mesh.halfedge(hePath[i].getIndex());
    }

    paths.emplace_back(new FlipEdgePath(*this, intPath, isClosed));
  }
}

bool FlipEdgeNetwork::edgeInPath(Edge e) { return !pathsAtEdge[e].empty(); }

bool FlipEdgeNetwork::halfedgeInPath(Halfedge he) { return !pathsAtEdge[he.edge()].empty(); }

FlipPathSegment FlipEdgeNetwork::getOutsideSegment(Halfedge he) {
  Edge e = he.edge();

  if (pathsAtEdge[e].empty()) {
    return {nullptr, INVALID_IND};
  }

  if (he == e.halfedge()) {
    return pathsAtEdge[e].front();
  } else {
    return pathsAtEdge[e].back();
  }
}

void FlipEdgeNetwork::popOutsideSegment(Halfedge he) {
  Edge e = he.edge();

  if (he == e.halfedge()) {
    pathsAtEdge[e].pop_front();
  } else {
    pathsAtEdge[e].pop_back();
  }
}

void FlipEdgeNetwork::pushOutsideSegment(Halfedge he, FlipPathSegment p) {
  Edge e = he.edge();

  if (he == e.halfedge()) {
    pathsAtEdge[e].push_front(p);
  } else {
    pathsAtEdge[e].push_back(p);
  }
}

FlipPathSegment FlipEdgeNetwork::getFirst() {

  for (auto& epPtr : paths) {
    FlipEdgePath& path = *epPtr;

    for (auto it : path.pathHeInfo) {

      // Gather values
      SegmentID currID = it.first;
      SegmentID prevID, nextID;
      Halfedge currHe;
      std::tie(currHe, prevID, nextID) = it.second;

      if (prevID == INVALID_IND) return FlipPathSegment{epPtr.get(), currID};
    }
  }

  throw std::runtime_error("could not find first segment");
}

FlipPathSegment FlipEdgeNetwork::getLast() {

  for (auto& epPtr : paths) {
    FlipEdgePath& path = *epPtr;

    for (auto it : path.pathHeInfo) {

      // Gather values
      SegmentID currID = it.first;
      SegmentID prevID, nextID;
      Halfedge currHe;
      std::tie(currHe, prevID, nextID) = it.second;

      if (nextID == INVALID_IND) return FlipPathSegment{epPtr.get(), currID};
    }
  }

  throw std::runtime_error("could not find last segment");
}

void FlipEdgeNetwork::validateHalfedgesOnly() {
  for (auto& epPtr : paths) {
    FlipEdgePath& path = *epPtr;

    int prevInvalidCount = 0;
    int nextInvalidCount = 0;
    for (auto it : path.pathHeInfo) {

      // Gather values
      SegmentID currID = it.first;
      SegmentID prevID, nextID;
      Halfedge currHe;
      std::tie(currHe, prevID, nextID) = it.second;

      // Check that the halfedge points to an actual halfedge
      if (currHe.getMesh() == nullptr) throw std::runtime_error("bad halfedge entry");
    }
  }
}

double FlipEdgeNetwork::length() {

  double length = 0.;

  for (auto& epPtr : paths) {
    FlipEdgePath& path = *epPtr;

    for (auto it : path.pathHeInfo) {

      // Gather values
      SegmentID currID = it.first;
      SegmentID prevID, nextID;
      Halfedge currHe;
      std::tie(currHe, prevID, nextID) = it.second;

      length += tri->edgeLengths[currHe.edge()];
    }
  }

  return length;
}

void FlipEdgeNetwork::validate() {

  // == Check that all paths are connected
  for (auto& epPtr : paths) {
    FlipEdgePath& path = *epPtr;

    int prevInvalidCount = 0;
    int nextInvalidCount = 0;
    for (auto it : path.pathHeInfo) {

      // Gather values
      SegmentID currID = it.first;
      SegmentID prevID, nextID;
      Halfedge currHe;
      std::tie(currHe, prevID, nextID) = it.second;

      // Check that the halfedge points to an actual halfedge
      if (currHe.getMesh() == nullptr) throw std::runtime_error("bad halfedge entry");

      // The only time prev and next should be invalid are on the endpoints of an open path
      if (prevID == INVALID_IND) {
        if (path.isClosed || prevInvalidCount > 0) throw std::runtime_error("too many invalid prev IDs");
        if (!isMarkedVertex[currHe.vertex()]) throw std::runtime_error("path begin endpoint isn't a marked vertex");
        prevInvalidCount++;
      }
      if (nextID == INVALID_IND) {
        if (path.isClosed || nextInvalidCount > 0) throw std::runtime_error("too many invalid next IDs");
        if (!isMarkedVertex[currHe.twin().vertex()])
          throw std::runtime_error("path end endpoint isn't a marked vertex");
        nextInvalidCount++;
      }

      // Check prev
      if (prevID != INVALID_IND) {
        Halfedge prevHe;
        SegmentID prevNextID, UNUSED;
        std::tie(prevHe, UNUSED, prevNextID) = path.pathHeInfo[prevID];
        if (prevNextID != currID) {
          throw std::runtime_error("prev next is not curr");
        }

        if (prevHe.twin().vertex() != currHe.vertex()) {
          std::cout << prevHe.twin().vertex().getIndex() << " " << currHe.vertex().getIndex() << std::endl;
          throw std::runtime_error("prev he not connected to curr");
        }
      }

      // Check next
      if (nextID != INVALID_IND) {
        Halfedge nextHe;
        SegmentID nextPrevID, UNUSED;
        std::tie(nextHe, nextPrevID, UNUSED) = path.pathHeInfo[nextID];
        if (nextPrevID != currID) {
          throw std::runtime_error("next prev is not curr");
        }

        if (currHe.twin().vertex() != nextHe.vertex()) {
          throw std::runtime_error("next he not connected to curr");
        }
      }
    }

    if (!path.isClosed) {
      if (prevInvalidCount != 1) throw std::runtime_error("open path does not have first halfedge");
      if (nextInvalidCount != 1) throw std::runtime_error("open path does not have last halfedge");
    }
  }

  // Check that all path endpoints are marked
  for (auto& epPtr : paths) {
    FlipEdgePath& path = *epPtr;

    if (path.isClosed) continue;

    for (auto it : path.pathHeInfo) {

      // Gather values
      SegmentID currID = it.first;
      SegmentID prevID, nextID;
      Halfedge currHe;
      std::tie(currHe, prevID, nextID) = it.second;

      if (prevID == INVALID_IND) {
        Vertex firstEndpoint = currHe.vertex();
        if (!isMarkedVertex[firstEndpoint]) {
          throw std::runtime_error("first endpoint of path is not a marked vertex");
        }
      }
      if (nextID == INVALID_IND) {
        Vertex lastEndpoint = currHe.twin().vertex();
        if (!isMarkedVertex[lastEndpoint]) {
          throw std::runtime_error("last endpoint of path is not a marked vertex");
        }
      }
    }
  }


  // Check that all path edges are noted as edges
  HalfedgeData<bool> halfedgeSeen(mesh, false);
  for (auto& epPtr : paths) {
    FlipEdgePath& path = *epPtr;

    for (auto it : path.pathHeInfo) {

      // Gather values
      SegmentID currID = it.first;
      SegmentID prevID, nextID;
      Halfedge currHe;
      std::tie(currHe, prevID, nextID) = it.second;

      if (!halfedgeInPath(currHe)) {
        throw std::runtime_error("path halfedge has halfedgeInPath() == false");
      }
      halfedgeSeen[currHe] = true;
    }
  }

  // Make sure no extra halfedges are marked which we did not see
  for (Edge e : mesh.edges()) {
    if (edgeInPath(e) && !(halfedgeSeen[e.halfedge()] || halfedgeSeen[e.halfedge().twin()])) {
      throw std::runtime_error("halfedge has entries but does not appear in path");
    }
  }
}


std::vector<std::vector<SurfacePoint>> FlipEdgeNetwork::getPathPolyline() {
  bool tmp;
  return getPathPolyline(tmp);
}
std::vector<std::vector<SurfacePoint>> FlipEdgeNetwork::getPathPolyline(bool& wasPerfectOut) {

  std::vector<std::vector<SurfacePoint>> result;
  wasPerfectOut = true;

  for (auto& epPtr : paths) {
    FlipEdgePath& path = *epPtr;

    // Build list of halfedges
    std::vector<Halfedge> heList = path.getHalfedgeList();

    // Trace out the halfedge
    result.emplace_back();
    std::vector<SurfacePoint>& thisResult = result.back();
    for (Halfedge he : heList) {
      std::vector<SurfacePoint> thisTrace = tri->traceIntrinsicHalfedgeAlongInput(he);

      // Check success
      SurfacePoint& lastP = thisTrace.back();
      bool thisTraceSuccess = onSameElement(lastP, tri->vertexLocations[he.twin().vertex()]);
      wasPerfectOut = wasPerfectOut && thisTraceSuccess;

      // Avoid repeat entries between consecutive traces
      if (!thisResult.empty() && onSameElement(thisResult.back(), thisTrace.front())) {
        thisResult.pop_back();
      }

      // Add the points to the list
      thisResult.insert(std::end(thisResult), std::begin(thisTrace), std::end(thisTrace));
    }
  }


  return result;
}

std::vector<std::vector<SurfacePoint>> FlipEdgeNetwork::getAllEdgePolyline() {

  std::vector<std::vector<SurfacePoint>> result;

  for (Edge e : tri->mesh.edges()) {

    // Trace out the halfedge
    result.emplace_back();
    std::vector<SurfacePoint>& thisResult = result.back();
    std::vector<SurfacePoint> thisTrace = tri->traceIntrinsicHalfedgeAlongInput(e.halfedge());

    // Add the points to the list
    thisResult.insert(std::end(thisResult), std::begin(thisTrace), std::end(thisTrace));
  }

  return result;
}

std::vector<std::vector<Vector3>> FlipEdgeNetwork::pathTo3D(const std::vector<std::vector<SurfacePoint>>& pathPoints) {
  std::vector<std::vector<Vector3>> pathTraces3D;

  if (posGeom == nullptr) {
    throw std::runtime_error("can't visualize, no position geometry registered. set the posGeom member");
    return pathTraces3D;
  }

  for (const std::vector<SurfacePoint>& edgePath : pathPoints) {
    pathTraces3D.emplace_back();
    for (const SurfacePoint& p : edgePath) {
      Vector3 p3d = p.interpolate(posGeom->vertexPositions);
      pathTraces3D.back().push_back(p3d);
    }
  }

  return pathTraces3D;
}

std::vector<std::vector<Vector3>> FlipEdgeNetwork::getPathPolyline3D() { return pathTo3D(getPathPolyline()); }

std::vector<std::vector<Vector3>> FlipEdgeNetwork::getAllEdgePolyline3D() { return pathTo3D(getAllEdgePolyline()); }

void FlipEdgeNetwork::savePathOBJLine(std::string filenamePrefix, bool withAll) {
  if (posGeom == nullptr) {
    throw std::runtime_error("can't visualize, no position geometry registered");
    return;
  }


  std::vector<std::vector<SurfacePoint>> polyline;
  if (withAll) {
    polyline = getAllEdgePolyline();
  } else {
    polyline = getPathPolyline();
  }


  // output file
  std::ofstream outFile(filenamePrefix + "lines_out.obj");
  if (!outFile) throw std::runtime_error("couldn't open");

  std::vector<std::vector<size_t>> lineInds;
  size_t iP = 0;
  for (auto& line : polyline) {
    lineInds.emplace_back();
    std::vector<size_t>& lineInd = lineInds.back();
    for (SurfacePoint& p : line) {
      Vector3 pos = p.interpolate(posGeom->vertexPositions);

      outFile << "v " << pos.x << " " << pos.y << " " << pos.z << "\n";
      lineInd.push_back(iP);
      iP++;
    }
  }

  for (auto& lineInd : lineInds) {
    outFile << "l";
    for (size_t i : lineInd) {
      outFile << " " << (i + 1);
    }
    outFile << "\n";
  }
}

bool FlipEdgeNetwork::intrinsicTriIsOriginal() {
  for (Edge e : mesh.edges()) {
    if (!tri->edgeIsOriginal[e]) {
      return false;
    }
  }
  return true;
}

void FlipEdgeNetwork::purgeStaleQueueEntries() {
  wedgeAngleQueue = std::priority_queue<WeightedAngle, std::vector<WeightedAngle>, std::greater<WeightedAngle>>();
  addAllWedgesToAngleQueue();
}


bool FlipPathSegment::isEndpoint() {
  SegmentID prevID, nextID;
  Halfedge he;
  std::tie(he, prevID, nextID) = path->pathHeInfo[id];
  return prevID == INVALID_IND || nextID == INVALID_IND;
}


FlipPathSegment FlipPathSegment::next() {
  SegmentID prevID, nextID;
  Halfedge he;
  std::tie(he, prevID, nextID) = path->pathHeInfo[id];
  return FlipPathSegment{path, nextID};
}
FlipPathSegment FlipPathSegment::prev() {
  SegmentID prevID, nextID;
  Halfedge he;
  std::tie(he, prevID, nextID) = path->pathHeInfo[id];
  return FlipPathSegment{path, prevID};
}
Halfedge FlipPathSegment::halfedge() {
  SegmentID prevID, nextID;
  Halfedge he;
  std::tie(he, prevID, nextID) = path->pathHeInfo[id];
  return he;
}

void FlipPathSegment::splitEdge(double tSplit) {
  SegmentID prevID, nextID;
  Halfedge he;
  std::tie(he, prevID, nextID) = path->pathHeInfo[id];
  Halfedge origHe = he;

  Halfedge newHeFront = path->network.tri->splitEdge(he, tSplit);
  path->network.updatePathAfterEdgeSplit(origHe, newHeFront);
}

void FlipEdgeNetwork::updatePathAfterEdgeSplit(Halfedge origHe, Halfedge newHeFront) {

  // nothing to do if no paths on this edge
  if (pathsAtEdge[origHe.edge()].empty()) {
    return;
  }

  // verify that there is just one path segment on the edge (this doesn't handle more general case yet)
  if (pathsAtEdge[origHe.edge()].size() != 1) {
    throw std::runtime_error("only tested for splitting edge with one path on it");
  }

  // Get a reference to the (one) path along the edge
  FlipPathSegment pathSeg = pathsAtEdge[origHe.edge()].front();
  FlipEdgePath* path = pathSeg.path;
  size_t id = pathSeg.id;

  // Get connectivity information for the path
  SegmentID prevID, nextID;
  Halfedge he;
  std::tie(he, prevID, nextID) = path->pathHeInfo[id];

  Halfedge newHeBack = newHeFront.prevOrbitFace().twin().prevOrbitFace();
  bool forwardDir = (he == he.edge().halfedge()); // does the segment point along the refence halfedge, or backwards?

  // ==== Diagram ====
  //
  // ... ---> [this segment] ---> [new segment] ---> ...
  // before:  --------------- he --------------->
  // after:   ----newHeBack---> ----newHeFront-->           (might be flipped if !forwardDir)
  //
  // =================

  // Create a new segment
  size_t newSegID = path->network.getNextUniquePathSegmentInd();

  // Move this segment to the prior (in terms of path orientation) of the two new halfedges, the new segment gets the
  // latter.
  Halfedge newSegHe;
  if (forwardDir) {
    he = newHeBack;
    newSegHe = newHeFront;
  } else {
    he = newHeFront.twin();
    newSegHe = newHeBack.twin();
  }
  std::get<0>(path->pathHeInfo[id]) = he;


  // Fix up the connectivity
  size_t newSegNext = nextID;
  size_t newSegPrev = id;
  std::get<2>(path->pathHeInfo[id]) = newSegID; // set next()
  if (nextID != INVALID_IND) {
    std::get<1>(path->pathHeInfo[nextID]) = newSegID; // set next().prev()
  }


  // Add to auxiliary data structures
  path->pathHeInfo[newSegID] = std::tie(newSegHe, newSegPrev, newSegNext);
  FlipPathSegment newSeg{path, newSegID};

  popOutsideSegment(origHe);
  pushOutsideSegment(he, pathSeg);
  pushOutsideSegment(newSegHe, newSeg);


  addToWedgeAngleQueue(pathSeg); // spliting the edge probably changed computed angle by numerical epsilon
  // NOTE: shouldn't need to add to wedge angle queue for straightening, the newly created wedge is straight, since
  // it's a subdivided (geodesic) edge, but do so anyway to think less about numerics
  addToWedgeAngleQueue(newSeg);

  // validate();
}

// Inherit operators from tuple
bool FlipPathSegment::operator==(const FlipPathSegment& other) const {
  return std::make_tuple(path, id) == std::make_tuple(other.path, other.id);
}
bool FlipPathSegment::operator!=(const FlipPathSegment& other) const {
  return std::make_tuple(path, id) != std::make_tuple(other.path, other.id);
}
bool FlipPathSegment::operator>(const FlipPathSegment& other) const {
  return std::make_tuple(path, id) > std::make_tuple(other.path, other.id);
}
bool FlipPathSegment::operator>=(const FlipPathSegment& other) const {
  return std::make_tuple(path, id) >= std::make_tuple(other.path, other.id);
}
bool FlipPathSegment::operator<(const FlipPathSegment& other) const {
  return std::make_tuple(path, id) < std::make_tuple(other.path, other.id);
}
bool FlipPathSegment::operator<=(const FlipPathSegment& other) const {
  return std::make_tuple(path, id) <= std::make_tuple(other.path, other.id);
}

} // namespace surface
} // namespace geometrycentral
