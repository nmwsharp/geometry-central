#include "geometrycentral/surface/signpost_intrinsic_triangulation.h"

#include "geometrycentral/surface/barycentric_coordinate_helpers.h"
#include "geometrycentral/surface/exact_polyhedral_geodesics.h"
#include "geometrycentral/surface/mesh_graph_algorithms.h"
#include "geometrycentral/surface/trace_geodesic.h"

#include <iomanip>
#include <queue>

using std::cout;
using std::endl;

namespace geometrycentral {
namespace surface {


SignpostIntrinsicTriangulation::SignpostIntrinsicTriangulation(ManifoldSurfaceMesh& mesh_,
                                                               IntrinsicGeometryInterface& inputGeom_)
    // Note: this initializer list does something slightly wacky: it creates the new mesh on the heap, then loses track
    // of pointer while setting the BaseGeometryInterface::mesh reference to it. Later, it picks the pointer back up
    // from the reference and wraps it in the intrinsicMesh unique_ptr<>. I believe that this is all valid, but its
    // probably a sign of bad design.
    : IntrinsicGeometryInterface(*mesh_.copy().release()), inputMesh(mesh_), inputGeom(inputGeom_),
      intrinsicMesh(dynamic_cast<ManifoldSurfaceMesh*>(&mesh)) {

  // Make sure the input mesh is triangular
  if (!mesh.isTriangular()) {
    throw std::runtime_error("signpost triangulation requires triangle mesh as input");
  }

  // == Initialize geometric data
  inputGeom.requireEdgeLengths();
  inputGeom.requireHalfedgeVectorsInVertex();
  inputGeom.requireVertexAngleSums();
  inputGeom.requireCornerAngles();

  // Just copy lengths
  intrinsicEdgeLengths = inputGeom.edgeLengths.reinterpretTo(mesh);

  // Prepare directions and angle sums
  intrinsicHalfedgeDirections = HalfedgeData<double>(mesh);
  intrinsicVertexAngleSums = VertexData<double>(mesh);

  // Walk around the vertex, constructing angular directions
  requireCornerAngles();
  for (Vertex v : mesh.vertices()) {
    double runningAngle = 0.;
    Halfedge firstHe = v.halfedge();
    Halfedge currHe = firstHe;
    do {
      double cornerAngle = cornerAngles[currHe.corner()];
      intrinsicHalfedgeDirections[currHe] = runningAngle;
      runningAngle += cornerAngle;

      if (!currHe.isInterior()) {
        break;
      }
      currHe = currHe.next().next().twin();
    } while (currHe != firstHe);

    intrinsicVertexAngleSums[v] = runningAngle;
  }

  // Initialize vertex locations
  vertexLocations = VertexData<SurfacePoint>(mesh);
  for (size_t iV = 0; iV < mesh.nVertices(); iV++) {
    vertexLocations[iV] = SurfacePoint(inputMesh.vertex(iV));
  }

  // Initialize all edges as original, but new ones should be false
  edgeIsOriginal = EdgeData<bool>(mesh, false);
  edgeIsOriginal.fill(true);

  requireHalfedgeVectorsInVertex();
  requireHalfedgeVectorsInFace();
  requireVertexAngleSums();

  // == Register the default callback which maintains marked edges
  auto updateMarkedEdges = [&](Edge oldE, Halfedge newHe1, Halfedge newHe2) {
    if (markedEdges.size() > 0 && markedEdges[oldE]) {
      markedEdges[newHe1.edge()] = true;
      markedEdges[newHe2.edge()] = true;
    }
  };
  edgeSplitCallbackList.push_back(updateMarkedEdges);
}

void SignpostIntrinsicTriangulation::setMarkedEdges(const EdgeData<bool>& markedEdges_) {
  markedEdges = markedEdges_;
  markedEdges.setDefault(false);
}


std::vector<SurfacePoint> SignpostIntrinsicTriangulation::traceHalfedge(Halfedge he, bool trimEnd) {

  // Optimization: don't both tracing original edges, just report them directly
  if (edgeIsOriginal[he.edge()]) {
    if (vertexLocations[he.vertex()].type != SurfacePointType::Vertex ||
        vertexLocations[he.twin().vertex()].type != SurfacePointType::Vertex) {
      throw std::runtime_error("edgeIsOriginal cache is out of date");
    }
    Vertex vA = vertexLocations[he.vertex()].vertex;
    Vertex vB = vertexLocations[he.twin().vertex()].vertex;
    std::vector<SurfacePoint> result{SurfacePoint(vA), SurfacePoint(vB)};
    return result;
  }

  // Gather values to trace
  SurfacePoint startP = vertexLocations[he.vertex()];
  Vector2 traceVec = halfedgeVector(he);


  // Do the actual tracing
  TraceOptions options;
  options.includePath = true;
  options.maxIters = mesh.nFaces() * 10;
  TraceGeodesicResult result = traceGeodesic(inputGeom, startP, traceVec, options);

  // Trim off end crumbs if applicable
  Vertex endVert = he.twin().vertex();
  if (trimEnd && vertexLocations[endVert].type == SurfacePointType::Vertex) {
    bool success = trimTraceResult(result, endVert);
    if (success) {
      // Append the endpoint
      result.pathPoints.push_back(vertexLocations[endVert]);
    } else {
      // If trimming failed (because the trace didn't even hit the 1-ring of target), just stick with whatever we go
      // initially
      result = traceGeodesic(inputGeom, startP, traceVec, options);
    }
  }

  return result.pathPoints;
}

EdgeData<std::vector<SurfacePoint>> SignpostIntrinsicTriangulation::traceEdges() {

  EdgeData<std::vector<SurfacePoint>> tracedEdges(mesh);

  for (Edge e : mesh.edges()) {
    Halfedge he = e.halfedge();
    tracedEdges[e] = traceHalfedge(he, false);
  }

  return tracedEdges;
}


// ======================================================
// ======== Queries & Accessors
// ======================================================


bool SignpostIntrinsicTriangulation::isDelaunay(Edge e) {
  if (!isFixed(e) && edgeCotanWeight(e) < -delaunayEPS) {
    return false;
  }
  return true;
}
bool SignpostIntrinsicTriangulation::isDelaunay() {
  for (Edge e : mesh.edges()) {
    if (!isDelaunay(e)) {
      return false;
    }
  }
  return true;
}

double SignpostIntrinsicTriangulation::minAngleDegrees() {
  double minAngle = std::numeric_limits<double>::infinity();
  for (Corner c : mesh.corners()) {
    minAngle = std::min(minAngle, cornerAngle(c));
  }
  return minAngle * 180. / M_PI;
}

// ======================================================
// ======== Mutators
// ======================================================

bool SignpostIntrinsicTriangulation::flipEdgeIfNotDelaunay(Edge e) {

  // Can't flip
  if (isFixed(e)) return false;

  // Don't want to flip
  double cWeight = edgeCotanWeight(e);
  if (cWeight > -delaunayEPS) return false;

  // Get geometric data
  Halfedge he = e.halfedge();
  std::array<Vector2, 4> layoutPositions = layoutDiamond(he);

  // Combinatorial flip
  bool flipped = intrinsicMesh->flip(e, false);

  // Should always be possible, something unusual is going on if we end up here
  if (!flipped) {
    return false;
  }

  // Compute the new edge length
  double newLength = (layoutPositions[1] - layoutPositions[3]).norm();

  // If we're going to create a non-finite edge length, abort the flip
  // (only happens if you're in a bad numerical place)
  if (!std::isfinite(newLength)) {
    intrinsicMesh->flip(e, false);
    return false;
  }

  // Assign the new edge lengths
  intrinsicEdgeLengths[e] = newLength;
  edgeLengths[e] = newLength;

  // Update edge angles
  updateAngleFromCWNeighor(e.halfedge());
  updateAngleFromCWNeighor(e.halfedge().twin());
  updateFaceBasis(e.halfedge().face());
  updateFaceBasis(e.halfedge().twin().face());

  edgeIsOriginal[e] = false;

  invokeEdgeFlipCallbacks(e);
  return true;
}

bool SignpostIntrinsicTriangulation::flipEdgeIfPossible(Edge e, double possibleEPS) {

  // Can't flip
  if (isFixed(e)) return false;

  // Get geometric data
  Halfedge he = e.halfedge();
  std::array<Vector2, 4> layoutPositions = layoutDiamond(he);

  // Test if geometryically flippable flippable (both signed areas of new triangles are positive)
  double A1 = cross(layoutPositions[1] - layoutPositions[0], layoutPositions[3] - layoutPositions[0]);
  double A2 = cross(layoutPositions[3] - layoutPositions[2], layoutPositions[1] - layoutPositions[2]);
  double areaEPS = possibleEPS * (A1 + A2);
  if (A1 < areaEPS || A2 < areaEPS) {
    return false;
  }


  // Combinatorial flip
  bool flipped = intrinsicMesh->flip(e, false);

  // Might not have been flippable for connectivity reasons
  if (!flipped) {
    return false;
  }

  // Compute the new edge length
  double newLength = (layoutPositions[1] - layoutPositions[3]).norm();

  // If we're going to create a non-finite edge length, abort the flip
  // (only happens if you're in a bad numerical place)
  if (!std::isfinite(newLength)) {
    intrinsicMesh->flip(e, false);
    return false;
  }

  // Assign the new edge lengths
  // TODO project to satisfy triangle inequality?
  intrinsicEdgeLengths[e] = newLength;
  edgeLengths[e] = newLength;

  // Update edge angles
  updateAngleFromCWNeighor(e.halfedge());
  updateAngleFromCWNeighor(e.halfedge().twin());
  updateFaceBasis(e.halfedge().face());
  updateFaceBasis(e.halfedge().twin().face());

  edgeIsOriginal[e] = false;

  invokeEdgeFlipCallbacks(e);
  return true;
}

void SignpostIntrinsicTriangulation::flipEdgeManual(Edge e, double newLength, double forwardAngle, double reverseAngle,
                                                    bool isOrig, bool reverseFlip) {

  int flipCount = reverseFlip ? 3 : 1; // three flips give opposite orientaiton
  for (int i = 0; i < flipCount; i++) {
    bool flipped = intrinsicMesh->flip(e, false);
    if (!flipped) throw std::runtime_error("could not perform manual flip");
  }

  intrinsicEdgeLengths[e] = newLength;
  edgeLengths[e] = newLength;

  // Update other derived geometric data
  intrinsicHalfedgeDirections[e.halfedge()] = forwardAngle;
  intrinsicHalfedgeDirections[e.halfedge().twin()] = reverseAngle;
  halfedgeVectorsInVertex[e.halfedge()] = halfedgeVector(e.halfedge());
  halfedgeVectorsInVertex[e.halfedge().twin()] = halfedgeVector(e.halfedge().twin());
  updateFaceBasis(e.halfedge().face());
  updateFaceBasis(e.halfedge().twin().face());

  edgeIsOriginal[e] = isOrig;
  invokeEdgeFlipCallbacks(e);
}

Vertex SignpostIntrinsicTriangulation::insertVertex(SurfacePoint newPositionOnIntrinsic) {
  switch (newPositionOnIntrinsic.type) {
  case SurfacePointType::Vertex: {
    throw std::logic_error("can't insert vertex at vertex");
    break;
  }
  case SurfacePointType::Edge: {
    return insertVertex_edge(newPositionOnIntrinsic).vertex();
    break;
  }
  case SurfacePointType::Face: {
    return insertVertex_face(newPositionOnIntrinsic);
    break;
  }
  }
  return Vertex();
}

Halfedge SignpostIntrinsicTriangulation::insertVertex_edge(SurfacePoint newP) {

  // === (1) Gather some data about the edge we're about to insert into

  Edge insertionEdge = newP.edge;
  Face fA = insertionEdge.halfedge().face();
  Face fB = insertionEdge.halfedge().twin().face();
  bool isOnBoundary = fB.isBoundaryLoop();

  // Find coordinates in (both) faces and compute the lengths of the new wedges
  double backLen, frontLen, Alen, Blen;

  // in A
  backLen = newP.tEdge * intrinsicEdgeLengths[insertionEdge];
  frontLen = (1. - newP.tEdge) * intrinsicEdgeLengths[insertionEdge];

  int iA = halfedgeIndexInTriangle(insertionEdge.halfedge());
  std::array<Vector2, 3> vertCoords = vertexCoordinatesInTriangle(fA);
  Vector2 posA = (1. - newP.tEdge) * vertCoords[iA] + newP.tEdge * vertCoords[(iA + 1) % 3];
  Alen = (posA - vertCoords[(iA + 2) % 3]).norm();


  if (!isOnBoundary) { // in B
    // WARNING: these code paths are not as well-tested, since they don't happen in the common insert-along-boundary
    // case
    int iB = halfedgeIndexInTriangle(insertionEdge.halfedge().twin());
    std::array<Vector2, 3> vertCoords = vertexCoordinatesInTriangle(fB);
    Vector2 posB = newP.tEdge * vertCoords[iB] + (1. - newP.tEdge) * vertCoords[(iB + 1) % 3];
    Blen = (posB - vertCoords[(iB + 2) % 3]).norm();
  } else {
    Blen = -777;
  }


  // === (2) Insert vertex

  // Put a new vertex inside of the proper intrinsic face
  Halfedge newHeFront = intrinsicMesh->splitEdgeTriangular(insertionEdge);
  edgeIsOriginal[insertionEdge] = false;
  Vertex newV = newHeFront.vertex();

  // = Update data arrays for the new vertex
  if (isOnBoundary) {
    intrinsicVertexAngleSums[newV] = M_PI;
    vertexAngleSums[newV] = M_PI;
  } else {
    intrinsicVertexAngleSums[newV] = 2. * M_PI;
    vertexAngleSums[newV] = 2. * M_PI;
  }


  // == (3) Assign edge lengths to the new edges
  Halfedge currHe = newHeFront;
  Halfedge newHeBack;
  std::array<double, 4> newLens = {frontLen, Alen, backLen, Blen};
  for (int i = 0; i < (isOnBoundary ? 3 : 4); i++) {
    intrinsicEdgeLengths[currHe.edge()] = newLens[i];
    edgeLengths[currHe.edge()] = newLens[i];
    if (i == 2) newHeBack = currHe;
    currHe = currHe.next().next().twin();
  }

  // === (4) Now that we have edge lengths, sort out tangent spaces and position on supporting.
  resolveNewVertex(newV, newP);

  invokeEdgeSplitCallbacks(insertionEdge, newHeFront, newHeBack);

  return newHeFront;
}

Vertex SignpostIntrinsicTriangulation::insertVertex_face(SurfacePoint newP) {

  // === (1) Gather some data about the face we're about to insert into
  Face insertionFace = newP.face;
  std::array<Vector2, 3> vertCoords = vertexCoordinatesInTriangle(insertionFace);
  Vector2 newPCoord = (newP.faceCoords[1] * vertCoords[1] + newP.faceCoords[2] * vertCoords[2]);
  std::array<double, 3> newEdgeLengths;
  std::array<Halfedge, 3> oldFaceHalfedges;
  size_t i = 0;
  for (Halfedge he : insertionFace.adjacentHalfedges()) {
    newEdgeLengths[i] = (newPCoord - vertCoords[i]).norm();
    if (!std::isfinite(newEdgeLengths[i])) {
      throw std::runtime_error("non finite edge length");
    }
    oldFaceHalfedges[i] = he;
    i++;
  }


  // === (2) Insert vertex

  // Put a new vertex inside of the proper intrinsic face
  Vertex newV = intrinsicMesh->insertVertex(insertionFace);

  // = Update data arrays for the new vertex
  intrinsicVertexAngleSums[newV] = 2. * M_PI;
  vertexAngleSums[newV] = 2. * M_PI;


  // == (3) Assign edge lengths to the new edges

  // Set edge lengths first by looking for the proper new edge
  for (size_t j = 0; j < 3; j++) {
    double thisLen = newEdgeLengths[j];
    Halfedge origHe = oldFaceHalfedges[j];

    // Find the new edge which this length belongs to
    for (Halfedge heV : newV.outgoingHalfedges()) {
      if (heV.next() == origHe) {
        intrinsicEdgeLengths[heV.edge()] = thisLen;
        edgeLengths[heV.edge()] = thisLen;
      }
    }
  }

  // === (4) Now that we have edge lengths, sort out tangent spaces and position on supporting.
  resolveNewVertex(newV, newP);

  invokeFaceInsertionCallbacks(insertionFace, newV);
  return newV;
}

Vertex SignpostIntrinsicTriangulation::insertCircumcenter(Face f) {

  // === Circumcenter in barycentric coordinates

  Halfedge he0 = f.halfedge();
  double a = intrinsicEdgeLengths[he0.next().edge()];
  double b = intrinsicEdgeLengths[he0.next().next().edge()];
  double c = intrinsicEdgeLengths[he0.edge()];
  double a2 = a * a;
  double b2 = b * b;
  double c2 = c * c;
  Vector3 circumcenterLoc = {a2 * (b2 + c2 - a2), b2 * (c2 + a2 - b2), c2 * (a2 + b2 - c2)};
  circumcenterLoc = normalizeBarycentric(circumcenterLoc);

  // Trace from the barycenter (have to trace from somewhere)
  Vector3 barycenter = Vector3::constant(1. / 3.);
  Vector3 vecToCircumcenter = circumcenterLoc - barycenter;

  // === Trace the ray to find the location of the new point on the intrinsic meshes

  // Data we need from the intrinsic trace
  TraceOptions options;
  if (markedEdges.size() > 0) {
    options.barrierEdges = &markedEdges;
  }
  TraceGeodesicResult intrinsicTraceResult = traceGeodesic(*this, f, barycenter, vecToCircumcenter, options);
  // intrinsicTracer->snapEndToEdgeIfClose(intrinsicCrumbs); TODO
  // SurfacePoint newPositionOnIntrinsic = intrinsicTraceResult.endPoint.inSomeFace();
  SurfacePoint newPositionOnIntrinsic = intrinsicTraceResult.endPoint;

  // If the circumcenter is blocked by an edge, insert the midpoint of that edge instead
  // (which happens to be just want is needed for Chew's 2nd algo).
  if (newPositionOnIntrinsic.type == SurfacePointType::Edge) {
    newPositionOnIntrinsic.tEdge = 0.5;
  }

  // === Phase 3: Add the new vertex
  return insertVertex(newPositionOnIntrinsic);
}

Vertex SignpostIntrinsicTriangulation::insertBarycenter(Face f) {
  SurfacePoint barycenterOnIntrinsic(f, Vector3::constant(1. / 3.));
  return insertVertex(barycenterOnIntrinsic);
}

Face SignpostIntrinsicTriangulation::removeInsertedVertex(Vertex v) {
  // Strategy: flip edges until the vertex has degree three, then remove by replacing with a single face
  // TODO needs a proof that this always works... what about self edges, etc? Seems to work well.

  // What about starting with degree < 3? Since this vertex necessarily has angle sum 2PI, this could only happen in the
  // case of degree 2, with exactly degenerate triangles. Since we assume non-degenerate triangles throughout, we'll
  // consider that to not happen.

  if (vertexLocations[v].type == SurfacePointType::Vertex) return Face(); // can't remove original vertices

  if (isOnFixedEdge(v)) {
    return Face(); // don't try to remove boundary vertices, for now at least
  }

  // Flip edges until
  size_t iterCount = 0;
  while (v.degree() != 3) {

    // Find any edge we can flip
    bool anyFlipped = false;
    for (Edge e : v.adjacentEdges()) {
      anyFlipped = flipEdgeIfPossible(e);
      if (anyFlipped) break;
    }

    // failsafe, in case we get numerically stuck, or there are too many fixed edges (or the algorithm is broken)
    if (!anyFlipped || iterCount > 10 * v.degree()) {
      return Face();
    }

    iterCount++;
  }

  // give up if something went wrong (eg. flipped edges)
  if (v.degree() != 3) return Face();

  // Remove the vertex
  Face newF = intrinsicMesh->removeVertex(v);
  updateFaceBasis(newF);
  return newF;
}

Halfedge SignpostIntrinsicTriangulation::splitEdge(Halfedge he, double tSplit) {
  return insertVertex_edge(SurfacePoint(he, tSplit));
}

void SignpostIntrinsicTriangulation::flipToDelaunay() {

  std::deque<Edge> edgesToCheck;
  EdgeData<bool> inQueue(mesh, true);
  for (Edge e : mesh.edges()) {
    edgesToCheck.push_back(e);
  }

  size_t nFlips = 0;
  while (!edgesToCheck.empty()) {

    // Get the top element from the queue of possibily non-Delaunay edges
    Edge e = edgesToCheck.front();
    edgesToCheck.pop_front();
    inQueue[e] = false;

    bool wasFlipped = flipEdgeIfNotDelaunay(e);

    if (!wasFlipped) continue;

    // Handle the aftermath of a flip
    nFlips++;

    // Add neighbors to queue, as they may need flipping now
    Halfedge he = e.halfedge();
    Halfedge heN = he.next();
    Halfedge heT = he.twin();
    Halfedge heTN = heT.next();
    std::vector<Edge> neighEdges = {heN.edge(), heN.next().edge(), heTN.edge(), heTN.next().edge()};
    for (Edge nE : neighEdges) {
      if (!inQueue[nE]) {
        edgesToCheck.push_back(nE);
        inQueue[nE] = true;
      }
    }
  }

  refreshQuantities();
}

void SignpostIntrinsicTriangulation::delaunayRefine(double angleThreshDegrees, double circumradiusThresh,
                                                    size_t maxInsertions) {

  // Relationship between angles and circumradius-to-edge
  double angleThreshRad = angleThreshDegrees * M_PI / 180.;
  double circumradiusEdgeRatioThresh = 1.0 / (2.0 * std::sin(angleThreshRad));

  // Build a function to test if a face violates the circumradius ratio condition
  auto needsCircumcenterRefinement = [&](Face f) {
    double c = circumradius(f);
    double l = shortestEdge(f);

    bool needsRefinementLength = c > circumradiusThresh;

    // Explicit check allows us to skip degree one vertices (can't make those angles smaller!)
    bool needsRefinementAngle = false;
    for (Halfedge he : f.adjacentHalfedges()) {

      double baseAngle = cornerAngle(he.corner());
      if (baseAngle < angleThreshRad) {

        // If it's already a degree one vertex, nothing we can do here
        bool isDegreeOneVertex = he.next().next() == he.twin();
        if (isDegreeOneVertex) {
          continue;
        }

        // If it's a fixed corner, can't make it smaller
        if (isFixed(he.edge()) && isFixed(he.prevOrbitFace().edge())) {
          continue;
        }

        needsRefinementAngle = true;
      }
    }

    return needsRefinementAngle || needsRefinementLength;
  };

  // Call the general version
  delaunayRefine(needsCircumcenterRefinement, maxInsertions);
}


void SignpostIntrinsicTriangulation::delaunayRefine(const std::function<bool(Face)>& shouldRefine,
                                                    size_t maxInsertions) {

  // Manages a check at the bottom to avoid infinite-looping when numerical baddness happens
  int recheckCount = 0;
  const int MAX_RECHECK_COUNT = 5;

  // Track statistics
  size_t nFlips = 0;
  size_t nInsertions = 0;

  // Initialize queue of (possibly) non-delaunay edges
  std::deque<Edge> delaunayCheckQueue;
  EdgeData<bool> inDelaunayQueue(mesh, false);
  for (Edge e : mesh.edges()) {
    delaunayCheckQueue.push_back(e);
    inDelaunayQueue[e] = true;
  }


  // Return a weight to use for sorting PQ. Usually sorts by biggest area, but also puts faces on boundary first with
  // weight inf.
  auto areaWeight = [&](Face f) {
    for (Edge e : f.adjacentEdges()) {
      if (isFixed(e)) return std::numeric_limits<double>::infinity();
    }
    return area(f);
  };

  // Initialize queue of (possibly) circumradius-violating faces, processing the largest faces first (good heuristic)
  typedef std::pair<double, Face> AreaFace;
  std::priority_queue<AreaFace, std::vector<AreaFace>, std::less<AreaFace>> circumradiusCheckQueue;
  for (Face f : mesh.faces()) {
    if (shouldRefine(f)) {
      circumradiusCheckQueue.push(std::make_pair(areaWeight(f), f));
    }
  }

  // Register a callback which checks the neighbors of an edge for further processing after a flip. It's useful to use a
  // callback, rather than just checking in the loop, because other internal subroutines might perform flips. In
  // particular, removeInsertedVertex() currently performs flips internally, which might trigger updates.
  auto checkNeighborsAfterFlip = [&](Edge e) {
    // std::cout << "  flipped edge " << e << std::endl;
    nFlips++;

    // Add neighboring faces, which might violate circumradius constraint
    std::vector<Face> neighFaces = {e.halfedge().face(), e.halfedge().twin().face()};
    for (Face nF : neighFaces) {
      if (shouldRefine(nF)) {
        circumradiusCheckQueue.push(std::make_pair(areaWeight(nF), nF));
      }
    }

    // Add neighbors to queue, as they may need flipping now
    Halfedge he = e.halfedge();
    Halfedge heN = he.next();
    Halfedge heT = he.twin();
    Halfedge heTN = heT.next();
    std::vector<Edge> neighEdges = {heN.edge(), heN.next().edge(), heTN.edge(), heTN.next().edge()};
    for (Edge nE : neighEdges) {
      if (!inDelaunayQueue[nE]) {
        delaunayCheckQueue.push_back(nE);
        inDelaunayQueue[nE] = true;
      }
    }
  };
  auto flipCallbackHandle = edgeFlipCallbackList.insert(std::end(edgeFlipCallbackList), checkNeighborsAfterFlip);

  // Register a callback, which will be invoked to delete previously-inserted vertices whenever refinment splits an edge
  auto deleteNearbyVertices = [&](Edge e, Halfedge he1, Halfedge he2) {
    // radius of the diametral ball
    double ballRad = std::max(intrinsicEdgeLengths[he1.edge()], intrinsicEdgeLengths[he2.edge()]);
    Vertex newV = he1.vertex();

    // Find all vertices within range.
    // Most properly, this should probably be a polyhedral geodesic ball search, but that creates a dependence on
    // polyhedral shortest paths which is bad for performance and robustness. Using a Dijsktra ball instead seems to
    // work fine in practice, and I think you could argue that the factor of 2 makes it provably correct, due to the
    // stretch factor of a Delaunay triangulation (recalling that deleting extra interior inserted vertices does not
    // affect correctness).
    //
    // std::unordered_map<Vertex, double> nearbyVerts = vertexGeodesicDistanceWithinRadius(*this, newV, ballRad);
    std::unordered_map<Vertex, double> nearbyVerts = vertexDijkstraDistanceWithinRadius(*this, newV, 2. * ballRad);

    // remove inserted vertices
    for (auto p : nearbyVerts) {
      Vertex v = p.first;
      if (v != newV && !isOnFixedEdge(v) && vertexLocations[v].type != SurfacePointType::Vertex) {
        // std::cout << "  removing inserted vertex " << v << std::endl;
        Face fReplace = removeInsertedVertex(v);

        if (fReplace != Face()) {

          // Add adjacent edges for Delaunay check
          for (Edge nE : fReplace.adjacentEdges()) {
            if (!inDelaunayQueue[nE]) {
              delaunayCheckQueue.push_back(nE);
              inDelaunayQueue[nE] = true;
            }
          }

          // Add face for refine check
          if (shouldRefine(fReplace)) {
            circumradiusCheckQueue.push(std::make_pair(areaWeight(fReplace), fReplace));
          }
        }
      }
    }
  };

  // Add our new callback at the end, so it gets invoked after any user-defined callbacks which ought to get called
  // right after the split before we mess with the mesh.
  auto splitCallbackHandle = edgeSplitCallbackList.insert(std::end(edgeSplitCallbackList), deleteNearbyVertices);

  // === Outer iteration: flip and insert until we have a mesh that satisfies both angle and circumradius goals
  do {

    // == First, flip to delaunay
    while (!delaunayCheckQueue.empty()) {

      // Get the top element from the queue of possibily non-Delaunay edges
      Edge e = delaunayCheckQueue.front();
      delaunayCheckQueue.pop_front();
      if (e.isDead()) continue;
      inDelaunayQueue[e] = false;

      flipEdgeIfNotDelaunay(e);

      // Remember that up we registered a callback up above which checks neighbors for subsequent processing after a
      // flip.
    }

    // == Second, insert one circumcenter

    // If we've already inserted the max number of points, call it a day
    if (maxInsertions != INVALID_IND && nInsertions == maxInsertions) {
      break;
    }

    // Try to insert just one circumcenter
    if (!circumradiusCheckQueue.empty()) {

      // Get the biggest face
      Face f = circumradiusCheckQueue.top().second;
      double A = circumradiusCheckQueue.top().first;
      circumradiusCheckQueue.pop();
      if (f.isDead()) continue;

      // Two things might have changed that would cause us to skip this entry:
      //   -If the area has changed since this face was inserted in to the queue, skip it. Note that we don't need to
      //    re-add it, because it must have been placed in the queue when its area was changed
      //   - This face might have been flipped to no longer violate constraint
      if (A == areaWeight(f) && shouldRefine(f)) {

        // std::cout << "  refining face " << f << std::endl;
        Vertex newVert = insertCircumcenter(f);
        nInsertions++;

        // Mark everything in the 1-ring as possibly non-Delaunay and possibly violating the circumradius constraint
        for (Face nF : newVert.adjacentFaces()) {

          // Check circumradius constraint
          if (shouldRefine(nF)) {
            circumradiusCheckQueue.push(std::make_pair(areaWeight(nF), nF));
          }

          // Check delaunay constraint
          for (Edge nE : nF.adjacentEdges()) {
            if (!inDelaunayQueue[nE]) {
              delaunayCheckQueue.push_back(nE);
              inDelaunayQueue[nE] = true;
            }
          }
        }
      }

      continue;
    }

    // If the circumradius queue is empty, make sure we didn't miss anything (can happen rarely due to numerics)
    // (but don't do this more than a few times, to avoid getting stuck in an infinite loop when numerical ultra-badness
    // happens)
    if (recheckCount < MAX_RECHECK_COUNT) {
      recheckCount++;
      bool anyFound = false;
      if (delaunayCheckQueue.empty() && circumradiusCheckQueue.empty()) {
        for (Face f : mesh.faces()) {
          if (shouldRefine(f)) {
            circumradiusCheckQueue.push(std::make_pair(areaWeight(f), f));
            anyFound = true;
          }
        }
        for (Edge e : mesh.edges()) {
          if (!isDelaunay(e)) {
            delaunayCheckQueue.push_back(e);
            inDelaunayQueue[e] = true;
            anyFound = true;
          }
        }
      }

      if (!anyFound) {
        // makes sure we don't recheck multiple times in a row
        break;
      }
    }

  } while (!delaunayCheckQueue.empty() || !circumradiusCheckQueue.empty() || recheckCount < MAX_RECHECK_COUNT);

  refreshQuantities();

  edgeSplitCallbackList.erase(splitCallbackHandle);
  edgeFlipCallbackList.erase(flipCallbackHandle);
}


void SignpostIntrinsicTriangulation::splitBentEdges(EmbeddedGeometryInterface& posGeom, double angleThreshDeg,
                                                    double relativeLengthEPS, size_t maxInsertions) {

  posGeom.requireVertexPositions();

  // == Process parameters

  // Compute shape length scale as bounding box diagonal
  Vector3 minP = Vector3::constant(std::numeric_limits<double>::infinity());
  Vector3 maxP = Vector3::constant(-std::numeric_limits<double>::infinity());
  for (Vertex v : posGeom.mesh.vertices()) {
    Vector3 p = posGeom.vertexPositions[v];
    minP = componentwiseMin(minP, p);
    maxP = componentwiseMax(maxP, p);
  }
  double lengthScale = norm(minP - maxP);
  double lengthEPS = lengthScale * relativeLengthEPS;

  double angleThresh = angleThreshDeg * PI / 180.;


  // === Make repeated passes through, splitting edges until no more bent edges remain
  bool anySplit = true;
  EdgeData<bool> edgeIsGood(mesh, false);
  size_t nSplit = 0;
  while (anySplit) {
    anySplit = false;
    for (Edge e : mesh.edges()) {

      if (maxInsertions != INVALID_IND && nSplit >= maxInsertions) break;

      if (edgeIsGood[e]) continue; // ONEDAY use a queue instead

      // Trace the edge
      std::vector<SurfacePoint> surfacePoints = traceHalfedge(e.halfedge(), false);

      // Detect the first sharp enough bend
      double tSplit = -1;
      double runningLen = 0.;
      for (size_t iP = 1; (iP + 1) < surfacePoints.size(); iP++) {
        SurfacePoint& prevS = surfacePoints[iP - 1];
        SurfacePoint& currS = surfacePoints[iP];
        SurfacePoint& nextS = surfacePoints[iP + 1];

        Vector3 prevP = prevS.interpolate(posGeom.vertexPositions);
        Vector3 currP = currS.interpolate(posGeom.vertexPositions);
        Vector3 nextP = nextS.interpolate(posGeom.vertexPositions);

        double lenFirst = (prevP - currP).norm();
        double lenSecond = (currP - nextP).norm();

        runningLen += lenFirst;

        // Skip if either segment is too short
        if (lenFirst < lengthEPS || lenSecond < lengthEPS) continue;

        // Measure the angle
        double angleBetween = angle(currP - prevP, nextP - currP);

        // Split if angle is too sharp
        if (angleBetween > angleThresh) {
          double thisTSplit = runningLen / intrinsicEdgeLengths[e];

          if (thisTSplit > relativeLengthEPS && thisTSplit < 1. - relativeLengthEPS) {
            tSplit = thisTSplit;
            break;
          }
        }
      }

      // If a bend was found, split the edge
      if (tSplit == -1) {
        edgeIsGood[e] = true;
      } else {
        anySplit = true;
        nSplit++;
        splitEdge(e.halfedge(), tSplit);
      }
    }
  }

  refreshQuantities();

  posGeom.unrequireVertexPositions();
}

// ======================================================
// ======== Geometry and Helpers
// ======================================================

void SignpostIntrinsicTriangulation::computeEdgeLengths() { edgeLengths = intrinsicEdgeLengths; }

void SignpostIntrinsicTriangulation::computeHalfedgeVectorsInVertex() {
  halfedgeVectorsInVertex = HalfedgeData<Vector2>(mesh);

  for (Halfedge he : mesh.halfedges()) {
    Vector2 traceVec = halfedgeVector(he);
    halfedgeVectorsInVertex[he] = traceVec;
  }
}

void SignpostIntrinsicTriangulation::updateAngleFromCWNeighor(Halfedge he) {

  // Handle boundary cases
  // NOTE: This makes sense because we preserve the invariant that intrinsic boundary vertices are always located along
  // the boundary of the original mesh, which has the convention that v.halfedge() begins a ccw arc along the interior.
  if (!he.isInterior()) {
    intrinsicHalfedgeDirections[he] = intrinsicVertexAngleSums[he.vertex()]; // last angle in boundary wedge
    halfedgeVectorsInVertex[he] = halfedgeVector(he);
    return;
  }
  if (!he.twin().isInterior()) {
    intrinsicHalfedgeDirections[he] = 0.; // first angle in boundary wedge
    halfedgeVectorsInVertex[he] = halfedgeVector(he);
    return;
  }

  // Get neighbor angle
  Halfedge cwHe = he.twin().next();
  double neighAngle = intrinsicHalfedgeDirections[cwHe];

  // Compute corner angle in between
  double cAngle = cornerAngle(cwHe.corner());

  // Set the updated angle
  double updatedAngle = standardizeAngle(he.vertex(), neighAngle + cAngle);
  intrinsicHalfedgeDirections[he] = updatedAngle;
  halfedgeVectorsInVertex[he] = halfedgeVector(he);
}

void SignpostIntrinsicTriangulation::updateFaceBasis(Face f) {
  Halfedge he = f.halfedge();
  double a = intrinsicEdgeLengths[he.edge()];
  he = he.next();
  double b = intrinsicEdgeLengths[he.edge()];
  he = he.next();
  double c = intrinsicEdgeLengths[he.edge()];

  Vector2 p0{0., 0.};
  Vector2 p1{a, 0.};
  Vector2 p2 = layoutTriangleVertex(p0, p1, b, c);

  he = f.halfedge();
  halfedgeVectorsInFace[he] = p1 - p0;
  he = he.next();
  halfedgeVectorsInFace[he] = p2 - p1;
  he = he.next();
  halfedgeVectorsInFace[he] = p0 - p2;
}


void SignpostIntrinsicTriangulation::resolveNewVertex(Vertex newV, SurfacePoint intrinsicPoint) {

  // == (1) Compute angular coordinates for the halfedges
  // Now that we have valid edge lengths, compute halfedge angular coordinates for our new vertex
  for (Halfedge heIn : newV.incomingHalfedges()) {
    updateAngleFromCWNeighor(heIn);
  }

  // == (2) Set up bases on the intrinsic faces
  for (Face f : newV.adjacentFaces()) {
    updateFaceBasis(f);
  }

  // == (3) Find the insertion point on the input mesh, use it to align tangent spaces

  // Heuristic: we have to choose some edge to trace from to resolve the vertex. Use the shortest one, as it will often
  // involve the fewest crossings. Furthermore, use an original vertex if possible to reduce accumulating numerical
  // error.
  /*
  Halfedge inputTraceHe = newV.halfedge().twin();
  double shortestTraceHeLen = intrinsicEdgeLengths[inputTraceHe.edge()];
  for (Halfedge heIn : newV.incomingHalfedges()) {
    double thisLen = intrinsicEdgeLengths[heIn.edge()];

    bool currVertIsOriginal = vertexLocations[inputTraceHe.vertex()].type == SurfacePointType::Vertex;
    bool newVertIsOriginal = vertexLocations[heIn.vertex()].type == SurfacePointType::Vertex;

    if (currVertIsOriginal && !newVertIsOriginal) continue;

    if ((newVertIsOriginal && !currVertIsOriginal) || thisLen < shortestTraceHeLen) {
      shortestTraceHeLen = thisLen;
      inputTraceHe = heIn;
    }
  }
  */


  // We have to choose some edge to trace from to resolve the vertex. Choose from neighbors according to the following
  // priority:
  //  (best)    original points
  //  (neutral) other points
  //  (worst)   boundary points
  //  [break ties by shortest edge length]
  Halfedge inputTraceHe = newV.halfedge().twin();
  std::tuple<int, double> priorityBest{9999, 0};
  for (Halfedge heIn : newV.incomingHalfedges()) {

    // length score
    double thisLen = intrinsicEdgeLengths[heIn.edge()];

    // type score
    int numScore = 2;
    SurfacePoint candidateLoc = vertexLocations[inputTraceHe.vertex()];
    if (candidateLoc.type == SurfacePointType::Vertex) {
      numScore = 1;
    }
    if (heIn.edge().isBoundary()) {
      numScore = 3;
    }

    // combined score
    std::tuple<int, double> priorityThis{numScore, thisLen};

    // keep best
    if (priorityThis < priorityBest) {
      priorityBest = priorityThis;
      inputTraceHe = heIn;
    }
  }

  // Trace from the adjacent vertex selected above to get the position/angle on the input mesh
  SurfacePoint newPositionOnInput;
  Vector2 outgoingVec;

  bool boundaryEdgeInsertion = (intrinsicPoint.type == SurfacePointType::Edge && intrinsicPoint.edge.isBoundary());
  if (boundaryEdgeInsertion) {
    // For boundary vertices, instead of tracing, directly interpolate along the edge
    inputTraceHe = newV.halfedge().twin();

    // Get t values along the edge for previous and next vertices (note that we must be along a single edge)
    double tPrev, tNext;
    Edge origEdge;

    Vertex nextV = newV.halfedge().twin().vertex();
    if (vertexLocations[nextV].type == SurfacePointType::Vertex) {
      tNext = 1.;
    } else /* SurfacePointType::Edge */ {
      tNext = vertexLocations[nextV].tEdge;
      origEdge = vertexLocations[nextV].edge;
    }

    Vertex prevV = newV.halfedge().twin().next().twin().vertex();
    if (vertexLocations[prevV].type == SurfacePointType::Vertex) {
      tPrev = 0.;
    } else /* SurfacePointType::Edge */ {
      tPrev = vertexLocations[prevV].tEdge;
      origEdge = vertexLocations[prevV].edge;
    }

    // If neither was along an edge, find the (should-be-unique) boundary edge connecting the vertices
    if (origEdge == Edge()) {
      Vertex nextVOrig = vertexLocations[nextV].vertex;
      Vertex prevVOrig = vertexLocations[prevV].vertex;

      for (Halfedge he : nextVOrig.incomingHalfedges()) {
        if (he.vertex() == prevVOrig && he.edge().isBoundary()) {
          origEdge = he.edge();
        }
      }
    }

    if (origEdge == Edge()) {
      throw std::runtime_error("edge split problem: no boundary edge connecting vertices in boundary insertion");
    }

    // Interpolate between the previous and next t values
    double localT = intrinsicPoint.tEdge;
    double thisT = (1. - localT) * tPrev + (localT)*tNext;

    newPositionOnInput = SurfacePoint(origEdge, thisT);
  } else {
    // Normal case: trace an edge inward, use the result to resolve position and tangent basis
    // std::cout << "tracing to resolve new vertex" << std::endl;
    TraceOptions options;

    TraceGeodesicResult inputTraceResult =
        traceGeodesic(inputGeom, vertexLocations[inputTraceHe.vertex()], halfedgeVector(inputTraceHe), options);
    // std::cout << " --> done tracing to resolve new vertex" << std::endl;
    // snapEndToEdgeIfClose(inputTraceResult); TODO
    newPositionOnInput = inputTraceResult.endPoint;
    outgoingVec = -inputTraceResult.endingDir;
  }


  // Set the location of our newly inserted vertex
  vertexLocations[newV] = newPositionOnInput;

  // Align the new vertex's tangent space to that of the input mesh.
  double incomingAngle = standardizeAngle(newV, outgoingVec.arg());
  if (!inputTraceHe.isInterior()) {
    incomingAngle = 0;
  }

  intrinsicHalfedgeDirections[inputTraceHe.twin()] = incomingAngle;
  // halfedgeVectorsInVertex[inputTraceHe.twin()] = outgoingVec.normalize() * intrinsicEdgeLengths[inputTraceHe.edge()];
  halfedgeVectorsInVertex[inputTraceHe.twin()] = halfedgeVector(inputTraceHe.twin());

  // Custom loop to orbit CCW from InputTraceHe
  Halfedge firstHe = inputTraceHe.twin();
  Halfedge currHe = firstHe.next().next().twin();
  do {
    updateAngleFromCWNeighor(currHe);
    if (!currHe.isInterior()) break;
    currHe = currHe.next().next().twin();
  } while (currHe != firstHe);
}


void SignpostIntrinsicTriangulation::invokeEdgeFlipCallbacks(Edge e) {
  for (auto& fn : edgeFlipCallbackList) {
    fn(e);
  }
}
void SignpostIntrinsicTriangulation::invokeFaceInsertionCallbacks(Face f, Vertex v) {
  for (auto& fn : faceInsertionCallbackList) {
    fn(f, v);
  }
}
void SignpostIntrinsicTriangulation::invokeEdgeSplitCallbacks(Edge e, Halfedge he1, Halfedge he2) {
  for (auto& fn : edgeSplitCallbackList) {
    fn(e, he1, he2);
  }
}

} // namespace surface
} // namespace geometrycentral
