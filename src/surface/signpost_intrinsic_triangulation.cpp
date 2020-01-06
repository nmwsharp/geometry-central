#include "geometrycentral/surface/signpost_intrinsic_triangulation.h"

#include "geometrycentral/surface/barycentric_coordinate_helpers.h"
#include "geometrycentral/surface/trace_geodesic.h"

#include <iomanip>
#include <queue>

using std::cout;
using std::endl;

namespace geometrycentral {
namespace surface {


SignpostIntrinsicTriangulation::SignpostIntrinsicTriangulation(IntrinsicGeometryInterface& inputGeom_)
    // Note: this initializer list does something slightly wacky: it creates the new mesh on the heap, then loses track
    // of pointer while setting the BaseGeometryInterface::mesh reference to it. Later, it picks the pointer back up
    // from the reference and wraps it in the intrinsicMesh unique_ptr<>. I believe that this is all valid, but its
    // probably a sign of bad design.
    : IntrinsicGeometryInterface(*inputGeom_.mesh.copy().release()), inputMesh(inputGeom_.mesh), inputGeom(inputGeom_),
      intrinsicMesh(&mesh) {

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
  edgeIsOriginal = EdgeData<char>(mesh, false);
  edgeIsOriginal.fill(true);

  requireHalfedgeVectorsInVertex();
  requireHalfedgeVectorsInFace();
  requireVertexAngleSums();
}

std::vector<SurfacePoint> SignpostIntrinsicTriangulation::traceHalfedge(Halfedge he) {
  // Gather values to trace
  SurfacePoint startP = vertexLocations[he.vertex()];
  Vector2 traceVec = halfedgeVector(he);

  // Do the actual tracing
  TraceGeodesicResult result = traceGeodesic(inputGeom, startP, traceVec, true);

  // Trim off end crumbs if applicable
  Vertex endVert = he.twin().vertex();
  if (vertexLocations[endVert].type == SurfacePointType::Vertex) {
    bool success = trimTraceResult(result, endVert);

    // If trimming failed (because the trace didn't even hit the 1-ring of target), just stick with whatever we go
    // initially
    if (!success) {
      result = traceGeodesic(inputGeom, startP, traceVec, true);
    }
  }

  // Append the endpoint
  result.pathPoints.push_back(vertexLocations[endVert]);

  return result.pathPoints;
}

EdgeData<std::vector<SurfacePoint>> SignpostIntrinsicTriangulation::traceEdges() {

  EdgeData<std::vector<SurfacePoint>> tracedEdges(mesh);

  for (Edge e : mesh.edges()) {
    Halfedge he = e.halfedge();
    tracedEdges[e] = traceHalfedge(he);
  }

  return tracedEdges;
}


// ======================================================
// ======== Queries & Accessors
// ======================================================


bool SignpostIntrinsicTriangulation::isDelaunay(Edge e) {
  if (!e.isBoundary() && edgeCotanWeight(e) < -delaunayEPS) {
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
  if (e.isBoundary()) return false;

  // Don't want to flip
  double cWeight = edgeCotanWeight(e);
  if (cWeight > -delaunayEPS) return false;

  // Get geometric data
  Halfedge he = e.halfedge();
  std::array<Vector2, 4> layoutPositions = layoutDiamond(he);

  // Combinatorial flip
  bool flipped = mesh.flip(e);

  // Should always be possible, something unusual is going on if we end up here
  if (!flipped) {
    return false;
  }

  // Compute the new edge length
  double newLength = (layoutPositions[1] - layoutPositions[3]).norm();

  // If we're going to create a non-finite edge length, abort the flip
  // (only happens if you're in a bad numerical place)
  if (!std::isfinite(newLength)) {
    mesh.flip(e);
    return false;
  }

  // Assign the new edge lengths
  intrinsicEdgeLengths[e] = newLength;

  // Update edge angles
  updateAngleFromCWNeighor(e.halfedge());
  updateAngleFromCWNeighor(e.halfedge().twin());
  updateFaceBasis(e.halfedge().face());
  updateFaceBasis(e.halfedge().twin().face());

  edgeIsOriginal[e] = false;

  return true;
}

bool SignpostIntrinsicTriangulation::flipEdgeIfPossible(Edge e, double possibleEPS) {

  // Can't flip
  if (e.isBoundary()) return false;

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
  bool flipped = mesh.flip(e);

  // Might not have been flippable for connectivity reasons
  if (!flipped) {
    return false;
  }

  // Compute the new edge length
  double newLength = (layoutPositions[1] - layoutPositions[3]).norm();

  // If we're going to create a non-finite edge length, abort the flip
  // (only happens if you're in a bad numerical place)
  if (!std::isfinite(newLength)) {
    mesh.flip(e);
    return false;
  }

  // Assign the new edge lengths
  // TODO project to satisfy triangle inequality?
  intrinsicEdgeLengths[e] = newLength;

  // Update edge angles
  updateAngleFromCWNeighor(e.halfedge());
  updateAngleFromCWNeighor(e.halfedge().twin());
  updateFaceBasis(e.halfedge().face());
  updateFaceBasis(e.halfedge().twin().face());

  edgeIsOriginal[e] = false;

  return true;
}

Vertex SignpostIntrinsicTriangulation::insertVertex(SurfacePoint newPositionOnIntrinsic) {
  switch (newPositionOnIntrinsic.type) {
  case SurfacePointType::Vertex: {
    throw std::logic_error("can't insert vertex at vertex");
    break;
  }
  case SurfacePointType::Edge: {
    return insertVertex_edge(newPositionOnIntrinsic);
    break;
  }
  case SurfacePointType::Face: {
    return insertVertex_face(newPositionOnIntrinsic);
    break;
  }
  }
  return Vertex();
}

Vertex SignpostIntrinsicTriangulation::insertVertex_edge(SurfacePoint newP) {

  std::cout << "inserting at edge " << newP << std::endl;

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
  Vector2 posA = (1. * newP.tEdge) * vertCoords[iA] + newP.tEdge * vertCoords[(iA + 1) % 3];
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
  Halfedge newHeFront = mesh.splitEdgeTriangular(insertionEdge);
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
  std::array<double, 4> newLens = {frontLen, Alen, backLen, Blen};
  for (int i = 0; i < (isOnBoundary ? 3 : 4); i++) {
    std::cout << "setting new edge len " << currHe.edge() << " = " << newLens[i] << std::endl;
    intrinsicEdgeLengths[currHe.edge()] = newLens[i];
    currHe = currHe.next().next().twin();
  }

  // === (4) Now that we have edge lengths, sort out tangent spaces and position on supporting.
  resolveNewVertex(newV);

  return newV;
}

Vertex SignpostIntrinsicTriangulation::insertVertex_face(SurfacePoint newP) {

  std::cout << "inserting at face " << newP << std::endl;

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
  Vertex newV = mesh.insertVertex(insertionFace);

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
      }
    }
  }

  // === (4) Now that we have edge lengths, sort out tangent spaces and position on supporting.
  resolveNewVertex(newV);

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
  TraceGeodesicResult intrinsicTraceResult = traceGeodesic(*this, f, barycenter, vecToCircumcenter, false);
  // intrinsicTracer->snapEndToEdgeIfClose(intrinsicCrumbs); TODO
  // SurfacePoint newPositionOnIntrinsic = intrinsicTraceResult.endPoint.inSomeFace();
  SurfacePoint newPositionOnIntrinsic = intrinsicTraceResult.endPoint;


  // === Phase 3: Add the new vertex
  return insertVertex(newPositionOnIntrinsic);
}

Vertex SignpostIntrinsicTriangulation::insertBarycenter(Face f) {
  SurfacePoint barycenterOnIntrinsic(f, Vector3::constant(1. / 3.));
  return insertVertex(barycenterOnIntrinsic);
}

void SignpostIntrinsicTriangulation::removeInsertedVertex(Vertex v) {
  // Strategy: flip edges until the vertex has degree three, then remove by replacing with a single face
  // TODO needs a proof that this always works... what about self edges, etc?

  // What about starting with degree < 3? Since this vertex necessarily has angle sum 2PI, this could only happen in the
  // case of degree 2, with exactly degenerate triangles. Since we assume non-degenerate triangles throughout, we'll
  // consider that to not happen.

  // Flip edges until
  size_t iterCount = 0;
  while (v.degree() != 3)  {

    // Find any edge we can flip
    bool anyFlipped = false;
    for(Edge e : v.adjacentEdges()) {
      anyFlipped = flipEdgeIfNotDelaunay(e);
      if(anyFlipped) break;
    }

    // failsafe, in case we get numerically stuck (or the algorithm is broken)
    if(!anyFlipped || iterCount > 10 * v.degree()) return;
  }

  // give up if something went wrong
  if(v.degree() != 3) return;

  //  
}

void SignpostIntrinsicTriangulation::flipToDelaunay() {

  std::deque<Edge> edgesToCheck;
  EdgeData<char> inQueue(mesh, true);
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


  // Initialize queue of (possibly) non-delaunay edges
  std::deque<Edge> delaunayCheckQueue;
  EdgeData<char> inDelaunayQueue(mesh, false);
  for (Edge e : mesh.edges()) {
    delaunayCheckQueue.push_back(e);
    inDelaunayQueue[e] = true;
  }


  // Initialize queue of (possibly) circumradius-violating faces, processing the largest faces first (good heuristic)
  typedef std::pair<double, Face> AreaFace;
  std::priority_queue<AreaFace, std::vector<AreaFace>, std::less<AreaFace>> circumradiusCheckQueue;
  for (Face f : mesh.faces()) {
    if (shouldRefine(f)) {
      circumradiusCheckQueue.push(std::make_pair(area(f), f));
    }
  }


  // === Outer iteration: flip and insert until we have a mesh that satisfies both angle and circumradius goals
  size_t nFlips = 0;
  size_t nInsertions = 0;
  do {

    // == First, flip to delaunay
    while (!delaunayCheckQueue.empty()) {

      // Get the top element from the queue of possibily non-Delaunay edges
      Edge e = delaunayCheckQueue.front();
      delaunayCheckQueue.pop_front();
      inDelaunayQueue[e] = false;

      bool wasFlipped = flipEdgeIfNotDelaunay(e);

      if (!wasFlipped) continue;

      // Handle the aftermath of a flip
      nFlips++;

      // Add neighboring faces, which might violate circumradius constraint
      std::vector<Face> neighFaces = {e.halfedge().face(), e.halfedge().twin().face()};
      for (Face nF : neighFaces) {
        if (shouldRefine(nF)) {
          circumradiusCheckQueue.push(std::make_pair(area(nF), nF));
        }
      }

      // Note: area of faces currently in queue is potentially changing due to flips, and we don't update the area.
      // That's fine; it's just a heuristic.

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
    }

    // == Second, insert one circumcenter

    // If we've already inserted the max number of points, empty the queue and call it a day
    if (maxInsertions != INVALID_IND && nInsertions == maxInsertions) {
      // circumradiusCheckQueue = std::priority_queue<AreaFace, std::vector<AreaFace>, std::less<AreaFace>>();
      break;
    }

    // Try to insert just one circumcenter
    if (!circumradiusCheckQueue.empty()) {

      // Get the biggest face
      Face f = circumradiusCheckQueue.top().second;
      double A = circumradiusCheckQueue.top().first;
      circumradiusCheckQueue.pop();

      // Two things might have changed that would cause us to skip this entry:
      //   -If the area has changed since this face was inserted in to the queue, skip it. Note that we don't need to
      //    re-add it, because it must have been placed in the queue when its area was changed
      //   - This face might have been flipped to no longer violate constraint
      if (A == area(f) && shouldRefine(f)) {

        Vertex newVert = insertCircumcenter(f);
        nInsertions++;

        // Mark everything in the 1-ring as possibly non-Delaunay and possibly violating the circumradius constraint
        for (Face nF : newVert.adjacentFaces()) {

          // Check circumradius constraint
          if (shouldRefine(nF)) {
            circumradiusCheckQueue.push(std::make_pair(area(nF), nF));
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
    }

    // If the circumradius queue is empty, make sure we didn't miss anything (can happen rarely due to numerics)
    // (but don't do this more than a few times, to avoid getting stuck in an infinite loop when numerical ultra-badness
    // happens)
    if (recheckCount < MAX_RECHECK_COUNT) {
      recheckCount++;
      if (delaunayCheckQueue.empty() && circumradiusCheckQueue.empty()) {
        for (Face f : mesh.faces()) {
          if (shouldRefine(f)) {
            circumradiusCheckQueue.push(std::make_pair(area(f), f));
          }
        }
        for (Edge e : mesh.edges()) {
          if (!isDelaunay(e)) {
            delaunayCheckQueue.push_back(e);
            inDelaunayQueue[e] = true;
          }
        }
      }
    }

  } while (!delaunayCheckQueue.empty() || !circumradiusCheckQueue.empty());

  refreshQuantities();
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


void SignpostIntrinsicTriangulation::resolveNewVertex(Vertex newV) {

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
  //  (best)  edge points along boundary
  //          original points
  //  (worst) other points
  //  [break ties by shortest edge length]
  Halfedge inputTraceHe = newV.halfedge().twin();
  std::tuple<int, double> priorityBest{9999, 0};
  for (Halfedge heIn : newV.incomingHalfedges()) {

    // length score
    double thisLen = intrinsicEdgeLengths[heIn.edge()];

    // type score
    int numScore = 3;
    SurfacePoint candidateLoc = vertexLocations[inputTraceHe.vertex()];
    if (candidateLoc.type == SurfacePointType::Vertex) {
      numScore = 2;
    }
    if (!heIn.isInterior()) { // NOTE: this is important! makes CCW loop below correct for boundary
      numScore = 1;
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
  TraceGeodesicResult inputTraceResult =
      traceGeodesic(inputGeom, vertexLocations[inputTraceHe.vertex()], halfedgeVector(inputTraceHe));


  // Set the location of our newly inserted vertex
  SurfacePoint newPositionOnInput = inputTraceResult.endPoint;
  vertexLocations[newV] = newPositionOnInput;

  // Align the new vertex's tangent space to that of the input mesh.
  Vector2 outgoingVec = -inputTraceResult.endingDir;
  double incomingAngle = standardizeAngle(newV, outgoingVec.arg());
  if (!inputTraceHe.isInterior()) {
    incomingAngle = 0;
  }

  intrinsicHalfedgeDirections[inputTraceHe.twin()] = incomingAngle;
  std::cout << "incoming angle = " << incomingAngle << std::endl;
  // halfedgeVectorsInVertex[inputTraceHe.twin()] = outgoingVec.normalize() * intrinsicEdgeLengths[inputTraceHe.edge()];
  halfedgeVectorsInVertex[inputTraceHe.twin()] = halfedgeVector(inputTraceHe.twin());

  // Custom loop to orbit CCW from InputTraceHe
  // NOTE: this is only correct because we always choose incoming boundary edge to trace above
  Halfedge firstHe = inputTraceHe.twin();
  Halfedge currHe = firstHe.next().next().twin();
  do {
    updateAngleFromCWNeighor(currHe);
    currHe = currHe.next().next().twin();
  } while (currHe != firstHe);
}

} // namespace surface
} // namespace geometrycentral
