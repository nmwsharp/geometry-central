#include "geometrycentral/surface/signpost_intrinsic_triangulation.h"

#include "geometrycentral/surface/trace_geodesic.h"

#include <queue>

using std::cout;
using std::endl;

namespace geometrycentral {
namespace surface {


SignpostIntrinsicTriangulation::SignpostIntrinsicTriangulation(IntrinsicGeometryInterface& inputGeom_)
    // Note: this initializer list does something slightly wacky: it creates the new mesh on the heap, then loses track
    // of pointer while setting the BaseGeometryInterface::mesh reference to it. Later, it picks the pointer back up and
    // wraps it in the intrinsicMesh unique_ptr<>. I believe that this is all valid, but its probably a sign of bad
    // design.
    : IntrinsicGeometryInterface(*inputGeom_.mesh.copy().release()), inputMesh(inputGeom_.mesh), inputGeom(inputGeom_),
      intrinsicMesh(&mesh) {

  // == Initialize geometric data
  inputGeom.requireEdgeLengths();
  inputGeom.requireHalfedgeVectorsInVertex();
  inputGeom.requireVertexAngleSums();

  // Just copy lengths and angle sums
  intrinsicEdgeLengths = inputGeom.edgeLengths.reinterpretTo(mesh);
  intrinsicVertexAngleSums = inputGeom.vertexAngleSums.reinterpretTo(mesh);

  // Convert directions to radians
  intrinsicHalfedgeDirections = HalfedgeData<double>(mesh);
  HalfedgeData<Vector2> halfedgeVecsOnIntrinsic = inputGeom.halfedgeVectorsInVertex.reinterpretTo(mesh);
  for (Halfedge he : mesh.halfedges()) {
    Vertex baseVert = he.vertex();
    double scaledAngle = halfedgeVecsOnIntrinsic[he].arg();
    double rawAngle = scaledAngle * vertexAngleScaling(baseVert);
    double rawAngleStand = standardizeAngle(baseVert, rawAngle);
    if (!he.isInterior()) {
      // mod can do bad things to last angle, make sure its right
      rawAngleStand = intrinsicVertexAngleSums[baseVert];
    }
    intrinsicHalfedgeDirections[he] = rawAngleStand;
  }

  // Initialize vertex locations
  vertexLocations = VertexData<SurfacePoint>(mesh);
  for (size_t iV = 0; iV < mesh.nVertices(); iV++) {
    vertexLocations[iV] = SurfacePoint(inputMesh.vertex(iV));
  }

  requireHalfedgeVectorsInVertex();
  requireHalfedgeVectorsInFace();
  requireVertexAngleSums();
}


EdgeData<std::vector<SurfacePoint>> SignpostIntrinsicTriangulation::traceEdges() {

  EdgeData<std::vector<SurfacePoint>> tracedEdges(mesh);

  for (Edge e : mesh.edges()) {
    Halfedge he = e.halfedge();

    // Gather values to trace
    SurfacePoint startP = vertexLocations[he.vertex()];
    Vector2 traceVec = halfedgeVector(he);

    // Do the actual tracing
    TraceGeodesicResult result = traceGeodesic(inputGeom, startP, traceVec, true);

    // Save result
    tracedEdges[e] = result.pathPoints;
  }

  return tracedEdges;
}


// ======================================================
// ======== Queries & Accessors
// ======================================================


bool SignpostIntrinsicTriangulation::isDelaunay() {
  for (Edge e : mesh.edges()) {
    if (!e.isBoundary() && edgeCotanWeight(e) < -delaunayEPS) {
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

Vertex SignpostIntrinsicTriangulation::insertVertex_edge(SurfacePoint newPositionOnIntrinsic) {
  throw std::logic_error("not yet implemented");
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

  // === Phase 1: Find the geodesic ray which takes us to the circumcenter
  Vertex circumcenterRaySourceVert;
  Vector2 circumcenterRayVector; // in vertex coords
  {

    // = First, find the vertex with the widest angle in the triangle.
    // Working based from that vertex will be useful, because then the vector towards the circumcenter is contained in
    // the wedge formed by that angle, which allows us to think less about tangent space calculations.
    double biggestAngle = -1.;
    Halfedge he0 = f.halfedge();
    for (Halfedge he : f.adjacentHalfedges()) {
      double baseAngle = cornerAngle(he.corner());
      if (baseAngle > biggestAngle) {
        biggestAngle = baseAngle;
        he0 = he;
      }
    }

    // Gather some values
    // (wedge vertex is at (0,0)
    Vector2 p1 = Vector2{intrinsicEdgeLengths[he0.edge()], 0.0};
    Vector2 p2 = intrinsicEdgeLengths[he0.next().next().edge()] * Vector2::fromAngle(biggestAngle);
    Vector2 midPoint1 = 0.5 * p1;
    Vector2 midPoint2 = 0.5 * p2;

    // Circumcenter is intersection of perpendicular bisectors
    Vector2 vec1{0.0, 1.0};
    Vector2 vec2 = -p2.rotate90().normalize();
    Vector2 diffV = midPoint2 - midPoint1;
    double t = cross(vec2, diffV) / cross(vec2, vec1);
    Vector2 circumcenterVecLocal = midPoint1 + t * vec1;

    // Vector which takes us to the circumcenter, in coordinates of vertex
    double circumcenterRayLength = circumcenterVecLocal.norm();
    circumcenterRaySourceVert = he0.vertex();
    circumcenterRayVector =
        halfedgeVector(he0) * circumcenterVecLocal.pow(vertexAngleScaling(circumcenterRaySourceVert));
    circumcenterRayVector = circumcenterRayVector.normalize() * circumcenterRayLength;
  }


  // === Phase 2: Trace the ray to find the location of the new point on the intrinsic meshes

  // Data we need from the intrinsic trace
  SurfacePoint newPositionOnIntrinsic;
  double geodesicAngleOnIntrinsic;

  { // Do the actual tracing

    TraceGeodesicResult intrinsicTraceResult =
        traceGeodesic(*this, SurfacePoint(circumcenterRaySourceVert), circumcenterRayVector, false);

    // intrinsicTracer->snapEndToEdgeIfClose(intrinsicCrumbs); TODO
    newPositionOnIntrinsic = intrinsicTraceResult.endPoint.inSomeFace();
  }

  // === Phase 3: Add the new vertex
  return insertVertex(newPositionOnIntrinsic);
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
}

void SignpostIntrinsicTriangulation::delaunyRefine(double angleThreshDegrees, double circumradiusThresh,
                                                   size_t maxInsertions) {

  // Relationship between angles and circumradius-to-edge
  double angleThreshRad = angleThreshDegrees * M_PI / 180.;
  double circumradiusEdgeRatioThresh = 1.0 / (2.0 * std::sin(angleThreshRad));


  // Helper to test if a face violates the circumradius ratio condition
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
    if (needsCircumcenterRefinement(f)) {
      circumradiusCheckQueue.push(std::make_pair(area(f), f));
    }
  }

  // === Outer iteration: flip and insert until we have a mesh that satisfies both angle and circumradius goals
  size_t nFlips = 0;
  size_t nInsertions = 0;
  do {

    //cout << "nFlips = " << nFlips << "nInsertions = " << nInsertions << endl;

    // == First, flip to delaunay
    while (!delaunayCheckQueue.empty()) {

      //cout << "nFlips = " << nFlips << "nInsertions = " << nInsertions << endl;

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
        if (needsCircumcenterRefinement(nF)) {
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
    // TODO handle boundary

    // If we've already inserted the max number of points, empty the queue and call it a day
    if (maxInsertions != INVALID_IND && nInsertions == maxInsertions) {
      circumradiusCheckQueue = std::priority_queue<AreaFace, std::vector<AreaFace>, std::less<AreaFace>>();
    }

    // Try to insert just one circumcenter
    if (!circumradiusCheckQueue.empty()) {

      // Get the biggest face
      Face f = circumradiusCheckQueue.top().second;
      double A = circumradiusCheckQueue.top().first;
      circumradiusCheckQueue.pop();

      // If the area has changed since this face was inserted in to the queue, skip it. Note that we don't need to
      // re-add it, because it must have been placed in the queue when its area was changed
      if (A != area(f)) {
        continue;
      }

      // This face might have been flipped to no longer violate constraint
      if (needsCircumcenterRefinement(f)) {

        Vertex newVert = insertCircumcenter(f);
        nInsertions++;

        if (nInsertions % 1000 == 0) {
          cout << "   ...inserted " << nInsertions << " points" << endl;
        }

        // Mark everything in the 1-ring as possibly non-Delaunay and possibly violating the circumradius constraint
        for (Face nF : newVert.adjacentFaces()) {

          // Check circumradius constraint
          if (needsCircumcenterRefinement(nF)) {
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

    // If the circumradius queue is empty, make sure we didn't miss anything (can happen due to numerics)
    /*
    for (Face f : mesh.faces()) {
      if (needsCircumcenterRefinement(f)) {
        circumradiusCheckQueue.push(std::make_pair(area(f), f));
      }
    }
    */

  } while (!delaunayCheckQueue.empty() || !circumradiusCheckQueue.empty());

  // TODO
  // Check all faces
  for (Face f : mesh.faces()) {
    if (needsCircumcenterRefinement(f)) {
      cout << "missed face!" << endl;
    }
  }


  cout << "  ... intrinsic Delaunay refinment finished. Took " << nFlips << " flips and " << nInsertions
       << " insertions." << endl;
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

  // Trace from the adjacent vertex selected above to get the position/angle on the input mesh
  TraceGeodesicResult inputTraceResult =
      traceGeodesic(inputGeom, vertexLocations[inputTraceHe.vertex()], halfedgeVector(inputTraceHe));


  // Set the location of our newly inserted vertex
  SurfacePoint newPositionOnInput = inputTraceResult.endPoint;
  vertexLocations[newV] = newPositionOnInput;

  // Align the new vertex's tangent space to that of the input mesh.
  Vector2 outgoingVec = -inputTraceResult.endingDir;
  double incomingAngle = standardizeAngle(newV, outgoingVec.arg());
  intrinsicHalfedgeDirections[inputTraceHe.twin()] = incomingAngle;
  halfedgeVectorsInVertex[inputTraceHe.twin()] = outgoingVec.normalize() * intrinsicEdgeLengths[inputTraceHe.edge()];

  // Custom loop to orbit CCW from InputTraceHe
  Halfedge firstHe = inputTraceHe.twin();
  Halfedge currHe = firstHe.next().next().twin();
  do {
    updateAngleFromCWNeighor(currHe);
    currHe = currHe.next().next().twin();
  } while (currHe != firstHe);
}


} // namespace surface
} // namespace geometrycentral
