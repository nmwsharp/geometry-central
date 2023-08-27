#include "geometrycentral/surface/signpost_intrinsic_triangulation.h"

#include "geometrycentral/surface/barycentric_coordinate_helpers.h"
#include "geometrycentral/surface/mesh_graph_algorithms.h"
#include "geometrycentral/surface/trace_geodesic.h"

#include <iomanip>
#include <queue>

namespace geometrycentral {
namespace surface {

SignpostIntrinsicTriangulation::SignpostIntrinsicTriangulation(ManifoldSurfaceMesh& mesh_,
                                                               IntrinsicGeometryInterface& inputGeom_)
    : IntrinsicTriangulation(mesh_, inputGeom_) {

  // == Initialize geometric data
  // TODO do we really need all these?
  inputGeom.requireEdgeLengths();
  inputGeom.requireHalfedgeVectorsInVertex();
  inputGeom.requireHalfedgeVectorsInFace();
  inputGeom.requireVertexAngleSums();

  // Prepare directions and angle sums
  signpostAngle = HalfedgeData<double>(mesh);

  // Walk around the vertex, constructing angular directions
  for (Vertex v : mesh.vertices()) {
    double runningAngle = 0.;
    Halfedge firstHe = v.halfedge();
    Halfedge currHe = firstHe;
    do {
      signpostAngle[currHe] = runningAngle;

      if (!currHe.isInterior()) {
        break;
      }

      double cornerAngleVal = cornerAngle(currHe.corner());
      runningAngle += cornerAngleVal;

      currHe = currHe.next().next().twin();
    } while (currHe != firstHe);
  }

  // Initialize all edges as original, but new ones should be false
  edgeIsOriginal = EdgeData<bool>(mesh, false);
  edgeIsOriginal.fill(true);
}

std::vector<SurfacePoint> SignpostIntrinsicTriangulation::traceIntrinsicHalfedgeAlongInput(Halfedge he) {
  return traceIntrinsicHalfedgeAlongInput(he, true);
}

std::vector<SurfacePoint> SignpostIntrinsicTriangulation::traceIntrinsicHalfedgeAlongInput(Halfedge he, bool trimEnd) {

  // Optimization: don't bother tracing original edges, just report them directly
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


std::vector<SurfacePoint> SignpostIntrinsicTriangulation::traceInputHalfedgeAlongIntrinsic(Halfedge inputHe) {
  return traceInputHalfedgeAlongIntrinsic(inputHe, true);
}

std::vector<SurfacePoint> SignpostIntrinsicTriangulation::traceInputHalfedgeAlongIntrinsic(Halfedge inputHe,
                                                                                           bool trimEnd) {

  // Optimization: don't bother tracing original edges, just report them directly
  // Here we need to check if the corresponding intrinsic edge is original
  Halfedge origHe = intrinsicMesh->halfedge(inputHe.getIndex()); // as usual, safe by dense construction
  if (edgeIsOriginal[origHe.edge()]) {
    Vertex vA = inputHe.tailVertex();
    Vertex vB = inputHe.tipVertex();
    std::vector<SurfacePoint> result{SurfacePoint(vA), SurfacePoint(vB)};
    return result;
  }


  Vertex vTail = inputHe.tailVertex();
  Vertex vTip = inputHe.tipVertex();

  // Get the corresponding tail vertex on the intrinsic triangulation
  Vertex vTailIntrinsic = intrinsicMesh->vertex(vTail.getIndex()); // as usual, safe by dense construction

  // Gather values to trace
  SurfacePoint startP = SurfacePoint(vTailIntrinsic);
  Vector2 traceVec = inputGeom.halfedgeVectorsInVertex[inputHe];

  // Do the actual tracing
  TraceOptions options;
  options.includePath = true;
  options.maxIters = inputMesh.nFaces() * 10;
  TraceGeodesicResult result = traceGeodesic(*this, startP, traceVec, options);

  // Trim off end crumbs if applicable
  Vertex vTipIntrinsic = intrinsicMesh->vertex(vTip.getIndex()); // as usual, safe by dense construction
  if (trimEnd) {
    bool success = trimTraceResult(result, vTipIntrinsic);
    if (success) {
      // Append the endpoint
      result.pathPoints.push_back(SurfacePoint(vTipIntrinsic));
    } else {
      // If trimming failed (because the trace didn't even hit the 1-ring of target), just stick with whatever we go
      // initially
      result = traceGeodesic(*this, startP, traceVec, options);
    }
  }

  return result.pathPoints;
}

SurfacePoint SignpostIntrinsicTriangulation::equivalentPointOnIntrinsic(const SurfacePoint& pointOnInput) {

  // If it's a vertex, just return the matching vertex
  if (pointOnInput.type == SurfacePointType::Vertex) {
    // Get the corresponding point on the intrinsic triangulation. The getIndex() is safe in this case: these vertices
    // will always be densely enumerated by construction.
    Vertex intrinsicVertex = intrinsicMesh->vertex(pointOnInput.vertex.getIndex());
    return SurfacePoint(intrinsicVertex);
  }


  // TODO: could be smarter here and explicitly return edge points for shared edges.
  // Currently the code below will generally return some extremely-nearby point inside of some face.


  // Handle the general case of a point inside a face or a not-necessarily-shared edge by tracing along the surface from
  // some adjacent vertex.
  SurfacePoint inputFacePoint = pointOnInput.inSomeFace();
  Vector3 bary = inputFacePoint.faceCoords;

  // Pick a vertex to trace from
  int traceI = 0;
  Halfedge traceHe;
  double smallestBary = 999;
  Halfedge testHe = inputFacePoint.face.halfedge();
  for (size_t iV = 0; iV < 3; iV++) {
    // Don't trace from a vertex if we are nearly on top of it, because it will mess up
    // the angle calculations below. Smallest bary coord is a good choice.
    if (bary[iV] < smallestBary) {
      smallestBary = bary[iV];
      traceI = iV;
      traceHe = testHe;
    }
    testHe = testHe.next();
  }
  Vertex traceVertex = traceHe.vertex();
  SurfacePoint traceVertexLoc =
      SurfacePoint(intrinsicMesh->vertex(traceVertex.getIndex())); // (see note abouve about getIndex())

  // Compute the direction and length of the vector pointing towards the vertex
  std::array<Vector2, 3> coords = {Vector2{0., 0.}, inputGeom.halfedgeVectorsInFace[inputFacePoint.face.halfedge()],
                                   -inputGeom.halfedgeVectorsInFace[inputFacePoint.face.halfedge().next().next()]};
  Vector2 vertPos = coords[traceI];
  Vector2 pointPos = coords[0] * bary[0] + coords[1] * bary[1] + coords[2] * bary[2];
  Vector2 vec = pointPos - vertPos;
  double relativeAngle = angle(vec, inputGeom.halfedgeVectorsInFace[traceHe]);
  double relativeScaledAngle;
  if (traceVertex.isBoundary()) {
    // follow convention that boundary verts are normalized to PI rather than 2PI
    relativeScaledAngle = relativeAngle * (M_PI / inputGeom.vertexAngleSums[traceVertex]);
  } else {
    relativeScaledAngle = relativeAngle * (2 * M_PI / inputGeom.vertexAngleSums[traceVertex]);
  }
  Vector2 vertexDir = unit(inputGeom.halfedgeVectorsInVertex[traceHe]) * Vector2::fromAngle(relativeScaledAngle);
  double len = norm(vec);
  Vector2 traceVec = vertexDir * len;

  // Trace out the vector along the original surface
  TraceOptions options;
  TraceGeodesicResult intrinsicTraceResult = traceGeodesic(*this, traceVertexLoc, traceVec, options);
  SurfacePoint intrinsicPos = intrinsicTraceResult.endPoint;

  return intrinsicPos;
}

SurfacePoint SignpostIntrinsicTriangulation::equivalentPointOnInput(const SurfacePoint& pointOnIntrinsic) {

  // If it's a vertex, just return the vertex location
  if (pointOnIntrinsic.type == SurfacePointType::Vertex) {
    return vertexLocations[pointOnIntrinsic.vertex];
  }

  // TODO: could be smarter here and explicitly return edge points for shared edges.
  // Currently the code below will generally return some extremely-nearby point inside of some face.


  // Handle the general case of a point inside a face or a not-necessarily-shared edge by tracing along the surface from
  // some adjacent vertex.
  SurfacePoint intrinsicFacePoint = pointOnIntrinsic.inSomeFace();
  Vector3 bary = intrinsicFacePoint.faceCoords;

  // Pick a vertex to trace from
  int traceI = 0;
  Halfedge traceHe;
  double smallestBary = 999;
  Halfedge testHe = intrinsicFacePoint.face.halfedge();
  for (size_t iV = 0; iV < 3; iV++) {
    // Don't trace from a vertex if we are nearly on top of it, because it will mess up
    // the angle calculations below. Smallest bary coord is a good choice.
    if (bary[iV] < smallestBary) {
      smallestBary = bary[iV];
      traceI = iV;
      traceHe = testHe;
    }
    testHe = testHe.next();
  }
  Vertex traceVertex = traceHe.vertex();

  // Compute the direction and length of the vector pointing towards the vertex
  std::array<Vector2, 3> coords = vertexCoordinatesInTriangle(intrinsicFacePoint.face);
  Vector2 vertPos = coords[traceI];
  Vector2 pointPos = coords[0] * bary[0] + coords[1] * bary[1] + coords[2] * bary[2];
  Vector2 vec = pointPos - vertPos;
  double relativeAngle = angle(vec, halfedgeVectorsInFace[traceHe]);
  if (relativeAngle < 0) relativeAngle += 2 * M_PI;
  double vertexAngle = signpostAngle[traceHe] + relativeAngle;
  vertexAngle = standardizeAngle(traceVertex, vertexAngle);
  double len = norm(vec);
  Vector2 traceVec = rescaledVertexVector(traceVertex, vertexAngle, len);

  // Trace out the vector along the original surface
  TraceOptions options;
  TraceGeodesicResult inputTraceResult = traceGeodesic(inputGeom, vertexLocations[traceVertex], traceVec, options);
  SurfacePoint inputPos = inputTraceResult.endPoint;

  return inputPos;
}


// ======================================================
// ======== Queries & Accessors
// ======================================================

bool SignpostIntrinsicTriangulation::checkEdgeOriginal(Edge e) const { return edgeIsOriginal[e]; }

// ======================================================
// ======== Mutators
// ======================================================

bool SignpostIntrinsicTriangulation::flipEdgeIfNotDelaunay(Edge e) {

  // Can't flip
  if (isFixed(e)) return false;

  // Don't want to flip
  double cWeight = edgeCotanWeight(e);
  if (cWeight > -triangleTestEPS) return false;

  // Get geometric data
  Halfedge he = e.halfedge();
  std::array<Vector2, 4> layoutPositions = layoutDiamond(he);

  // Compute the new edge length
  double newLength = (layoutPositions[1] - layoutPositions[3]).norm();

  // If we're going to create a non-finite edge length, abort the flip
  // (only happens if you're in a bad numerical place)
  if (!std::isfinite(newLength)) {
    return false;
  }

  // Combinatorial flip
  bool flipped = intrinsicMesh->flip(e, false);

  // Should always be possible, something unusual is going on if we end up here
  if (!flipped) {
    return false;
  }

  // Assign the new edge lengths
  edgeLengths[e] = newLength;

  // Update edge angles
  updateAngleFromCWNeighor(e.halfedge());
  updateAngleFromCWNeighor(e.halfedge().twin());
  updateFaceBasis(e.halfedge().face());
  updateFaceBasis(e.halfedge().twin().face());

  edgeIsOriginal[e] = false;

  triangulationChanged();
  invokeEdgeFlipCallbacks(e);
  return true;
}

bool SignpostIntrinsicTriangulation::flipEdgeIfPossible(Edge e) {

  // Can't flip
  if (isFixed(e)) return false;

  // Get geometric data
  Halfedge he = e.halfedge();
  std::array<Vector2, 4> layoutPositions = layoutDiamond(he);

  // Test if geometryically flippable flippable (both signed areas of new triangles are positive)
  double A1 = cross(layoutPositions[1] - layoutPositions[0], layoutPositions[3] - layoutPositions[0]);
  double A2 = cross(layoutPositions[3] - layoutPositions[2], layoutPositions[1] - layoutPositions[2]);
  double areaEPS = triangleTestEPS * (A1 + A2);
  if (A1 < areaEPS || A2 < areaEPS) {
    return false;
  }

  // Compute the new edge length
  double newLength = (layoutPositions[1] - layoutPositions[3]).norm();

  // If we're going to create a non-finite edge length, abort the flip
  // (only happens if you're in a bad numerical place)
  if (!std::isfinite(newLength)) {
    return false;
  }

  // Combinatorial flip
  bool flipped = intrinsicMesh->flip(e, false);

  // Might not have been flippable for connectivity reasons
  if (!flipped) {
    return false;
  }

  // Assign the new edge lengths
  // TODO project to satisfy triangle inequality?
  edgeLengths[e] = newLength;

  // Update edge angles
  updateAngleFromCWNeighor(e.halfedge());
  updateAngleFromCWNeighor(e.halfedge().twin());
  updateFaceBasis(e.halfedge().face());
  updateFaceBasis(e.halfedge().twin().face());

  edgeIsOriginal[e] = false;

  triangulationChanged();
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

  edgeLengths[e] = newLength;

  // Update other derived geometric data
  signpostAngle[e.halfedge()] = forwardAngle;
  signpostAngle[e.halfedge().twin()] = reverseAngle;
  halfedgeVectorsInVertex[e.halfedge()] = halfedgeVector(e.halfedge());
  halfedgeVectorsInVertex[e.halfedge().twin()] = halfedgeVector(e.halfedge().twin());
  updateFaceBasis(e.halfedge().face());
  updateFaceBasis(e.halfedge().twin().face());

  edgeIsOriginal[e] = isOrig;
  triangulationChanged();
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
  backLen = newP.tEdge * edgeLengths[insertionEdge];
  frontLen = (1. - newP.tEdge) * edgeLengths[insertionEdge];

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
    vertexAngleSums[newV] = M_PI;
  } else {
    vertexAngleSums[newV] = 2. * M_PI;
  }


  // == (3) Assign edge lengths to the new edges
  Halfedge currHe = newHeFront;
  Halfedge newHeBack;
  std::array<double, 4> newLens = {frontLen, Alen, backLen, Blen};
  for (int i = 0; i < (isOnBoundary ? 3 : 4); i++) {
    edgeLengths[currHe.edge()] = newLens[i];
    if (i == 2) newHeBack = currHe;
    currHe = currHe.next().next().twin();
  }

  // === (4) Now that we have edge lengths, sort out tangent spaces and position on supporting.
  resolveNewVertex(newV, newP);

  triangulationChanged();
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
  vertexAngleSums[newV] = 2. * M_PI;


  // == (3) Assign edge lengths to the new edges

  // Set edge lengths first by looking for the proper new edge
  for (size_t j = 0; j < 3; j++) {
    double thisLen = newEdgeLengths[j];
    Halfedge origHe = oldFaceHalfedges[j];

    // Find the new edge which this length belongs to
    for (Halfedge heV : newV.outgoingHalfedges()) {
      if (heV.next() == origHe) {
        edgeLengths[heV.edge()] = thisLen;
      }
    }
  }

  // === (4) Now that we have edge lengths, sort out tangent spaces and position on supporting.
  resolveNewVertex(newV, newP);

  triangulationChanged();
  invokeFaceInsertionCallbacks(insertionFace, newV);
  return newV;
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
  triangulationChanged();
  return newF;
}

Halfedge SignpostIntrinsicTriangulation::splitEdge(Halfedge he, double tSplit) {
  return insertVertex_edge(SurfacePoint(he, tSplit));
}


// ======================================================
// ======== Geometry and Helpers
// ======================================================

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
    signpostAngle[he] = vertexAngleSums[he.vertex()]; // last angle in boundary wedge
    halfedgeVectorsInVertex[he] = halfedgeVector(he);
    return;
  }
  if (!he.twin().isInterior()) {
    signpostAngle[he] = 0.; // first angle in boundary wedge
    halfedgeVectorsInVertex[he] = halfedgeVector(he);
    return;
  }

  // Get neighbor angle
  Halfedge cwHe = he.twin().next();
  double neighAngle = signpostAngle[cwHe];

  // Compute corner angle in between
  double cAngle = cornerAngle(cwHe.corner());

  // Set the updated angle
  double updatedAngle = standardizeAngle(he.vertex(), neighAngle + cAngle);
  signpostAngle[he] = updatedAngle;
  halfedgeVectorsInVertex[he] = halfedgeVector(he);
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
    double thisLen = edgeLengths[heIn.edge()];

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

  signpostAngle[inputTraceHe.twin()] = incomingAngle;
  // halfedgeVectorsInVertex[inputTraceHe.twin()] = outgoingVec.normalize() * edgeLengths[inputTraceHe.edge()];
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

void SignpostIntrinsicTriangulation::constructCommonSubdivision() {

  intrinsicMesh->compress();

  // Do all the tracing
  EdgeData<std::vector<SurfacePoint>> traces = traceAllIntrinsicEdgesAlongInput();

  // Create the new common subdivision object
  commonSubdivision.reset(new CommonSubdivision(inputMesh, *intrinsicMesh));
  CommonSubdivision& cs = *commonSubdivision;

  // Construct CommonSubdivisionPoints corresponding to shared vertices
  VertexData<CommonSubdivisionPoint*> aVtx(inputMesh);
  VertexData<CommonSubdivisionPoint*> bVtx(*intrinsicMesh);

  for (size_t iV = 0; iV < intrinsicMesh->nVertices(); iV++) {
    Vertex vB = intrinsicMesh->vertex(iV);
    SurfacePoint pA = vertexLocations[vB];

    switch (pA.type) {
    case SurfacePointType::Vertex: {
      Vertex vA = pA.vertex;

      cs.subdivisionPoints.push_back(
          CommonSubdivisionPoint{CSIntersectionType::VERTEX_VERTEX, pA, SurfacePoint(vB), true});

      aVtx[vA] = &cs.subdivisionPoints[cs.subdivisionPoints.size() - 1];
      break;
    }
    case SurfacePointType::Edge:
      cs.subdivisionPoints.push_back(
          CommonSubdivisionPoint{CSIntersectionType::EDGE_VERTEX, pA, SurfacePoint(vB), true});

      // also add the point to the edge that it crosses
      cs.pointsAlongA[pA.edge].push_back(&cs.subdivisionPoints.back());

      break;
    case SurfacePointType::Face:
      cs.subdivisionPoints.push_back(
          CommonSubdivisionPoint{CSIntersectionType::FACE_VERTEX, pA, SurfacePoint(vB), true});
      break;
    }

    bVtx[vB] = &cs.subdivisionPoints.back();
  }

  // For any input mesh edges which got edge points along them above (aka there have been vertices along the edge), sort


  for (Edge e : intrinsicMesh->edges()) {
    std::vector<SurfacePoint>& edgeTrace = traces[e];

    // First point
    cs.pointsAlongB[e].push_back(bVtx[e.halfedge().tailVertex()]);

    if (edgeIsOriginal[e] || edgeTrace.size() == 2) {
      // For edges which are coincident, don't do anything right now. We will recover these points below.
    } else {
      // General case: record all intersections

      // We don't pass back these from the trace routine, so reconstruct them from the geometry
      std::vector<double> tVals = recoverTraceTValues(edgeTrace);

      for (size_t iP = 1; iP + 1 < edgeTrace.size(); iP++) {
        SurfacePoint& point = edgeTrace[iP];

        if (point.type != SurfacePointType::Edge) {
          throw std::runtime_error(
              "signpost common subdivision construction failure: edge traces should contain only edge crossings");
        }

        // Create new intersection records
        cs.subdivisionPoints.emplace_back();
        CommonSubdivisionPoint& csPoint = cs.subdivisionPoints.back();

        csPoint.posA = point;
        csPoint.posB = SurfacePoint(e, tVals[iP]);
        csPoint.intersectionType = CSIntersectionType::EDGE_TRANSVERSE;
        csPoint.orientation = false; // TODO not properly populated. We don't propagate this info.

        // Store the the pont along the intrinsic edge
        cs.pointsAlongB[e].push_back(&csPoint);

        // Also store them for the input mesh. These will be out of order now, but we will sort them below.
        cs.pointsAlongA[point.edge].push_back(&csPoint);
      }
    }

    // Last point
    cs.pointsAlongB[e].push_back(bVtx[e.halfedge().tipVertex()]);
  }

  // Get all the intersections right along the input mesh edges
  for (Edge e : inputMesh.edges()) {
    std::vector<CommonSubdivisionPoint*>& vec = cs.pointsAlongA[e];

    // Sort the existing points
    std::sort(vec.begin(), vec.end(), [](CommonSubdivisionPoint*& a, CommonSubdivisionPoint*& b) -> bool {
      return a->posA.tEdge < b->posA.tEdge;
    });

    // Prepend the first vertex point to the front
    CommonSubdivisionPoint* firstP = aVtx[e.halfedge().tailVertex()];
    vec.insert(vec.begin(), firstP);

    // Append the last vertex point to the back
    CommonSubdivisionPoint* lastP = aVtx[e.halfedge().tipVertex()];
    vec.push_back(lastP);

    // Walk along the edge, any time we see two consecutive vertex points, there must be a
    // parallel edge there.
    for (size_t iP = 1; iP + 1 < vec.size(); iP++) {

      if (vec[iP - 1]->posB.type == SurfacePointType::Vertex && vec[iP]->posB.type == SurfacePointType::Vertex) {

        // Create a new edge-parallel intersection
        cs.subdivisionPoints.emplace_back();
        CommonSubdivisionPoint& csPoint = cs.subdivisionPoints.back();

        csPoint.posA = SurfacePoint(e, 0.5);
        csPoint.intersectionType = CSIntersectionType::EDGE_PARALLEL;
        csPoint.orientation = false; // TODO not properly populated. We don't propagate this info.

        // Find the matching edge on meshB
        // Note: I'm 90% sure that this is fine, and this search should always find exactly one matching edge.
        // The reasoning why is that the common subdivision is always a simplicial complex, and this segment is
        // and edge of the common subdivision, so it must be unique between its endpoints.
        bool found = false;
        for (Halfedge he : vec[iP - 1]->posB.vertex.outgoingHalfedges()) {
          if (he.tipVertex() == vec[iP]->posB.vertex) {
            csPoint.posB = SurfacePoint(he.edge(), 0.5);
            found = true;
            break;
          }
        }
        GC_SAFETY_ASSERT(found, "could not align shared edges in common subdivision");

        // Add the new point to both lists
        // this list, meshA
        vec.insert(vec.begin() + iP, &csPoint);
        // meshB
        Edge meshBEdge = csPoint.posB.edge;
        GC_SAFETY_ASSERT(cs.pointsAlongB[meshBEdge].size() == 2, "should be exactly two points (endpoints)");
        cs.pointsAlongB[meshBEdge].insert(cs.pointsAlongB[meshBEdge].begin() + 1, &csPoint);
      }
    }
  }
}


} // namespace surface
} // namespace geometrycentral
