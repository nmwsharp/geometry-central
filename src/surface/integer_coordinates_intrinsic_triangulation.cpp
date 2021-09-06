#include "geometrycentral/surface/integer_coordinates_intrinsic_triangulation.h"

#include <ctime>

namespace geometrycentral {
namespace surface {

// helper functions
namespace {
template <typename T>
inline T positivePart(const T& t) {
  return fmax(t, 0.);
}

template <typename T>
inline T negativePart(const T& t) {
  return fmin(t, 0.);
}

inline Vertex src(Edge e) { return e.halfedge().vertex(); }
inline Vertex dst(Edge e) { return e.halfedge().next().vertex(); }

template <typename T>
std::array<T, 3> rotate(const std::array<T, 3>& data) {
  return {data[1], data[2], data[0]};
}

inline Vector3 rotate(const Vector3& v) { return {v.y, v.z, v.x}; }

} // namespace

IntegerCoordinatesIntrinsicTriangulation::IntegerCoordinatesIntrinsicTriangulation(
    ManifoldSurfaceMesh& mesh_, IntrinsicGeometryInterface& inputGeom_, double mollifyEPS)
    : IntrinsicTriangulation(mesh_, inputGeom_), normalCoordinates(*intrinsicMesh) {

  normalCoordinates.setCurvesFromEdges(*intrinsicMesh);

  // TODO document/expose this somehow, rather than just doing it silently
  if(mollifyEPS > 0) {
    mollifyIntrinsic(*intrinsicMesh, edgeLengths, mollifyEPS);
  }

  vertexLocations = VertexData<SurfacePoint>(*intrinsicMesh);
  for (size_t iV = 0; iV < intrinsicMesh->nVertices(); iV++) {
    Vertex vIntrinsic = intrinsicMesh->vertex(iV);
    Vertex vInput = inputMesh.vertex(iV);
    vertexLocations[vIntrinsic] = vInput;
  }
}

// ======================================================
//                 Queries & Accesses
// ======================================================

EdgeData<std::vector<SurfacePoint>> IntegerCoordinatesIntrinsicTriangulation::traceAllIntrinsicEdgesAlongInput() {

  CommonSubdivision& cs = getCommonSubdivision();

  EdgeData<std::vector<SurfacePoint>> tracedEdges(*intrinsicMesh);
  for (Edge eB : intrinsicMesh->edges()) {
    for (CommonSubdivisionPoint* pt : cs.pointsAlongB[eB]) {
      if (pt->intersectionType != CSIntersectionType::EDGE_PARALLEL) tracedEdges[eB].push_back(pt->posA);
    }
  }

  return tracedEdges;
}

std::vector<SurfacePoint>
IntegerCoordinatesIntrinsicTriangulation::traceIntrinsicHalfedgeAlongInput(Halfedge intrinsicHe) {
  CommonSubdivision& cs = getCommonSubdivision();

  std::vector<SurfacePoint> trajectory;
  if (intrinsicHe == intrinsicHe.edge().halfedge()) {
    // Halfedge points in same direction as edge, so points are already in correct order
    for (CommonSubdivisionPoint* pt : cs.pointsAlongB[intrinsicHe.edge()]) {
      if (pt->intersectionType != CSIntersectionType::EDGE_PARALLEL) trajectory.push_back(pt->posA);
    }
  } else {
    // Halfedge points in opposite direction as edge, so order must be reversed
    const std::vector<CommonSubdivisionPoint*>& pts = cs.pointsAlongB[intrinsicHe.edge()];
    for (auto rit = pts.rbegin(); rit != pts.rend(); ++rit) { // loop backwards
      if ((*rit)->intersectionType != CSIntersectionType::EDGE_PARALLEL) trajectory.push_back((*rit)->posA);
    }
  }

  return trajectory;
}

EdgeData<std::vector<SurfacePoint>> IntegerCoordinatesIntrinsicTriangulation::traceAllInputEdgesAlongIntrinsic() {
  // Might as well compute full common subdivision
  CommonSubdivision& cs = getCommonSubdivision();

  EdgeData<std::vector<SurfacePoint>> tracedEdges(inputMesh);
  for (Edge eA : inputMesh.edges()) {
    for (CommonSubdivisionPoint* pt : cs.pointsAlongA[eA]) {
      if (pt->intersectionType != CSIntersectionType::EDGE_PARALLEL) tracedEdges[eA].push_back(pt->posB);
    }
  }
  return tracedEdges;
}

std::vector<SurfacePoint> IntegerCoordinatesIntrinsicTriangulation::traceInputHalfedgeAlongIntrinsic(Halfedge inputHe) {
  if (commonSubdivision) {
    // If the common subdivision has already been computed, just read off values

    CommonSubdivision& cs = getCommonSubdivision();

    std::vector<SurfacePoint> trajectory;
    if (inputHe == inputHe.edge().halfedge()) {
      // Halfedge points in same direction as edge, so points are already in correct order
      for (CommonSubdivisionPoint* pt : cs.pointsAlongA[inputHe.edge()]) {
        if (pt->intersectionType != CSIntersectionType::EDGE_PARALLEL) trajectory.push_back(pt->posB);
      }
    } else {
      // Halfedge points in opposite direction as edge, so order must be reversed
      const std::vector<CommonSubdivisionPoint*>& pts = cs.pointsAlongA[inputHe.edge()];
      for (auto rit = pts.rbegin(); rit != pts.rend(); ++rit) { // loop backwards
        if ((*rit)->intersectionType != CSIntersectionType::EDGE_PARALLEL) trajectory.push_back((*rit)->posB);
      }
    }

    return trajectory;

  } else {
    // Otherwise do the tracing manually
    std::vector<SurfacePoint> trajectory;

    NormalCoordinatesCompoundCurve ncCurves = traceInputHalfedge(inputHe);

    for (size_t iC = 0; iC < ncCurves.components.size(); iC++) {
      const NormalCoordinatesCurve& curve = ncCurves.components[iC];
      if (iC == 0) {
        // Assumes that the vertices of inputMesh all appear in intrinsicMesh with the same indices
        Vertex intrinsicTail = intrinsicMesh->vertex(inputHe.tailVertex().getIndex());
        trajectory.push_back(SurfacePoint(intrinsicTail));
      }


      auto& path = curve.crossings;
      if (path.size() == 1 && std::get<0>(path[0]) < 0) { // shared edge
        Halfedge heB = std::get<1>(path[0]);
        trajectory.push_back(SurfacePoint(heB.tipVertex()));
      } else { // transverse crossings
        std::vector<std::pair<SurfacePoint, double>> geodesicPath =
            generateFullSingleGeodesicGeometry(*intrinsicMesh, *this, curve);

        Halfedge hB;
        for (size_t iC = 0; iC < path.size(); ++iC) {
          hB = std::get<1>(path[iC]);
          // geodesicPath stores the start and end point, which
          // path doesn't do, so we offset by 1 here
          SurfacePoint ptB = std::get<0>(geodesicPath[iC + 1]);
          trajectory.push_back(ptB);
        }

        Vertex bDst = hB.twin().next().tipVertex();
        trajectory.push_back(SurfacePoint(bDst));
      }
    }

    return trajectory;
  }
}

SurfacePoint IntegerCoordinatesIntrinsicTriangulation::equivalentPointOnIntrinsic(const SurfacePoint& pointOnInput) {
  // TODO
  throw std::runtime_error("not implemented");
  return SurfacePoint(Vertex());
}

SurfacePoint IntegerCoordinatesIntrinsicTriangulation::equivalentPointOnInput(const SurfacePoint& pointOnIntrinsic) {
  switch (pointOnIntrinsic.type) {
  case SurfacePointType::Vertex:
    return vertexLocations[pointOnIntrinsic.vertex];
  case SurfacePointType::Edge: {
    SurfacePoint facePoint = pointOnIntrinsic.inSomeFace();
    return computeFaceSplitData(facePoint.face, facePoint.faceCoords).first;
  }
  case SurfacePointType::Face:
    return computeFaceSplitData(pointOnIntrinsic.face, pointOnIntrinsic.faceCoords).first;
  }
  return SurfacePoint(); // unreachable
}


void IntegerCoordinatesIntrinsicTriangulation::constructCommonSubdivision() {

  intrinsicMesh->compress();

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
      break;
    case SurfacePointType::Face:
      cs.subdivisionPoints.push_back(
          CommonSubdivisionPoint{CSIntersectionType::FACE_VERTEX, pA, SurfacePoint(vB), true});
      break;
    }

    bVtx[vB] = &cs.subdivisionPoints.back();
  }

  // Allocate space for points along B's edges, and also fill in the endpoints
  // of B's edges
  // Note that for shared edges, we will actually need 3 entries instead of 2
  // (due to a HACK), but we'll fix that later
  for (Edge eB : intrinsicMesh->edges()) {
    int n = positivePart(normalCoordinates[eB]);
    cs.pointsAlongB[eB].resize(n + 2);
    cs.pointsAlongB[eB][0] = bVtx[src(eB)];
    cs.pointsAlongB[eB][n + 1] = bVtx[dst(eB)];
  }

  // Trace the edges of mesh A (inputMesh) over mesh B (intrinsicMesh)
  for (Edge eA : inputMesh.edges()) {
    NormalCoordinatesCompoundCurve compoundPath = traceInputEdge(eA);

    for (size_t iC = 0; iC < compoundPath.components.size(); iC++) {

      const NormalCoordinatesCurve& curve = compoundPath.components[iC];
      bool first = iC == 0;

      auto& path = curve.crossings;

      if (path.size() == 1 && std::get<0>(path[0]) < 0) {

        // Shared edge
        Halfedge heB = std::get<1>(path[0]);
        Edge eB = heB.edge();
        bool positiveOrientation = heB == eB.halfedge();

        GC_SAFETY_ASSERT(normalCoordinates[eB] < 0,
                         "eB should have negative n, but has n = " + std::to_string(normalCoordinates[eB]));

        // make pointsAlongB[eB] one longer
        cs.pointsAlongB[eB].push_back(cs.pointsAlongB[eB][cs.pointsAlongB[eB].size() - 1]);

        // Construct intersection point, flagging it with barycentric
        // coordinate 0.5 to indicate that the edges are parallel
        cs.subdivisionPoints.push_back(CommonSubdivisionPoint{CSIntersectionType::EDGE_PARALLEL, SurfacePoint(eA, 0.5),
                                                              SurfacePoint(heB.edge(), 0.5), positiveOrientation});

        CommonSubdivisionPoint* crPt = &cs.subdivisionPoints[cs.subdivisionPoints.size() - 1];

        if (first) cs.pointsAlongA[eA].push_back(aVtx[src(eA)]);
        cs.pointsAlongA[eA].push_back(crPt);
        cs.pointsAlongA[eA].push_back(bVtx[dst(eB)]);

        cs.pointsAlongB[eB][1] = crPt;
      } else {
        // Indirect path

        // if (geodesic) { // Lay out geodesic
        // std::cout << path << std::endl;

        std::vector<std::pair<SurfacePoint, double>> geodesicPath =
            generateFullSingleGeodesicGeometry(*intrinsicMesh, *this, curve);

        if (first) cs.pointsAlongA[eA].push_back(aVtx[src(eA)]);

        Halfedge hB;
        for (size_t iC = 0; iC < path.size(); ++iC) {
          hB = std::get<1>(path[iC]);
          Edge eB = hB.edge();
          bool positiveOrientation = hB == eB.halfedge();

          // geodesicPath stores the start and end point, which
          // path doesn't do, so we offset by 1 here
          SurfacePoint ptB = std::get<0>(geodesicPath[iC + 1]);
          double tA = std::get<1>(geodesicPath[iC + 1]);

          int iB = std::get<0>(path[iC]);
          if (!positiveOrientation) iB = positivePart(normalCoordinates[eB]) - iB - 1;

          cs.subdivisionPoints.push_back(CommonSubdivisionPoint{CSIntersectionType::EDGE_TRANSVERSE,
                                                                SurfacePoint(eA, tA), ptB, positiveOrientation});

          CommonSubdivisionPoint* crPt = &cs.subdivisionPoints[cs.subdivisionPoints.size() - 1];

          cs.pointsAlongA[eA].push_back(crPt);
          cs.pointsAlongB[eB][iB + 1] = crPt;
        }

        Vertex bDst = hB.twin().next().tipVertex();
        cs.pointsAlongA[eA].push_back(bVtx[bDst]);
        //
        /* Disabled option: space points evenly
         *
         * } else {

          //std::cout << "PATH" << endl;

          if (first) cs.pointsAlongA[eA].push_back(aVtx[src(eA)]);
          Halfedge hB;
          for (const std::pair<int, Halfedge>& pt : path) {
            hB = std::get<1>(pt);
            Edge eB = hB.edge();
            bool positiveOrientation = hB == eB.halfedge();

            int iB = std::get<0>(pt);
            if (!positiveOrientation) iB = positivePart(normalCoordinates[eB]) - iB - 1;

            double tB = (iB + 1) / (double)(positivePart(normalCoordinates[eB]) + 1);

            cs.subdivisionPoints.push_back(CommonSubdivisionPoint{CSIntersectionType::EDGE_TRANSVERSE,
                                                                  SurfacePoint(eA, 0.5), SurfacePoint(hB.edge(), tB),
                                                                  positiveOrientation});

            CommonSubdivisionPoint* crPt = &cs.subdivisionPoints[cs.subdivisionPoints.size() - 1];

            cs.pointsAlongA[eA].push_back(crPt);
            cs.pointsAlongB[eB][iB + 1] = crPt;
          }
          Vertex bDst = hB.twin().next().tipVertex();
          cs.pointsAlongA[eA].push_back(bVtx[bDst]);

          // Space out points on mesh A
          for (size_t iA = 1; iA + 1 < cs.pointsAlongA[eA].size(); iA++) {
            double tA = iA / (double)(cs.pointsAlongA[eA].size() - 1);
            cs.pointsAlongA[eA][iA]->posA.tEdge = tA;
          }
        }
        */
      }
    }
  }

  // Check that we've accounted for all promised crossings along B's edges
  for (Edge eB : intrinsicMesh->edges()) {
    for (const auto& p : cs.pointsAlongB[eB]) {
      if (!p) {
        std::cout << "Missed a crossing on edge " << eB << " which is supposed to have " << normalCoordinates[eB]
                  << std::endl;
        std::cout << "Endpoints " << vertexLocations[eB.halfedge().tailVertex()] << std::endl
                  << "\t->" << vertexLocations[eB.halfedge().tipVertex()] << std::endl;
        GC_SAFETY_ASSERT(p, "oops, missed a crossing");
      }
    }
    for (size_t iC = 1; iC + 1 < cs.pointsAlongB[eB].size(); ++iC) {
      if (cs.pointsAlongB[eB][iC]->intersectionType == CSIntersectionType::VERTEX_VERTEX) {
        std::cout << "Error at crossing " << iC << " of " << cs.pointsAlongB[eB].size() << " on edge " << eB
                  << std::endl;
        for (const auto& pt : cs.pointsAlongB[eB]) std::cout << "\t" << *pt << std::endl;
        throw std::runtime_error("encountered vertex intersection "
                                 "in the middle of an edge ?!");
      }
    }
  }
}

// ======================================================
//                Low-Level Mutators
// ======================================================

// If the edge is not Delaunay, flip it. Returns true if flipped.
// TODO
bool IntegerCoordinatesIntrinsicTriangulation::flipEdgeIfNotDelaunay(Edge e) {
  // Can't flip
  if (isFixed(e)) return false;

  // Don't need to
  if (isDelaunay(e)) return false;

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
  auto nUpdate = normalCoordinates.computeFlippedData(e);
  bool flipped = intrinsicMesh->flip(e, false);

  // Might not have been flippable for connectivity reasons
  if (!flipped) {
    return false;
  }

  // Assign the new edge lengths
  // TODO project to satisfy triangle inequality?
  edgeLengths[e] = newLength;
  normalCoordinates.applyFlippedData(e, nUpdate);

  // === Update various quantities

  // depends on edgeLengths
  updateFaceBasis(e.halfedge().face());
  updateFaceBasis(e.halfedge().twin().face());

  std::array<Corner, 6> incidentCorners{
      e.halfedge().corner(),        e.halfedge().next().corner(),        e.halfedge().next().next().corner(),
      e.halfedge().twin().corner(), e.halfedge().twin().next().corner(), e.halfedge().twin().next().next().corner()};

  std::array<Vertex, 4> incidentVertices{e.halfedge().vertex(), e.halfedge().next().vertex(),
                                         e.halfedge().next().next().vertex(),
                                         e.halfedge().twin().next().next().vertex()};
  for (Vertex v : incidentVertices) {
    // depends on vertexAngleSums
    updateHalfedgeVectorsInVertex(v);
  }

  // Do callbacks
  triangulationChanged();
  invokeEdgeFlipCallbacks(e);

  return true;
}

// If the edge can be flipped, flip it (must be combinatorially flippable
// and inside a convex quad). Returns true if flipped.
// TODO
bool IntegerCoordinatesIntrinsicTriangulation::flipEdgeIfPossible(Edge e) {
  // Can't flip
  if (isFixed(e)) return false;

  // Get geometric data
  Halfedge he = e.halfedge();
  std::array<Vector2, 4> layoutPositions = layoutDiamond(he);

  // Test if geometryically flippable flippable (both signed areas of new
  // triangles are positive)
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
  auto nUpdate = normalCoordinates.computeFlippedData(e);
  bool flipped = intrinsicMesh->flip(e, false);

  // Might not have been flippable for connectivity reasons
  if (!flipped) {
    return false;
  }

  // Assign the new edge lengths
  // TODO project to satisfy triangle inequality?
  edgeLengths[e] = newLength;
  normalCoordinates.applyFlippedData(e, nUpdate);

  // === Update various quantities

  // depends on edgeLengths
  updateFaceBasis(e.halfedge().face());
  updateFaceBasis(e.halfedge().twin().face());

  std::array<Corner, 6> incidentCorners{
      e.halfedge().corner(),        e.halfedge().next().corner(),        e.halfedge().next().next().corner(),
      e.halfedge().twin().corner(), e.halfedge().twin().next().corner(), e.halfedge().twin().next().next().corner()};

  std::array<Vertex, 4> incidentVertices{e.halfedge().vertex(), e.halfedge().next().vertex(),
                                         e.halfedge().next().next().vertex(),
                                         e.halfedge().twin().next().next().vertex()};
  for (Vertex v : incidentVertices) {
    // depends on vertexAngleSums
    updateHalfedgeVectorsInVertex(v);
  }


  // Do callbacks
  triangulationChanged();
  invokeEdgeFlipCallbacks(e);

  return true;
}



double IntegerCoordinatesIntrinsicTriangulation::checkFlip(Edge e) {
  // Can't flip
  if (isFixed(e)) return std::numeric_limits<double>::infinity();

  // Check topologically flippable
  {
    Halfedge ha1 = e.halfedge();
    Halfedge ha2 = ha1.next();
    Halfedge ha3 = ha2.next();
    Halfedge hb1 = ha1.sibling();
    Halfedge hb2 = hb1.next();
    Halfedge hb3 = hb2.next();

    // incident on degree 1 vertex
    if (ha2 == hb1 || hb2 == ha1) {
      return std::numeric_limits<double>::infinity();
    }
  }

  // Get geometric data
  Halfedge he = e.halfedge();
  std::array<Vector2, 4> layoutPositions = layoutDiamond(he);

  // Test if geometryically flippable flippable (both signed areas of new
  // triangles are positive)
  double A1 = cross(layoutPositions[1] - layoutPositions[0], layoutPositions[3] - layoutPositions[0]);
  double A2 = cross(layoutPositions[3] - layoutPositions[2], layoutPositions[1] - layoutPositions[2]);
  double areaSum = (A1 + A2);

  return std::min(A1 / areaSum, A2 / areaSum);
}

Vertex IntegerCoordinatesIntrinsicTriangulation::insertCircumcenterOrSplitSegment(Face f, bool verbose) {
  // === Circumcenter in barycentric coordinates

  Halfedge he0 = f.halfedge();
  double a = edgeLengths[he0.next().edge()];
  double b = edgeLengths[he0.next().next().edge()];
  double c = edgeLengths[he0.edge()];
  double a2 = a * a;
  double b2 = b * b;
  double c2 = c * c;
  Vector3 circumcenterLoc = {a2 * (b2 + c2 - a2), b2 * (c2 + a2 - b2), c2 * (a2 + b2 - c2)};
  circumcenterLoc = normalizeBarycentric(circumcenterLoc);

  // Trace from the barycenter (have to trace from somewhere)
  Vector3 barycenter = Vector3::constant(1. / 3.);
  Vector3 vecToCircumcenter = circumcenterLoc - barycenter;

  // === Trace the ray to find the location of the new point on the intrinsic
  // meshes

  // Data we need from the intrinsic trace
  TraceOptions options;
  if (markedEdges.size() > 0) {
    options.barrierEdges = &markedEdges;
  }
  TraceGeodesicResult intrinsicTraceResult = traceGeodesic(*this, f, barycenter, vecToCircumcenter, options);
  SurfacePoint newPositionOnIntrinsic = intrinsicTraceResult.endPoint;
  if (newPositionOnIntrinsic.type == SurfacePointType::Vertex && newPositionOnIntrinsic.vertex == Vertex()) {
    // tracing failed
    return Vertex();
  }

  // If the circumcenter is blocked by an edge, insert the midpoint of that
  // edge instead (which happens to be just what is needed for Chew's 2nd
  // algo).
  // TODO: also do this if we hit locked edges
  // TODO: delete inserted vertices inside the diametral circle?
  if (newPositionOnIntrinsic.type == SurfacePointType::Edge) {
    newPositionOnIntrinsic.tEdge = 0.5;
  }

  if (verbose) {
    std::cout << "insertion point: " << newPositionOnIntrinsic << std::endl;
  }

  // === Phase 3: Add the new vertex
  Vertex newV = insertVertex(newPositionOnIntrinsic);
  return newV;
}


// Assumes vertexAngleSums exist and are up to date
void IntegerCoordinatesIntrinsicTriangulation::updateHalfedgeVectorsInVertex(Vertex v) {
  // stolen from intrinsic_geometry_interface.cpp

  auto cornerScaledAngle = [&](Corner c) -> double {
    if (c.vertex().isBoundary()) {
      double s = PI / vertexAngleSums[c.vertex()];
      return s * cornerAngle(c);
    } else {
      double s = 2.0 * PI / vertexAngleSums[c.vertex()];
      return s * cornerAngle(c);
    }
  };

  double coordSum = 0.0;

  // Custom loop to orbit CCW
  Halfedge firstHe = v.halfedge();
  Halfedge currHe = firstHe;
  do {
    halfedgeVectorsInVertex[currHe] = Vector2::fromAngle(coordSum) * edgeLengths[currHe.edge()];
    if (!currHe.isInterior()) break;
    coordSum += cornerScaledAngle(currHe.corner());
    currHe = currHe.next().next().twin();
  } while (currHe != firstHe);
}

// ======================================================
//                Low-Level Queries
// ======================================================


// Takes in a halfedge of the intrinsic mesh whose edge's normal coordinate
// is negative (meaning that it lies along an edge of the input mesh) and
// returns the halfedge in the input mesh pointing in the same direction
// e.vertex() must live in both meshes
Halfedge IntegerCoordinatesIntrinsicTriangulation::getSharedInputEdge(Halfedge he) const {
  GC_SAFETY_ASSERT(vertexLocations[he.tailVertex()].type == SurfacePointType::Vertex,
                   "I can only identify edges which come out of shared vertices");

  int iE = normalCoordinates.roundabouts[he];
  while (iE < 0) iE += normalCoordinates.roundaboutDegrees[he.vertex()];

  Vertex inputVertex = vertexLocations[he.tailVertex()].vertex;
  Halfedge inputHe = inputVertex.halfedge();

  // step counterclockise iE times
  for (int i = 0; i < iE; i++) inputHe = inputHe.next().next().twin();

  return inputHe;
}

std::pair<SurfacePoint, std::array<int, 3>>
IntegerCoordinatesIntrinsicTriangulation::computeFaceSplitData(Face f, Vector3 bary, bool verbose) {
  Face insertionFace;
  Vector3 insertionBary;
  std::array<int, 3> counts;

  if (normalCoordinates[f.halfedge().edge()] < 0 && normalCoordinates[f.halfedge().next().edge()] < 0 &&
      normalCoordinates[f.halfedge().next().next().edge()] < 0) {

    // Case 0a: face is empty, all edges shared
    // Note that this means all vertices must be shared
    //

    Halfedge inputHalfedge = identifyInputEdge(f.halfedge());
    insertionFace = inputHalfedge.face();
    insertionBary = bary;

    if (inputHalfedge == insertionFace.halfedge()) {
      // good!
    } else if (inputHalfedge == insertionFace.halfedge().next()) {
      insertionBary = rotate(rotate(insertionBary));
    } else if (inputHalfedge == insertionFace.halfedge().next().next()) {
      insertionBary = rotate(insertionBary);
    } else {
      throw std::runtime_error("Face " + std::to_string(insertionFace) + " is not triangular");
    }
    counts = {0, 0, 0};
  } else if (normalCoordinates[f.halfedge().edge()] <= 0 && normalCoordinates[f.halfedge().next().edge()] <= 0 &&
             normalCoordinates[f.halfedge().next().next().edge()] <= 0) {

    // TODO: use getParentFace
    insertionFace = Face();
    // First, look for a FacePoint - that's the easiest thing to deal with
    for (Vertex v : f.adjacentVertices()) {
      if (vertexLocations[v].type == SurfacePointType::Face) {
        if (insertionFace != Face()) {
          // all faces should agree
          GC_SAFETY_ASSERT(insertionFace == vertexLocations[v].face, "if there are no crossings, this face must "
                                                                     "be contained inside a single input face");
        } else {
          insertionFace = vertexLocations[v].face;
        }
      }
    }

    if (insertionFace == Face()) {
      // If we couldn't find a FacePoint, look for an EdgePoint
      // Note that not all vertices can be shared - if they are, then we
      // can't have any boundary crossings. But this happens sometimes
      // anyway, and is handled next

      auto containsVertex = [](Face f, Vertex v) -> bool {
        for (Vertex vF : f.adjacentVertices()) {
          if (vF == v) return true;
        }
        return false;
      };

      auto containsEdge = [](Face f, Edge e) -> bool {
        for (Edge eF : f.adjacentEdges()) {
          if (eF == e) return true;
        }
        return false;
      };

      auto compatible = [&](const SurfacePoint& pt, Face f) -> bool {
        switch (pt.type) {
        case SurfacePointType::Vertex:
          return containsVertex(f, pt.vertex);
        case SurfacePointType::Edge:
          return containsEdge(f, pt.edge);
        case SurfacePointType::Face:
          return pt.face == f;
        }
        return false; // unreachable
      };

      for (Vertex v : f.adjacentVertices()) {
        if (vertexLocations[v].type == SurfacePointType::Edge) {
          Edge e = vertexLocations[v].edge;
          Face f1 = e.halfedge().face();
          Face f2 = e.halfedge().twin().face();

          bool f1Okay = e.halfedge().isInterior();
          bool f2Okay = e.halfedge().twin().isInterior();

          for (Vertex w : f.adjacentVertices()) {
            f1Okay = f1Okay && compatible(vertexLocations[w], f1);
            f2Okay = f2Okay && compatible(vertexLocations[w], f2);
          }

          if ((f1 != f2) && (f1Okay && f2Okay)) {
            std::cerr << "splitFace err: There are two options for an "
                         "input "
                         "face to insert "
                         "in, and they both look fine. What should I "
                         "do?. (In the meantime, I'm picking option 1)"
                      << std::endl;
          }

          if (f1Okay) {
            insertionFace = f1;
          } else if (f2Okay) {
            insertionFace = f2;
          } else {
            std::cout << " ==== Failed to find valid face to insert" << std::endl;
            std::cout << "f1: " << f1 << " (isBoundaryLoop: " << f1.isBoundaryLoop() << ")" << std::endl;
            std::cout << "f2: " << f2 << " (isBoundaryLoop: " << f2.isBoundaryLoop() << ")" << std::endl;
            std::cout << "Underlying edge: " << e << " (isBoundary: " << e.isBoundary() << ")" << std::endl;

            throw std::runtime_error("Could not find a valid face in inputMesh to "
                                     "insert the new vertex into");
          }
        }
      }
    }

    if (insertionFace == Face()) {
      // In really degenerate situations, we can have triangles between
      // three
      bool allSharedVertices = true;
      for (Vertex v : f.adjacentVertices()) {
        allSharedVertices = allSharedVertices && vertexLocations[v].type == SurfacePointType::Vertex;
      }
      if (allSharedVertices) {
        // Hope that we find a shared edge. If not, something terrible
        // has happened
        std::cout << vertexLocations[intrinsicMesh->vertex(2336)] << std::endl;
        for (Halfedge he : f.adjacentHalfedges()) {
          std::cout << std::endl;
        }

        for (Halfedge he : f.adjacentHalfedges()) {
          if (normalCoordinates[he.edge()] < 0) {
            Halfedge inputHalfedge = identifyInputEdge(he);
            insertionFace = inputHalfedge.face();
            insertionBary = bary;

            if (inputHalfedge == insertionFace.halfedge()) {
              // good!
            } else if (inputHalfedge == insertionFace.halfedge().next()) {
              insertionBary = rotate(rotate(insertionBary));
            } else if (inputHalfedge == insertionFace.halfedge().next().next()) {
              insertionBary = rotate(insertionBary);
            } else {
              throw std::runtime_error("Face " + std::to_string(insertionFace) + " is not triangular");
            }
            counts = {0, 0, 0};
            break;
          }
        }
      }
    }

    GC_SAFETY_ASSERT(insertionFace != Face(), "failed to find FacePoint.");

    std::array<Vector3, 3> vertexBary;
    size_t iV = 0;
    for (Vertex v : f.adjacentVertices()) {
      vertexBary[iV] = vertexLocations[v].inFace(insertionFace).faceCoords;
      iV++;
    }

    insertionBary = bary.x * vertexBary[0] + bary.y * vertexBary[1] + bary.z * vertexBary[2];

    if (verbose) {
      for (Vertex v : f.adjacentVertices()) {
        std::cout << vertexLocations[v] << std::endl;
        std::cout << "\t" << vertexLocations[v].inFace(insertionFace) << std::endl;
      }
    }
    counts = {0, 0, 0};

  } else {
    // Populate the crossing locations for the edges of the triangle
    size_t iHe = 0;
    std::array<std::vector<NormalCoordinatesCurve>, 3> curves;
    std::array<std::vector<std::vector<std::pair<SurfacePoint, double>>>, 3> geodesics;
    std::array<std::vector<double>, 3> transverseCrossingTimes;
    std::array<std::vector<double>, 3> boundaryCrossings;
    std::array<std::vector<size_t>, 3> crossingID;
    for (Halfedge he : f.adjacentHalfedges()) {
      for (int ind = 0; ind < normalCoordinates[he.edge()]; ind++) {

        // Get the topological crossings for the curve
        NormalCoordinatesCurve crossings;
        int centerCrossInd;
        std::tie(crossings, centerCrossInd) = normalCoordinates.topologicalTraceBidirectional(he, ind);

        curves[iHe].push_back(crossings);
        geodesics[iHe].push_back(generateFullSingleGeodesicGeometry(*intrinsicMesh, *this, crossings));
        std::vector<std::pair<SurfacePoint, double>>& geodesic = geodesics[iHe].back();
        SurfacePoint& thisCross = geodesic[centerCrossInd + 1].first;
        double tCross = (he.edge().halfedge() == he) ? thisCross.tEdge : 1 - thisCross.tEdge;
        boundaryCrossings[iHe].push_back(tCross);
        transverseCrossingTimes[iHe].push_back(geodesic[centerCrossInd + 1].second);
        crossingID[iHe].push_back(centerCrossInd + 1);
      }
      iHe++;
    }

    // tStart = std::clock();

    // TODO apply some sanity policies, like that the crossings should be
    // correctly ordered

    counts = computeVertexInsertionCrossingCounts(bary, boundaryCrossings);

    // Find a halfedge bounding the region that our inserted vertex will
    // live inserted in
    std::array<Corner, 3> faceCorners{f.halfedge().corner(), f.halfedge().next().corner(),
                                      f.halfedge().next().next().corner()};

    auto heIndex = [](Halfedge he) -> int {
      if (he == he.face().halfedge()) {
        return 0;
      } else if (he == he.face().halfedge().next()) {
        return 1;
      } else {
        return 2;
      }
    };

    auto next = [](int i) { return (i + 1) % 3; };
    auto heBary = [&](Halfedge he, double t) -> Vector3 {
      int i = heIndex(he);
      int j = next(i);

      Vector3 bary = Vector3::zero();
      bary[i] = (1 - t);
      bary[j] = (t);

      return bary;
    };

    // Compute the intrinsic and input coordinates of a crossing
    // along a boundary halfedge of intrinsic face f.
    // If pos == 0, returns the intrinsic and input coordinates for the
    // vertex
    auto computeIntrinsicAndInputPoints = [&](Halfedge he, int pos) -> std::pair<Vector3, SurfacePoint> {
      if (pos == 0) {
        // source vertex
        return {heBary(he, 0), vertexLocations[he.vertex()]};
      } else if (pos == positivePart(normalCoordinates[he.edge()]) + 1) {
        // target vertex
        return {heBary(he, 1), vertexLocations[he.next().vertex()]};
      } else {
        // intermediate crossing
        double tInput = transverseCrossingTimes[heIndex(he)][pos - 1];
        Halfedge heInput = identifyInputEdge(curves[heIndex(he)][pos - 1]);
        double tIntrinsic = boundaryCrossings[heIndex(he)][pos - 1];
        return {heBary(he, tIntrinsic), SurfacePoint(heInput, tInput)};
      }
    };

    // Identify 2 points in the input face and their corresponding positions
    // in the intrinsic face
    // TODO: pick 2 points to be spread out. We get numerical errors if the
    // chosen points are too close together

    // Identify as many points as possible on the input face and their
    // corresponding positions on the intrinsic face. In theory 2 suffice,
    // but using more improves stability

    Face inputFace;
    std::vector<std::pair<Vector3, SurfacePoint>> intrinsicInputPairs;
    bool foundRegion = false;

    // Used to print error messages later
    size_t myCase = 0;

    // Case 1: in a corner, or just past a corner
    for (size_t iC = 0; iC < 3 && !foundRegion; iC++) {
      int cornerCoord = static_cast<int>(normalCoordinates.strictCornerCoord(faceCorners[iC]));

      if (cornerCoord < 1 || counts[iC] > cornerCoord) continue;

      if (verbose) std::cout << "Case I" << std::endl;
      myCase = 1;

      bool justPast = (counts[iC] == static_cast<int>(normalCoordinates.strictCornerCoord(faceCorners[iC])));

      // new vertex is contained in corner iC
      // Bounding halfedge is the counts[iC]'th crossing along
      // corner.halfedge()

      Halfedge firstHedge = faceCorners[iC].halfedge();
      Halfedge secondHedge = firstHedge.next().next();

      // Grab "outside" halfedge
      int firstInputPointIndex = counts[iC];
      if (justPast) firstInputPointIndex -= 1;
      int secondInputPointIndex = normalCoordinates[secondHedge.edge()] - (firstInputPointIndex + 1);

      auto& inputHedgeCurve = curves[heIndex(firstHedge)][firstInputPointIndex];

      // Identify halfedge in input mesh
      // TODO: repeated effort in computeIntrinsicAndInputPoints
      // Have to take twin due to weird orientation conventions
      Halfedge inputHalfedge = identifyInputEdge(inputHedgeCurve).twin();

      if (justPast) { // flip orientation
        inputHalfedge = inputHalfedge.twin();
      }

      inputFace = inputHalfedge.face();

      intrinsicInputPairs.push_back(computeIntrinsicAndInputPoints(firstHedge, firstInputPointIndex + 1));
      intrinsicInputPairs.push_back(computeIntrinsicAndInputPoints(secondHedge, secondInputPointIndex + 1));

      if (justPast) {
        // If we're just past the end of a corner cell, add in the next
        // crossing along firstHedge and the previous crossing along
        // secondHedge
        // TODO: special case for central polygon?
        intrinsicInputPairs.push_back(computeIntrinsicAndInputPoints(firstHedge, firstInputPointIndex + 2));
        intrinsicInputPairs.push_back(computeIntrinsicAndInputPoints(secondHedge, secondInputPointIndex));
      } else {
        // If we're in a corner cell, add in the previous crossing along
        // firstHedge
        intrinsicInputPairs.push_back(computeIntrinsicAndInputPoints(firstHedge, firstInputPointIndex));

        // If the previous crossing along firstHedge wasn't a vertex,
        // we're in a quad and can also add the next crossing along
        // secondHedge
        if (firstInputPointIndex > 0) {
          intrinsicInputPairs.push_back(computeIntrinsicAndInputPoints(secondHedge, secondInputPointIndex + 2));
        }
      }

      for (size_t iP = 0; iP < intrinsicInputPairs.size(); iP++) {
        if (!checkAdjacent(intrinsicInputPairs[iP].second, SurfacePoint(inputFace, Vector3::zero()))) {
          std::cout << "Case: I, justPast: " << justPast << std::endl;
          std::cout << "\tfirstInputPointIndex: " << firstInputPointIndex << std::endl;
          std::cout << "\tsecondInputPointIndex: " << secondInputPointIndex << std::endl;
          std::cout << "\tpoint: " << intrinsicInputPairs[iP].second << std::endl;
          std::cout << "\tface: " << inputFace << std::endl;
          std::cout << "\tiP: " << iP << std::endl;
          for (size_t jP = 0; jP < intrinsicInputPairs.size(); jP++) {
            std::cout << "\t\t point " << jP << " : " << intrinsicInputPairs[jP].second << std::endl;
          }
        }
      }

      foundRegion = true;
    }

    // Case 2: in a corner whose corner coordinate is 0
    // Note that we must be in the fan configuration. Otherwise we would
    // also be one-past the final curve crossing some opposite corner
    for (int iC = 0; iC < 3 && !foundRegion; iC++) {
      if (counts[iC] != 0) continue;

      Halfedge firstHedge;
      bool hasFanEdge = normalCoordinates.triangleInequalityViolation(f, firstHedge);

      if (!hasFanEdge) {
        for (Edge e : f.adjacentEdges()) {
          std::cout << normalCoordinates[e] << std::endl;
        }
        GC_SAFETY_ASSERT(hasFanEdge, "Triangles in Case 2 must be fan triangles");
      }
      int longHe = heIndex(firstHedge);

      if (iC == next(longHe)) {
        // Covered by case 3
        continue;
      } else if (iC == next(next(longHe))) {
        // next(next(longHe)) is the fan corner, which doesn't work.
        // One of the other corners must be better
        continue;
      }
      if (verbose) std::cout << "Case II" << std::endl;
      myCase = 2;

      // We must have iC == longHe
      // So the new point is in the first section along longHe

      auto& inputHedgeCurve = curves[heIndex(firstHedge)][0];

      // Identify halfedge in input mesh
      Halfedge inputHalfedge = identifyInputEdge(inputHedgeCurve).twin();

      inputFace = inputHalfedge.face();

      intrinsicInputPairs.push_back(computeIntrinsicAndInputPoints(firstHedge, 1));
      intrinsicInputPairs.push_back(computeIntrinsicAndInputPoints(firstHedge, 0));
      intrinsicInputPairs.push_back(computeIntrinsicAndInputPoints(firstHedge.next().next(), 0));

      foundRegion = true;
    }

    // Case 3: in the center of the fan region
    if (!foundRegion) {
      if (verbose) std::cout << "Case III" << std::endl;
      myCase = 3;

      Halfedge firstHedge;
      bool hasFanEdge = normalCoordinates.triangleInequalityViolation(f, firstHedge);
      GC_SAFETY_ASSERT(hasFanEdge, "Triangles in Case 3 must be fan triangles");
      int longHe = heIndex(firstHedge);

      // Fan configuration
      int inputPointIndex = counts[longHe] - 1;

      auto& inputHedgeCurve = curves[heIndex(firstHedge)][inputPointIndex];

      // Identify halfedge in input mesh
      Halfedge inputHalfedge = identifyInputEdge(inputHedgeCurve);
      inputFace = inputHalfedge.face();

      intrinsicInputPairs.push_back(computeIntrinsicAndInputPoints(firstHedge, inputPointIndex + 1));
      intrinsicInputPairs.push_back(computeIntrinsicAndInputPoints(firstHedge.next().next(), 0));

      // Also add in previous crossing
      intrinsicInputPairs.push_back(computeIntrinsicAndInputPoints(firstHedge, inputPointIndex + 2));
    }

    // === Compute input position via least squares in barycentric
    // coordinates

    // Use least squares to express intrinsic barycentric coordinate of
    // inserted vertex as a linear combination of barycentric coordinates of
    // intrinsic points on boundary

    // Let P be the matrix if barycentric coordinates for the boundary
    // points
    // We want to find a vector of coefficients b such that P * b = bary,
    // the barycentric coordinate of the input point
    // This is underdetermined, so we find the minimal-norm solution, i.e.
    // min |b|^2 such that Pb = bary

    // Intrinsic barycentric coordinate matrix
    Eigen::MatrixXd P(3, intrinsicInputPairs.size());
    for (size_t iP = 0; iP < intrinsicInputPairs.size(); iP++) {
      Vector3 intBary = intrinsicInputPairs[iP].first;
      P(0, iP) = intBary.x;
      P(1, iP) = intBary.y;
      P(2, iP) = intBary.z;
    }
    Eigen::VectorXd rhs(3);
    rhs << bary.x, bary.y, bary.z;

    Eigen::MatrixXd PPT = P * P.transpose();

    Eigen::VectorXd lambda = PPT.colPivHouseholderQr().solve(rhs);
    Eigen::VectorXd b = P.transpose() * lambda;

    Eigen::VectorXd err = P * b - rhs;

    insertionBary = Vector3::zero();

    for (size_t iP = 0; iP < intrinsicInputPairs.size(); iP++) {

      // if (!checkAdjacent(intrinsicInputPairs[iP].second,
      //                    SurfacePoint(inputFace, Vector3::zero()))) {
      //    std::cout << "Case: " << myCase << std::endl;
      //     if (myCase == 1) {
      //        std::cout << "\tjustPast: " << isJustPast << std::endl;
      //     }
      // }

      Vector3 inputCrossingBary = intrinsicInputPairs[iP].second.inFace(inputFace).faceCoords;
      insertionBary += b(iP) * inputCrossingBary;
    }


    // Very rarely, due to floating point problems, we get barycentric
    // coordinates with one big value and one small value. In this case,
    // we just set them each to be (1-last value)/2
    double roundingTol = 0.5;
    if (insertionBary.x < -roundingTol || insertionBary.y < -roundingTol || insertionBary.z < -roundingTol ||
        insertionBary.x > 1 + roundingTol || insertionBary.y > 1 + roundingTol || insertionBary.z > 1 + roundingTol) {
      if (insertionBary.x > 1 + roundingTol) {
        if (insertionBary.y < 0) {
          double s = (1. - insertionBary.z) / 2.;
          insertionBary.x = s;
          insertionBary.y = s;
        } else if (insertionBary.z < 0) {
          double s = (1. - insertionBary.y) / 2.;
          insertionBary.x = s;
          insertionBary.z = s;
        }
      } else if (insertionBary.y > 1 + roundingTol) {
        if (insertionBary.x < 0) {
          double s = (1. - insertionBary.z) / 2.;
          insertionBary.x = s;
          insertionBary.y = s;
        } else if (insertionBary.z < 0) {
          double s = (1. - insertionBary.x) / 2.;
          insertionBary.y = s;
          insertionBary.z = s;
        }
      } else if (insertionBary.z > 1 + roundingTol) {
        if (insertionBary.x < 0) {
          double s = (1. - insertionBary.y) / 2.;
          insertionBary.x = s;
          insertionBary.z = s;
        } else if (insertionBary.y < 0) {
          double s = (1. - insertionBary.x) / 2.;
          insertionBary.y = s;
          insertionBary.z = s;
        }
      }
      std::cerr << "insertionBary rounded to " << insertionBary << std::endl;
    }

    // HACK: clamp very small values
    // TODO: edge split?
    insertionBary.x = clamp(insertionBary.x, 0., 1.);
    insertionBary.y = clamp(insertionBary.y, 0., 1.);
    insertionBary.z = clamp(insertionBary.z, 0., 1.);
    insertionBary /= (insertionBary.x + insertionBary.y + insertionBary.z);


    insertionFace = inputFace;
  }

  return {SurfacePoint(insertionFace, insertionBary), counts};
}

Vertex IntegerCoordinatesIntrinsicTriangulation::insertVertex(SurfacePoint pt) {
  Vertex newVertex;
  switch (pt.type) {
  case SurfacePointType::Vertex:
    newVertex = pt.vertex;
    break;
  case SurfacePointType::Edge:
    newVertex = splitEdge(pt.edge, pt.tEdge, false);
    break;
  case SurfacePointType::Face:
    newVertex = splitFace(pt.face, pt.faceCoords, false);
    break;
  }

  return newVertex;
}

Vertex IntegerCoordinatesIntrinsicTriangulation::splitFace(Face f, Vector3 bary, bool verbose) {
  // std::clock_t tStart = std::clock();

  std::array<Vector2, 3> vertCoords = vertexCoordinatesInFace(f);
  Vector2 newPCoord = (bary.x * vertCoords[0] + bary.y * vertCoords[1] + bary.z * vertCoords[2]);

  Vector3 fEdgeLengths{edgeLengths[f.halfedge().next().edge()], edgeLengths[f.halfedge().next().next().edge()],
                       edgeLengths[f.halfedge().edge()]};

  std::array<double, 3> newEdgeLengths;
  newEdgeLengths[0] = displacementLength(bary - Vector3{1, 0, 0}, fEdgeLengths);
  newEdgeLengths[1] = displacementLength(bary - Vector3{0, 1, 0}, fEdgeLengths);
  newEdgeLengths[2] = displacementLength(bary - Vector3{0, 0, 1}, fEdgeLengths);

  SurfacePoint inputMeshPosition;
  std::array<int, 3> counts;
  std::tie(inputMeshPosition, counts) = computeFaceSplitData(f, bary);


  // reorder to fit order of edges incident on newVertex
  std::swap(newEdgeLengths[1], newEdgeLengths[2]);

  auto data = normalCoordinates.computeVertexInsertionData(f, counts);

  Vertex newVertex = intrinsicMesh->insertVertex(f);
  normalCoordinates.applyVertexInsertionData(newVertex, data);

  vertexLocations[newVertex] = inputMeshPosition;

  size_t iE = 0;
  for (Halfedge he : newVertex.outgoingHalfedges()) {
    Edge e = he.edge();
    edgeLengths[e] = newEdgeLengths[iE];

    normalCoordinates.setRoundaboutFromPrevRoundabout(he.twin());
    normalCoordinates.roundabouts[he] = 0;

    iE++;
  }
  normalCoordinates.roundaboutDegrees[newVertex] = 0;

  // Update quantities
  // refreshQuantities();

  for (Face f : newVertex.adjacentFaces()) {
    // depends on edgeLengths, faceAreas
    updateFaceBasis(f);
  }

  for (Vertex v : newVertex.adjacentVertices()) {
    // depends on vertexAngleSums
    updateHalfedgeVectorsInVertex(v);
  }

  vertexAngleSums[newVertex] = 2 * M_PI;
  updateHalfedgeVectorsInVertex(newVertex);

  invokeFaceInsertionCallbacks(f, newVertex);

  return newVertex;
}

Halfedge IntegerCoordinatesIntrinsicTriangulation::splitEdge(Halfedge he, double tSplit) {
  // FIXME need to  return correct halfedge, so can't just call function below
  throw std::runtime_error("not implemented");
  return Halfedge();
}

Vertex IntegerCoordinatesIntrinsicTriangulation::splitEdge(Edge e, double bary, bool verbose) {
  return (e.isBoundary()) ? splitBoundaryEdge(e, bary, verbose) : splitInteriorEdge(e, bary, verbose);
}

Vertex IntegerCoordinatesIntrinsicTriangulation::splitBoundaryEdge(Edge e, double bary, bool verbose) {
  auto inEdge = [](Edge e, SurfacePoint p) -> SurfacePoint {
    switch (p.type) {
    case SurfacePointType::Vertex:
      if (p.vertex == e.halfedge().tailVertex()) {
        return SurfacePoint(e, 0);
      } else if (p.vertex == e.halfedge().tipVertex()) {
        return SurfacePoint(e, 1);
      }
      break;
    case SurfacePointType::Edge:
      if (p.edge == e) {
        return p;
      }
      break;
    default:
      break;
    }
    throw std::runtime_error("SurfacePoint not in edge");
  };

  auto heBary = [&](Halfedge he, double t) -> Vector3 {
    int i = halfedgeIndexInTriangle(he);
    int j = (i + 1) % 3;

    Vector3 bary = Vector3::zero();
    bary[i] = (1 - t);
    bary[j] = (t);

    return bary;
  };
  auto faceEdgeLengths = [&](Face f) -> Vector3 {
    // lengths[i] is the length of the edge opposite the i'th vertex
    return Vector3{edgeLengths[f.halfedge().next().edge()], edgeLengths[f.halfedge().next().next().edge()],
                   edgeLengths[f.halfedge().edge()]};
  };

  if (normalCoordinates[e] >= 0) {
    if (verbose) std::cout << "Easy Edge Split" << std::endl;

    throw std::runtime_error("NSHARP: I think we can get rid of this whole case, since we will "
                             "never have normalCoordinates[e] > 0 at boundary. The one case "
                             "below for the == 0 case should be sufficient. But it's untested "
                             "for == 0, so uncomment this and try it.");

    return Vertex();

    /*
    // Easy case - edge not shared
    // TODO: use normal coordinate edge split code explicitly - it
    // needs fewer geodesic crossing points

    // throw std::runtime_error("Not fully implemented yet - need to mark
    // split segments as fixed if e is");

    Vertex vTipBefore  = e.halfedge().tipVertex();
    Vertex vTailBefore = e.halfedge().tailVertex();
    bool fixedBefore   = isFixed[e];
    isFixed[e]         = false;

    Vertex newVertex =
        splitFace(e.halfedge().face(), heBary(e.halfedge(), bary));
    flipEdgeIfPossible(e);


    // Mark new edges as fixed
    // TODO FIXME this search could potentially fail on a gnarly
    // Delta-complex if the input was a self edge? Fix by passing back the
    // appropriate halfedge from the function below
    Halfedge newHalfedge;
    for (Halfedge he : newVertex.outgoingHalfedges()) {
        if (he.tipVertex() == vTipBefore) {
            newHalfedge = he;
        }
        if (he.tipVertex() == vTipBefore || he.tipVertex() == vTailBefore) {
            if (fixedBefore) {
                isFixed[he.edge()] = true;
                // fixedEdges.push_back(e);
            }
        }
    }


    // TODO: invokeEdgeSplitCallbacks here
    // throw std::runtime_error("Didn't call edge split callbacks");

    invokeEdgeSplitCallbacks(
        e, newHalfedge,
        newHalfedge.next().next().twin().next().next().twin());

    return newVertex;
    */

  } else {
    if (verbose) std::cout << "Shared Edge Split" << std::endl;
    // Hard case - edge also exists in input mesh
    Edge inputEdge;
    SurfacePoint src = vertexLocations[e.halfedge().tailVertex()];
    SurfacePoint dst = vertexLocations[e.halfedge().tipVertex()];

    if (src.type == SurfacePointType::Edge) {
      inputEdge = src.edge;
    } else if (dst.type == SurfacePointType::Edge) {
      inputEdge = dst.edge;
    } else {
      inputEdge = getSharedInputEdge(e.halfedge()).edge();
    }

    double tSrc = inEdge(inputEdge, src).tEdge;
    double tDst = inEdge(inputEdge, dst).tEdge;

    double tInsertion = (1 - bary) * tSrc + bary * tDst;

    SurfacePoint inputPoint(inputEdge, tInsertion);

    std::array<int, 3> newNormalCoordinates = normalCoordinates.computeBoundaryEdgeSplitDataGeodesic(*this, e, bary);

    // Compute new edge lengths
    double oldLen = edgeLengths[e];
    std::array<double, 3> newEdgeLengths{(1 - bary) * oldLen, 0, bary * oldLen};
    newEdgeLengths[1] = displacementLength(heBary(e.halfedge(), bary) - heBary(e.halfedge().next(), 1),
                                           faceEdgeLengths(e.halfedge().face()));

    std::array<bool, 3> newEdgeFixed{isFixed(e), false, isFixed(e)};

    // "edge 2"
    Halfedge newHalfedge = intrinsicMesh->splitEdgeTriangular(e); // TODO: use mutation Manager
    Vertex newVertex = newHalfedge.vertex();
    vertexLocations[newVertex] = inputPoint;

    GC_SAFETY_ASSERT(newHalfedge.isInterior() && !newHalfedge.twin().isInterior(),
                     "I'm wrong about orientation conventions");

    // TODO: write applyVertexInsertionData in NormalCoordinates
    // class

    size_t iE = 0; // TODO: check indexing convention
    // Explicit loop to go counterclockwise
    Halfedge he = newHalfedge;
    while (true) {
      Edge e = he.edge();
      edgeLengths[e] = newEdgeLengths[iE];
      normalCoordinates.edgeCoords[e] = newNormalCoordinates[iE];

      // nsharp: don't need, because we keep this updated with an edge split callback
      // if (!isFixed(e) && newEdgeFixed[iE]) {
      // isFixed[e] = true;
      //// fixedEdges.push_back(e);
      //} else if (isFixed(e) && !newEdgeFixed[iE]) {
      // throw std::runtime_error("Need to remove e from the fixed list");
      //}

      if (!he.isInterior()) break;
      he = he.next().next().twin();
      iE++;
    }

    for (Halfedge he : newVertex.outgoingHalfedges()) {
      // depends on normalCoordinates.edgeCoords
      normalCoordinates.setRoundaboutFromPrevRoundabout(he.twin());
      normalCoordinates.roundabouts[he] = 0;
    }
    normalCoordinates.roundaboutDegrees[newVertex] = 0;

    // Update quantities

    for (Face f : newVertex.adjacentFaces()) {
      // depends on edgeLengths, faceAreas
      updateFaceBasis(f);
    }

    for (Vertex v : newVertex.adjacentVertices()) {
      // depends on vertexAngleSums
      updateHalfedgeVectorsInVertex(v);
    }

    vertexAngleSums[newVertex] = 2 * M_PI;
    updateHalfedgeVectorsInVertex(newVertex);

    triangulationChanged();
    invokeEdgeSplitCallbacks(e, newHalfedge, newHalfedge.next().next().twin().next().next().twin());

    return newVertex;
  }
}

Vertex IntegerCoordinatesIntrinsicTriangulation::splitInteriorEdge(Edge e, double bary, bool verbose) {
  auto inEdge = [](Edge e, SurfacePoint p) -> SurfacePoint {
    switch (p.type) {
    case SurfacePointType::Vertex:
      if (p.vertex == e.halfedge().tailVertex()) {
        return SurfacePoint(e, 0);
      } else if (p.vertex == e.halfedge().tipVertex()) {
        return SurfacePoint(e, 1);
      }
      break;
    case SurfacePointType::Edge:
      if (p.edge == e) {
        return p;
      }
      break;
    default:
      break;
    }
    throw std::runtime_error("SurfacePoint not in edge");
  };

  auto heBary = [&](Halfedge he, double t) -> Vector3 {
    int i = halfedgeIndexInTriangle(he);
    int j = (i + 1) % 3;

    Vector3 bary = Vector3::zero();
    bary[i] = (1 - t);
    bary[j] = (t);

    return bary;
  };
  auto faceEdgeLengths = [&](Face f) -> Vector3 {
    // lengths[i] is the length of the edge opposite the i'th vertex
    return Vector3{edgeLengths[f.halfedge().next().edge()], edgeLengths[f.halfedge().next().next().edge()],
                   edgeLengths[f.halfedge().edge()]};
  };

  if (normalCoordinates[e] >= 0) {
    if (verbose) std::cout << "Easy Edge Split" << std::endl;
    // Easy case - edge not shared
    // TODO: use normal coordinate edge split code explicitly - it
    // needs fewer geodesic crossing points

    Vertex vTipBefore = e.halfedge().tipVertex();
    Vertex vTailBefore = e.halfedge().tailVertex();
    bool fixedBefore = isFixed(e);

    Vertex newVertex = splitFace(e.halfedge().face(), heBary(e.halfedge(), bary));
    bool flipHappened = flipEdgeIfPossible(e);

    if (!flipHappened) throw std::runtime_error("pos-split flip failed!");

    // Find two of the halfedges that together make up the old edge
    Halfedge newHalfedge = e.halfedge().prevOrbitFace().twin();
    Halfedge otherNewHalfedge = newHalfedge.next().next().twin().next().next().twin();


    // Update is-fixed array
    // nsharp: don't need, we keep updated with edge split callback
    // if (fixedBefore) {
    // isFixed[newHalfedge.edge()] = true;
    // isFixed[otherNewHalfedge.edge()] = true;
    //}

    triangulationChanged();
    invokeEdgeSplitCallbacks(e, newHalfedge, otherNewHalfedge);

    return newVertex;

  } else {
    if (verbose) std::cout << "Shared Edge Split" << std::endl;
    // Hard case - edge also exists in input mesh
    Edge inputEdge;
    SurfacePoint src = vertexLocations[e.halfedge().tailVertex()];
    SurfacePoint dst = vertexLocations[e.halfedge().tipVertex()];

    if (src.type == SurfacePointType::Edge) {
      inputEdge = src.edge;
    } else if (dst.type == SurfacePointType::Edge) {
      inputEdge = dst.edge;
    } else {
      inputEdge = getSharedInputEdge(e.halfedge()).edge();
    }

    double tSrc = inEdge(inputEdge, src).tEdge;
    double tDst = inEdge(inputEdge, dst).tEdge;

    double tInsertion = (1 - bary) * tSrc + bary * tDst;

    SurfacePoint inputPoint(inputEdge, tInsertion);

    std::array<int, 4> newNormalCoordinates = normalCoordinates.computeInteriorEdgeSplitDataGeodesic(*this, e, bary);

    // Compute new edge lengths
    double oldLen = edgeLengths[e];
    std::array<double, 4> newEdgeLengths{0, (1 - bary) * oldLen, 0, bary * oldLen};
    newEdgeLengths[0] =
        displacementLength(heBary(e.halfedge().twin(), 1 - bary) - heBary(e.halfedge().twin().next(), 1),
                           faceEdgeLengths(e.halfedge().twin().face()));
    newEdgeLengths[2] = displacementLength(heBary(e.halfedge(), bary) - heBary(e.halfedge().next(), 1),
                                           faceEdgeLengths(e.halfedge().face()));

    std::array<bool, 4> newEdgeFixed{false, isFixed(e), false, isFixed(e)};

    // "edge 2"
    Halfedge newHalfedge = intrinsicMesh->splitEdgeTriangular(e); // TODO: use mutation Manager
    Vertex newVertex = newHalfedge.vertex();
    vertexLocations[newVertex] = inputPoint;

    // TODO: write applyVertexInsertionData in NormalCoordinates
    // class

    size_t iE = 1; // TODO: fix indexing convention
    for (Halfedge he : newVertex.outgoingHalfedges()) {
      Edge e = he.edge();

      edgeLengths[e] = newEdgeLengths[iE];
      normalCoordinates.edgeCoords[e] = newNormalCoordinates[iE];

      // nsharp: don't need this, managed by callback
      // if (!isFixed[e] && newEdgeFixed[iE]) {
      // isFixed[e] = true;
      //// fixedEdges.push_back(e);
      //} else if (isFixed[e] && !newEdgeFixed[iE]) {
      // throw std::runtime_error("Need to remove e from the fixed list");
      //}

      // indexing goes counterclockwise, but loop goes clockwise
      iE = (iE + 3) % 4;
    }

    for (Halfedge he : newVertex.outgoingHalfedges()) {
      // depends on normalCoordinates.edgeCoords
      normalCoordinates.setRoundaboutFromPrevRoundabout(he.twin());
      normalCoordinates.roundabouts[he] = 0;
    }
    normalCoordinates.roundaboutDegrees[newVertex] = 0;

    // Update quantities

    for (Face f : newVertex.adjacentFaces()) {

      // depends on edgeLengths, faceAreas
      updateFaceBasis(f);
    }

    for (Vertex v : newVertex.adjacentVertices()) {
      // depends on vertexAngleSums
      updateHalfedgeVectorsInVertex(v);
    }

    vertexAngleSums[newVertex] = 2 * M_PI;
    updateHalfedgeVectorsInVertex(newVertex);

    triangulationChanged();
    invokeEdgeSplitCallbacks(e, newHalfedge, newHalfedge.next().next().twin().next().next().twin());

    return newVertex;
  }
}

Face IntegerCoordinatesIntrinsicTriangulation::removeInsertedVertex(Vertex v) {
  // Stolen from geometrycentral/signpost_intrinsic_triangulation.cpp
  // Strategy: flip edges until the vertex has degree three, then remove by
  // replacing with a single face
  // TODO needs a proof that this always works... what about self edges, etc?
  // Seems to work well.

  // What about starting with degree < 3? Since this vertex necessarily has
  // angle sum 2PI, this could only happen in the case of degree 2, with
  // exactly degenerate triangles. Since we assume non-degenerate triangles
  // throughout, we'll consider that to not happen.

  if (vertexLocations[v].type == SurfacePointType::Vertex) return Face(); // can't remove original vertices

  if (isOnFixedEdge(v)) {
    return Face(); // don't try to remove boundary vertices, for now at
                   // least
  }

  if (v.isBoundary()) {
    throw std::runtime_error("boundary vertex removal not implemented");
  }

  // === Flip edges until v has degree 3


  size_t iterCount = 0;
  while (v.degree() != 3 && iterCount < 10 * v.degree()) {

    // Find the highest priority edge to flip
    Edge bestFlipEdge;
    double bestFlipScore = -std::numeric_limits<double>::infinity();
    bool bestFlipIsLoop = false;
    for (Edge e : v.adjacentEdges()) {

      double flipScore = checkFlip(e);
      bool isLoop = e.firstVertex() == e.secondVertex();

      // This logic picks the most-preferred edge to flip. The policy is
      // basically "pick the edge whith the highest flipScore", except
      // that we prefer loop edges if there are any.
      if (isLoop) {
        if (bestFlipIsLoop) {
          // if the one we currently have is a loop, only take this
          // one if it is better
          if (flipScore > bestFlipScore) {
            bestFlipScore = flipScore;
            bestFlipEdge = e;
          }

        } else {
          // if the one we currently have is not a loop, always take
          // this one if it is valid
          if (flipScore > 0.) {
            bestFlipScore = flipScore;
            bestFlipEdge = e;
          }
        }

        bestFlipIsLoop = true;
      } else {
        if (!bestFlipIsLoop) { // only overwrite if the best is not a
                               // loop
          if (flipScore > bestFlipScore) {
            bestFlipScore = flipScore;
            bestFlipEdge = e;
            bestFlipIsLoop = false;
          }
        }
      }
    }

    if (bestFlipEdge == Edge()) {
      return Face();
      // throw std::runtime_error("failed to remove vertex " +
      //                             std::to_string(v) +
      //                             ".  Could not find any edge to
      //                             flip");
    }

    // Passing -inf as the tolerance forces us to always do the flip, since
    // we've already verified it above
    // flipEdgeIfPossible(bestFlipEdge,
    //                    -std::numeric_limits<double>::infinity());
    flipEdgeIfPossible(bestFlipEdge);

    iterCount++;
  }

  // give up if something went wrong (eg. flipped edges)
  if (v.degree() != 3) {
    // throw std::runtime_error(
    //     "failed to remove vertex " + std::to_string(v) +
    //     ".  Somehow vertex degree is not 3. Was it 2 beforehand? (which "
    //     "can only be degenerate)");
    return Face();
  }

  // ==== Remove the vertex
  Face newF = intrinsicMesh->removeVertex(v);

  // ==== Update cached data
  // Edge lengths, normal coordinates, and roundabouts should be okay
  updateFaceBasis(newF);
  for (Halfedge he : newF.adjacentHalfedges()) {
    updateHalfedgeVectorsInVertex(he.vertex());
  }

  triangulationChanged();

  return newF;
}


Vertex IntegerCoordinatesIntrinsicTriangulation::moveVertex(Vertex v, Vector2 vec) {

  // Find the insertion location
  TraceOptions options;
  SurfacePoint startP(v);
  TraceGeodesicResult intrinsicTraceResult = traceGeodesic(*this, startP, vec, options);
  SurfacePoint newPositionOnIntrinsic = intrinsicTraceResult.endPoint;
  if (newPositionOnIntrinsic.type == SurfacePointType::Vertex && newPositionOnIntrinsic.vertex == Vertex()) {
    // tracing failed
    return Vertex();
  }

  // Insert a new vertex
  Vertex newV = insertVertex(newPositionOnIntrinsic);

  // Delete the old vertex
  removeInsertedVertex(v);

  triangulationChanged();

  return newV;
}


size_t IntegerCoordinatesIntrinsicTriangulation::nSubdividedVertices() const {
  size_t nNewVertices = 0;
  for (Edge e : intrinsicMesh->edges()) nNewVertices += positivePart(normalCoordinates[e]);

  return intrinsicMesh->nVertices() + nNewVertices;
}

// HACK: represents arcs parallel to a mesh edge with a single pair {-n,
// he} where n is the number of arcs parallel to he.edge() Trace an edge
// of the input mesh over the intrinsic triangulation
NormalCoordinatesCompoundCurve IntegerCoordinatesIntrinsicTriangulation::traceInputEdge(Edge e, bool verbose) const {
  return traceInputHalfedge(e.halfedge(), verbose);
}

NormalCoordinatesCompoundCurve IntegerCoordinatesIntrinsicTriangulation::traceInputHalfedge(Halfedge inputHe,
                                                                                            bool verbose) const {
  // verbose = verbose || e.getIndex() == 18027;

  auto vertexHalfedge = [&](Vertex v, size_t iH) {
    Halfedge he = v.halfedge();

    // Iterate counterclockwise
    for (size_t i = 0; i < iH; ++i) {
      he = he.next().next().twin();
    }
    return he;
  };

  auto traceFollowingComponents = [&](const NormalCoordinatesCurve& firstComponent) -> NormalCoordinatesCompoundCurve {
    NormalCoordinatesCompoundCurve cc;
    bool valid = true;
    NormalCoordinatesCurve nextComponent = firstComponent;
    while (valid) {
      cc.components.push_back(nextComponent);
      std::tie(valid, nextComponent) = traceNextCurve(cc.components.back(), verbose);
    }

    return cc;
  };

  auto directPath = [&](Halfedge he) { return NormalCoordinatesCurve{{std::make_pair(-1, he)}}; };

  Vertex vTrace = inputHe.tailVertex();

  // TODO: Allow different vertex sets
  Vertex v = intrinsicMesh->vertex(vTrace.getIndex());

  Halfedge he = v.halfedge();


  // Loop over all halfedges of intrinsicMesh coming out of v until we
  // find the one whose corner contains the halfedge e.halfedge() of
  // inputMesh
  // TODO: this could be optimized. Calling vertexHalfedge repeatedly
  // repeats a lot of work
  do {
    size_t iStart = normalCoordinates.roundabouts[he];

    size_t em = normalCoordinates.strictDegree(he.corner());
    size_t startInd = -negativePart(normalCoordinates[he.edge()]);
    size_t endInd = -negativePart(normalCoordinates[he.next().next().edge()]);

    size_t width = em + startInd + endInd;

    // Loop over all halfedges of inputMesh coming out of this
    // corner
    for (size_t iH = 0; iH < width; ++iH) {
      Halfedge heTrace = vertexHalfedge(vTrace, iStart + iH);
      if (heTrace != inputHe) continue;

      if (verbose) {
        std::cout << "Tracing from vertex " << v << std::endl;
        std::cout << "Found edge " << inputHe.edge() << " (halfedge " << inputHe << ")" << std::endl;
        std::cout << "\t iH = " << iH << " of " << width << " = " << startInd << " + " << em << " + " << endInd
                  << std::endl;
        std::cout << "\t\t index is " << iStart + iH << " from source vertex " << vTrace << std::endl;
        std::cout << "\t\t past halfedge  " << he << " / edge " << he.edge() << std::endl;
        std::cout << "\t ALL HEDGES: " << std::endl;
        Halfedge he = v.halfedge();

        // Iterate counterclockwise
        for (size_t i = 0; i < v.degree(); ++i) {
          std::cout << "\t\t " << he << "\t (" << he.edge() << ") : " << normalCoordinates.roundabouts[he]
                    << " | edge coord: " << normalCoordinates.edgeCoords[he.edge()]
                    << " | corner deg: " << normalCoordinates.strictDegree(he.corner()) << std::endl;

          he = he.next().next().twin();
        }
      }

      if (iH >= width - endInd) {
        Halfedge pathHe = he.next().next().twin();
        return traceFollowingComponents(directPath(pathHe));
      } else if (iH < startInd) {
        return traceFollowingComponents(directPath(he));
      } else {
        size_t idx = iH - startInd;

        return traceFollowingComponents(
            normalCoordinates.topologicalTrace(he.next(), positivePart(normalCoordinates[he.edge()]) + idx));
      }
    }
    // orbit counterclockwise
    he = he.next().next().twin();
  } while (he != v.halfedge());


  std::cerr << "Something somewhere went horribly wrong" << std::endl;
  he = v.halfedge();
  do {
    size_t iStart = normalCoordinates.roundabouts[he];

    size_t em = positivePart(normalCoordinates[he.next().edge()] - normalCoordinates[he.edge()] -
                             normalCoordinates[he.next().next().edge()]);
    size_t startInd = (normalCoordinates[he.edge()] == 0) ? 1 : 0;
    size_t endInd = (normalCoordinates[he.next().next().edge()] == 0) ? 1 : 0;

    size_t width = em + startInd + endInd;
    std::cerr << "iStart: " << iStart << "\tem: " << em << "\tstartInd: " << startInd << "\t endInd: " << endInd
              << std::endl;
    for (size_t iH = 0; iH < width; ++iH) {
      std::cerr << "\ttrying halfedge " << iStart + iH << " of " << vTrace.degree() << " (aka "
                << normalCoordinates.roundaboutDegrees[v] << ")" << std::endl;
    }
    // orbit counterclockwise
    he = he.next().next().twin();
  } while (he != v.halfedge());

  return {};
}

std::pair<bool, NormalCoordinatesCurve>
IntegerCoordinatesIntrinsicTriangulation::traceNextCurve(const NormalCoordinatesCurve& oldCurve, bool verbose) const {
  GC_SAFETY_ASSERT(!oldCurve.crossings.empty(), "curves must contain some crossings. Even shared edges "
                                                "contain a 'parallel' crossing");
  Halfedge finalHe;
  int finalPos;
  std::tie(finalPos, finalHe) = oldCurve.crossings.back();

  if (finalPos < 0) {

    // Shared edge
    Vertex v = finalHe.tipVertex();
    SurfacePoint inputPos = vertexLocations[v];
    if (inputPos.type == SurfacePointType::Edge) {
      // keep tracing

      // Look for another shared edge
      for (Halfedge he : v.outgoingHalfedges()) {
        if (he != finalHe.twin() && normalCoordinates[he.edge()] < 0) {
          return {true, NormalCoordinatesCurve{{{-1, he}}}};
        }
      }

      // Look for a non-shared curve
      for (Corner c : v.adjacentCorners()) {
        if (normalCoordinates.strictDegree(c) > 0) {
          GC_SAFETY_ASSERT(normalCoordinates.strictDegree(c) == 1, "there can only be one");
          // we found a curve. And it's different than
          // oldCurve since it's not a shared edge
          return {true, normalCoordinates.topologicalTrace(c, 0)};
        }
      }
    } else {
      return {false, NormalCoordinatesCurve{}};
    }
  } else {
    if (verbose) std::cout << "Normal Edge" << std::endl;
    // Ordinary (transverse) edge
    Corner incomingCorner = finalHe.twin().next().next().corner();
    Vertex v = incomingCorner.vertex();
    SurfacePoint inputPos = vertexLocations[v];
    if (inputPos.type == SurfacePointType::Edge) {
      // keep tracing

      // Look for a shared edge
      for (Halfedge he : v.outgoingHalfedges()) {
        if (normalCoordinates[he.edge()] < 0) {
          return {true, NormalCoordinatesCurve{{{-1, he}}}};
        }
      }

      // Look for a non-shared curve
      for (Corner c : v.adjacentCorners()) {
        if (normalCoordinates.strictDegree(c) > 0 && c != incomingCorner) {
          GC_SAFETY_ASSERT(normalCoordinates.strictDegree(c) == 1, "there can only be one");
          // we found a curve. And it's different than
          // oldCurve since it's not a shared edge
          return {true, normalCoordinates.topologicalTrace(c, 0)};
        } else if (normalCoordinates.strictDegree(c) == 2) {
          throw std::runtime_error("Your geometry is really really bad");
          GC_SAFETY_ASSERT(c == incomingCorner, "there can only be one");

          return {true, normalCoordinates.topologicalTrace(c.halfedge().next(), finalPos)};

        } else {
          throw std::runtime_error("this shouldn't be possible");
          return {false, NormalCoordinatesCurve{}};
        }
      }
    } else {
      return {false, NormalCoordinatesCurve{}};
    }
  }
  throw std::runtime_error("This should be unreachable");
}

// Inverse to traceInputEdge
Halfedge IntegerCoordinatesIntrinsicTriangulation::identifyInputEdge(const NormalCoordinatesCurve& path,
                                                                     bool verbose) const {
  GC_SAFETY_ASSERT(path.crossings.size() >= 1, "paths of length 0 do not correspond to edges");

  auto vertexHalfedge = [&](Vertex v, size_t iH) {
    Halfedge he = v.halfedge();

    // Iterate counterclockwise
    for (size_t i = 0; i < iH; ++i) {
      he = he.next().next().twin();
    }
    return he;
  };

  if (path.crossings[0].first < 0) {
    // If edge is shared, use roundabouts to identify shared edge

    Halfedge he = path.crossings[0].second;
    Vertex inputSrc = vertexLocations[he.vertex()].vertex;
    return vertexHalfedge(inputSrc, normalCoordinates.roundabouts[he]);
  } else {
    Halfedge he = path.crossings[0].second;
    int p = path.crossings[0].first - normalCoordinates.strictCornerCoord(he.corner());

    if (verbose) {
      std::cerr << "halfedge: " << he << "\t crossing No: " << path.crossings[0].first
                << "\t cornerCoord: " << normalCoordinates.strictCornerCoord(he.corner()) << "\tp: " << p << std::endl;
    }
    GC_SAFETY_ASSERT(p >= 0, "crossing must have nonnegative index");

    SurfacePoint src = vertexLocations[he.next().next().vertex()];
    GC_SAFETY_ASSERT(src.type == SurfacePointType::Vertex, "edge must start at vertex");
    Vertex inputSrc = src.vertex;

    Halfedge hePrev = he.next().next();
    int rIdx = p + normalCoordinates.roundabouts[hePrev] - negativePart(normalCoordinates[hePrev.edge()]);

    return vertexHalfedge(inputSrc, rIdx);
  }
}

// Identify shared halfedge, throw exception if halfedge is not shared
// (i.e. edgeCoords[he.edge()] must be negative)
Halfedge IntegerCoordinatesIntrinsicTriangulation::identifyInputEdge(Halfedge he) const {
  GC_SAFETY_ASSERT(normalCoordinates[he.edge()] < 0, "shared edge must have edgeCoord -1");

  auto vertexHalfedge = [&](Vertex v, size_t iH) {
    Halfedge he = v.halfedge();

    // Iterate counterclockwise
    for (size_t i = 0; i < iH; ++i) {
      he = he.next().next().twin();
    }
    return he;
  };

  Vertex src = vertexLocations[he.vertex()].vertex;
  return vertexHalfedge(src, normalCoordinates.roundabouts[he]);
}

std::array<Vector2, 3> IntegerCoordinatesIntrinsicTriangulation::vertexCoordinatesInFace(Face face) const {
  // stolen from
  // gc/intrinsic_geometry_interface.cpp:computeHalfedgeVectorsInFace

  // Gather some values
  Halfedge heAB = face.halfedge();
  Halfedge heBC = heAB.next();
  Halfedge heCA = heBC.next();
  GC_SAFETY_ASSERT(heCA.next() == heAB, "faces must be triangular");

  double lAB = edgeLengths[heAB.edge()];
  double lBC = edgeLengths[heBC.edge()];
  double lCA = edgeLengths[heCA.edge()];

  // Assign positions to all three vertices
  Vector2 pA{0., 0.}; // used implicitly
  Vector2 pB{lAB, 0.};
  // pC is the hard one:

  // Herons formula
  // stolen from
  // gc/intrinsic_geometry_interface.cpp:computeFaceAreas
  double s = (lAB + lBC + lCA) / 2.0;
  double arg = s * (s - lAB) * (s - lBC) * (s - lCA);
  arg = std::fmax(0., arg);
  double tArea = std::sqrt(arg);

  // Compute width and height of right triangle formed via altitude
  // from C
  double h = 2. * tArea / lAB;
  double w = std::sqrt(std::max(0., lCA * lCA - h * h));

  // Take the closer of the positive and negative solutions
  if (lBC * lBC > (lAB * lAB + lCA * lCA)) w *= -1.0;

  // Project some vectors to get the actual position
  Vector2 pC{w, h};

  return {pA, pB, pC};
}

// If f is entirely contained in some face of the input mesh, return that
// face Otherwise return Face()
Face IntegerCoordinatesIntrinsicTriangulation::getParentFace(Face f) const {
  auto containsVertex = [](Face f, Vertex v) -> bool {
    for (Vertex vF : f.adjacentVertices()) {
      if (vF == v) return true;
    }
    return false;
  };

  auto containsEdge = [](Face f, Edge e) -> bool {
    for (Edge eF : f.adjacentEdges()) {
      if (eF == e) return true;
    }
    return false;
  };

  auto compatible = [&](const SurfacePoint& pt, Face f) -> bool {
    switch (pt.type) {
    case SurfacePointType::Vertex:
      return containsVertex(f, pt.vertex);
    case SurfacePointType::Edge:
      return containsEdge(f, pt.edge);
    case SurfacePointType::Face:
      return pt.face == f;
    }
    return false; // unreachable
  };

  // Look for a FacePoint
  for (Vertex v : f.adjacentVertices()) {
    SurfacePoint vP = vertexLocations[v];
    if (vP.type == SurfacePointType::Face) {
      Face parentFace = vP.face;

      // Check if this works for everyone else
      for (Vertex w : f.adjacentVertices()) {
        if (!compatible(vertexLocations[w], parentFace)) {
          return Face();
        }
      }

      return parentFace;
    }
  }

  // Look for an EdgePoint
  for (Vertex v : f.adjacentVertices()) {
    if (vertexLocations[v].type == SurfacePointType::Edge) {
      Edge e = vertexLocations[v].edge;
      Face f1 = e.halfedge().face();
      Face f2 = e.halfedge().twin().face();

      bool f1Okay = e.halfedge().isInterior();
      bool f2Okay = e.halfedge().twin().isInterior();

      for (Vertex w : f.adjacentVertices()) {
        f1Okay = f1Okay && compatible(vertexLocations[w], f1);
        f2Okay = f2Okay && compatible(vertexLocations[w], f2);
      }

      return (f1Okay) ? f1 : (f2Okay) ? f2 : Face();
    }
  }

  // Give up
  return Face();
}

// ======================================================
//          Geometry and Helpers
// ======================================================


FaceData<Vector2> interpolateTangentVectorsB(const IntegerCoordinatesIntrinsicTriangulation& tri,
                                             const CommonSubdivision& cs, const FaceData<Vector2>& dataB) {

  FaceData<Vector2> interp(*cs.mesh);

  for (Face f : cs.mesh->faces()) {
    Face fB = cs.sourceFaceB[f];
    if (fB == Face()) {
      std::cerr << "Encountered Face() as parent?" << std::endl;
      continue;
    }

    // Find the position's of fB's vertices in its tangent space
    // Vector2 p0{0, 0};
    Vector2 p1 = tri.halfedgeVectorsInFace[fB.halfedge()];
    Vector2 p2 = -tri.halfedgeVectorsInFace[fB.halfedge().next().next()];

    // Find the position's of f's vertices in barycentric coords on fB
    Vertex v0 = f.halfedge().tailVertex();
    Vector3 v0Bary = cs.sourcePoints[v0]->posB.inFace(fB).faceCoords;
    Vertex v1 = f.halfedge().tipVertex();
    Vector3 v1Bary = cs.sourcePoints[v1]->posB.inFace(fB).faceCoords;

    // Find the position's of f's vertices in fB's tangent space
    Vector2 q0 = v0Bary.y * p1 + v0Bary.z * p2;
    Vector2 q1 = v1Bary.y * p1 + v1Bary.z * p2;

    //
    Vector2 fBasisInFB = (q1 - q0).normalize();

    interp[f] = dataB[fB] / fBasisInFB;
  }

  return interp;
}

} // namespace surface
} // namespace geometrycentral
