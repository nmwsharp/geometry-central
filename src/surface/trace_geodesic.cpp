#include "geometrycentral/surface/trace_geodesic.h"

#include "geometrycentral/surface/barycentric_coordinate_helpers.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include <Eigen/Dense>

#include <iomanip>

using std::cout;
using std::endl;

namespace geometrycentral {
namespace surface {

// The default trace options
const TraceOptions defaultTraceOptions;

// Helper functions which support tracing
namespace {


// === Geometric subroutines for tracing

// a few parameters
const double TRACE_EPS_TIGHT = 1e-12;
const double TRACE_EPS_LOOSE = 1e-9;
const bool TRACE_PRINT = false;

inline std::array<Vector2, 3> vertexCoordinatesInTriangle(IntrinsicGeometryInterface& geom, Face face) {
  return {Vector2{0., 0.}, geom.halfedgeVectorsInFace[face.halfedge()],
          -geom.halfedgeVectorsInFace[face.halfedge().next().next()]};
}

inline Vector3 faceCoordsToBaryCoords(const std::array<Vector2, 3>& vertCoords, Vector2 faceCoord) {

  // Warning: bakes in assumption that vertCoords[0] == (0,0)

  // Invert the system to solve for coordinates
  double b2 = faceCoord.y / vertCoords[2].y;
  b2 = clamp(b2, 0.0, 1.0); // comment these to get useful errors rather than clamping
  double b1 = (faceCoord.x - b2 * vertCoords[2].x) / vertCoords[1].x;
  b1 = clamp(b1, 0.0, 1.0 - b2);
  double b0 = 1.0 - b1 - b2;
  b0 = clamp(b0, 0.0, 1.0);
  Vector3 result{b0, b1, b2};

  return result;
}

inline Vector2 baryCoordsToFaceCoords(const std::array<Vector2, 3>& vertCoords, Vector3 baryCoord) {
  return vertCoords[0] * baryCoord.x + vertCoords[1] * baryCoord.y + vertCoords[2] * baryCoord.z;
}

inline Vector3 cartesianVectorToBarycentric(const std::array<Vector2, 3>& vertCoords, Vector2 faceVec) {


  // Build matrix for linear transform problem
  // (last constraint comes from chosing the displacement vector with sum = 0)
  Eigen::Matrix3d A;
  Eigen::Vector3d rhs;
  const std::array<Vector2, 3>& c = vertCoords; // short name
  A << c[0].x, c[1].x, c[2].x, c[0].y, c[1].y, c[2].y, 1., 1., 1.;
  rhs << faceVec.x, faceVec.y, 0.;

  // Solve
  Eigen::Vector3d result = A.colPivHouseholderQr().solve(rhs);
  Vector3 resultBary{result(0), result(1), result(2)};

  resultBary = normalizeBarycentricDisplacement(resultBary);

  if (TRACE_PRINT) {
    cout << "       cartesianVectorToBarycentric() " << endl;
    cout << "         input = " << faceVec << endl;
    cout << "         positions = " << vertCoords[0] << " " << vertCoords[1] << " " << vertCoords[2] << endl;
    cout << "         transform result = " << resultBary << endl;
    cout << "         transform back = " << barycentricDisplacementToCartesian(vertCoords, resultBary) << endl;
    cout << " A = " << endl << A << endl;
    cout << " rhs = " << endl << rhs << endl;
    cout << " Ax = " << endl << A * result << endl;
  }

  return resultBary;
}

// Converts tCross from halfedge to edge coordinates, handling sign conventions
inline double convertTToEdge(Halfedge he, double tCross) {
  if (he == he.edge().halfedge()) return tCross;
  return 1.0 - tCross;
}
// Converts vectors in halfedge basis from halfedge to edge coordinates, handling sign conventions
inline Vector2 convertVecToEdge(Halfedge he, Vector2 halfedgeVec) {
  if (he == he.edge().halfedge()) return halfedgeVec;
  return -halfedgeVec;
}


// === Tracing subroutines


// Return type from tracing subroutines
// When trace ends, will set newFace = Face() and newVector = Vector3::zero(). only then is incomingDirToPoint populated
struct TraceSubResult {
  // Did the trace end?
  bool terminated;

  // One of the two sets of values will be defined:

  // If the trace continues (terminated == false)
  Halfedge crossHe;                 // halfedge we crossed over from (crossHe.twin().face() is new face)
  double tCross;                    // t along crossHe
  Vector2 traceVectorInHalfedgeDir; // vector to keep tracing, measured against crossHe
  double traceVectorInHalfedgeLen;  // vector to keep tracing, measured against crossHe

  // If the trace ends (terminated == true)
  SurfacePoint endPoint;      // ending location
  Vector2 incomingDirToPoint; // final incoming direction
};


// General form for tracing barycentrically within a face
// Assumes that approriate projects have already been performed such that startPoint and vectors are valid (inside
// triangle and pointing in the right direction)
//
// This function is tightly coupled with the routines which call it. They prepare the values startPoint, vecBary, and
// vecCartesian, ensuring that those values satisify basic properties (essentially that the trace points in a vallid
// direction).
//
// Note that this expects to be given the trace vector in both barycentric _and_ cartesian coordinates. These are two
// different representations of the same data! This is useful the barycentric representation is good for relilably
// performing tracing, while the cartesian representation is good for transforming the trace vector between triangles.
inline TraceSubResult traceInFaceBarycentric(IntrinsicGeometryInterface& geom, Face face, Vector3 startPoint,
                                             Vector3 vecBary, Vector2 vecCartesianDir, double vecCartesianLen,
                                             std::array<bool, 3> edgeIsHittable, const TraceOptions& traceOptions) {

  // Gather values
  std::array<Vector2, 3> vertexCoords = vertexCoordinatesInTriangle(geom, face);
  Vector3 triangleLengths{geom.edgeLengths[face.halfedge().edge()], geom.edgeLengths[face.halfedge().next().edge()],
                          geom.edgeLengths[face.halfedge().next().next().edge()]};

  if (sum(startPoint) < 0.5) {
    if (TRACE_PRINT) {
      cout << "  bad bary point: " << startPoint << endl;
    }
    if (traceOptions.errorOnProblem) {
      throw std::runtime_error("bad bary point");
    }
  }

  if (TRACE_PRINT) {
    cout << "  general trace in face: " << endl;
    cout << "  face: " << face << " startPoint " << startPoint << " vecBary = " << vecBary << " vecCartesian "
         << vecCartesianDir << " " << vecCartesianLen << endl;
  }

  if (TRACE_PRINT) {
    cout << "  vec bary  = " << vecBary << endl;
    cout << "  reconvert = " << cartesianVectorToBarycentric(vertexCoords, vecCartesianDir) * vecCartesianLen << endl;
  }

  // Test if the vector ends in the triangle
  Vector3 endPoint = startPoint + vecBary;
  if (TRACE_PRINT) {
    cout << "    endpoint: " << endPoint << endl;
  }

  /*
  if (isInsideTriangle(endPoint)) {
    // Fancy test if ending on edges
    // (to detect when we're within eps of edge)
    bool foundEpsEdge = false;
    if (traceOptions.allowEndOnEdge) {
      double A = 0.5 * cross(vertexCoords[1] - vertexCoords[0], vertexCoords[2] - vertexCoords[0]);
      Halfedge endHe = face.halfedge();
      for (int i = 0; i < 3; i++) {
        double dist = 2. * endPoint[i] * A / triangleLengths[i]; // perp distance to edge
        if (endPoint[i] > 0 && dist < traceOptions.allowEndOnEdgeEps) {
          // force the process below to hit this edge, let other logic proceed
          edgeIsHittable[i] = true;
          edgeIsHittable[(i + 1) % 3] = false;
          edgeIsHittable[(i + 2) % 3] = false;
          foundEpsEdge = true;
          break;
        }
      }
    }
    // Simple test if not ending on edges
    if (!foundEpsEdge) {
      // The trace ended! Call it a day.
      TraceSubResult result;
      result.terminated = true;
      result.endPoint = SurfacePoint(face, endPoint);
      result.incomingDirToPoint = vecCartesian;
      return result;
    }
  }
  */

  if (isInsideTriangle(endPoint)) {
    // The trace ended! Call it a day.
    TraceSubResult result;
    result.terminated = true;
    result.endPoint = SurfacePoint(face, endPoint);
    result.incomingDirToPoint = vecCartesianDir;
    result.traceVectorInHalfedgeLen = 0;
    return result;
  }


  // The vector did not end in this triangle. Pick an appropriate point along some edge
  double tRay = std::numeric_limits<double>::infinity();
  Halfedge crossHe = Halfedge();
  int iOppVertEnd = -777;
  Halfedge currHe = face.halfedge();
  for (int i = 0; i < 3; i++) {
    currHe = currHe.next(); // always opposite the i'th vertex

    // Check the crossing
    double tRayThisRaw = -startPoint[i] / vecBary[i];
    double tRayThis = clamp(tRayThisRaw, 0., 1. - TRACE_EPS_LOOSE);

    if (TRACE_PRINT) {
      cout << "    considering intersection:" << endl;
      cout << std::boolalpha;
      cout << "      hittable[(i+1)%3]: " << edgeIsHittable[(i + 1) % 3] << endl;
      cout << "      vecBary[i]: " << vecBary[i] << endl;
      cout << "      startPoint[i]: " << startPoint[i] << endl;
      cout << "      tRayThisRaw: " << tRayThisRaw << endl;
      cout << "      tRayThis: " << tRayThis << endl;
    }


    if (!edgeIsHittable[(i + 1) % 3] || vecBary[i] >= 0) {
      // note should ALWAYS satisfy precondition that vecBary[i] is negative for at least one hittable edge.
      // if not, fix projection of inputs in caller
      continue;
    }

    if (tRayThisRaw < tRay) {
      // This is the new closest intersection
      tRay = tRayThisRaw;
      crossHe = currHe;
      iOppVertEnd = i;
    }
  }

  if (TRACE_PRINT) {
    cout << "    selected intersection:" << endl;
    cout << "      crossHe: " << crossHe << endl;
    cout << "      tRay: " << tRay << endl;
    cout << "      iOppVertEnd: " << iOppVertEnd << endl;
  }

  // Clamp to a sane range
  tRay = clamp(tRay, 0., 1. - TRACE_EPS_LOOSE);


  if (crossHe == Halfedge()) {
    if (traceOptions.errorOnProblem) {
      throw std::logic_error("no halfedge intersection was selected, precondition problem?");
    }
    if (TRACE_PRINT) {
      cout << "    PROBLEM PROBLEM NO INTERSECTION:" << endl;
    }

    // End immediately
    TraceSubResult result;
    result.terminated = true;
    result.endPoint = SurfacePoint(face, startPoint);
    result.incomingDirToPoint = vecCartesianDir;
    result.traceVectorInHalfedgeLen = 0;
    return result;
  }

  // Compute some useful info about the endpoint
  Vector3 endPointOnEdge = startPoint + tRay * vecBary;
  double tCross = endPointOnEdge[(iOppVertEnd + 2) % 3] /
                  (endPointOnEdge[(iOppVertEnd + 1) % 3] + endPointOnEdge[(iOppVertEnd + 2) % 3]);
  if (TRACE_PRINT) {
    cout << "    end point on edge: " << endPointOnEdge << endl;
    cout << "    tCross raw: " << tCross << endl;
  }
  tCross = clamp(tCross, 0., 1.);

  // Rotate the vector in to the frame of crossHe and shorten it
  double lenRemaining = (1.0 - tRay) * vecCartesianLen;
  Vector2 crossingEdgeVec = (vertexCoords[(iOppVertEnd + 2) % 3] - vertexCoords[(iOppVertEnd + 1) % 3]);
  Vector2 remainingDirInHalfedge = vecCartesianDir / crossingEdgeVec.normalize();
  if (!isfinite(remainingDirInHalfedge)) {
    if (TRACE_PRINT) {
      cout << "    NON FINITE REMAINING TRACE" << endl;
      cout << "    lenRemaining = " << lenRemaining << endl;
      cout << "    crossingEdgeVec = " << crossingEdgeVec << endl;
      cout << "    remainingDirInHalfedge = " << remainingDirInHalfedge << endl;
    }

    if (traceOptions.errorOnProblem) {
      throw std::runtime_error("bad value transforming to new edge. is there a zero-length edge?");
    }
  }


  // Stop tracing if we hit a boundary
  if (!crossHe.twin().isInterior() || (traceOptions.barrierEdges && (*traceOptions.barrierEdges)[crossHe.edge()])) {
    // Build the result
    TraceSubResult result;
    result.terminated = true;
    result.endPoint = SurfacePoint(crossHe.edge(), convertTToEdge(crossHe, tCross));
    result.incomingDirToPoint = remainingDirInHalfedge;
    result.traceVectorInHalfedgeLen = lenRemaining;

    return result;
  }

  // Build the result
  TraceSubResult result;
  result.terminated = false;
  result.crossHe = crossHe;
  result.tCross = tCross;
  result.traceVectorInHalfedgeDir = remainingDirInHalfedge;
  result.traceVectorInHalfedgeLen = lenRemaining;
  return result;
}

// Trace within a face towards a given edge. The trace is assumed to start at the vertex opposite towardsHe.
//   - towardsHe: the halfedge we are tracing towards (opposite the source vertex)
//   - vecCartesian: vector to trace, in the cartesian basis of the face
inline TraceSubResult traceInFaceTowardsEdge(IntrinsicGeometryInterface& geom, Halfedge towardsHe,
                                             Vector2 vecCartesianDir, double vecCartesianLen,
                                             const TraceOptions& traceOptions) {

  // Gather some values
  Face face = towardsHe.face();
  Halfedge rootHe = towardsHe.next().next();
  std::array<Vector2, 3> vertexCoords = vertexCoordinatesInTriangle(geom, face);

  if (TRACE_PRINT) {
    cout << "  face trace towards edge " << towardsHe << " vec = " << vecCartesianDir << endl;
    cout << "  wedge vec right = " << geom.halfedgeVectorsInFace[towardsHe.next().next()] << endl;
    cout << "  wedge vec left  = " << -geom.halfedgeVectorsInFace[towardsHe.next()] << endl;
    cout << "  wedge vec opp = " << geom.halfedgeVectorsInFace[towardsHe] << endl;
  }

  // TODO do some reasonable angular projection on the cartesian vector

  // Convert to barycentric
  Vector3 vecBaryCanonical = cartesianVectorToBarycentric(vertexCoords, vecCartesianDir) * vecCartesianLen;
  Vector3 vecBaryFromRoot = permuteBarycentricFromCanonical(vecBaryCanonical, towardsHe.next().next());

  if (TRACE_PRINT) {
    cout << "  canonical bary vec" << vecBaryCanonical << endl;
    cout << "  bary vec before projection " << vecBaryFromRoot << endl;
  }

  { // Project to ensure the vector is inside the triangle
    vecBaryFromRoot.x = std::fmin(vecBaryFromRoot.x, TRACE_EPS_TIGHT);
    vecBaryFromRoot.y = std::fmax(vecBaryFromRoot.y, 0.);
    vecBaryFromRoot.z = std::fmax(vecBaryFromRoot.z, 0.);

    // Manual displacement projection to sum to 0 while perserving above properties
    double diff = -sum(vecBaryFromRoot);
    if (diff > 0) {
      vecBaryFromRoot.y += diff / 2;
      vecBaryFromRoot.z += diff / 2;
    } else {
      vecBaryFromRoot.x += diff;
    }
  }

  if (TRACE_PRINT) {
    cout << "  bary vec after projection " << vecBaryFromRoot << endl;
  }

  // Assemble data to call the general trace function
  int iHe = halfedgeIndexInTriangle(towardsHe.next().next());
  Vector3 startPoint{0., 0., 0.};
  startPoint[iHe] = 1.0;
  Vector3 vecBaryCanonicalFixed = permuteBarycentricToCanonical(vecBaryFromRoot, rootHe);
  std::array<bool, 3> hittable = {{false, false, false}};
  hittable[(iHe + 1) % 3] = true;

  return traceInFaceBarycentric(geom, face, startPoint, vecBaryCanonicalFixed, vecCartesianDir, vecCartesianLen,
                                hittable, traceOptions);
}


// Trace within a face away from a given edge. The trace must hit one of the two opposite edges
//   - fromHe: the halfedge we enter from along the face
//   - tCrossFrom: t value in [0, 1] along fromHe we we enter the face
//   - traceVecInHalfedge: vector to trace, in the basis of fromHe
inline TraceSubResult traceInFaceFromEdge(IntrinsicGeometryInterface& geom, Halfedge fromHe, double tCrossFrom,
                                          Vector2 traceVecInHalfedgeDir, double traceVecInHalfedgeLen,
                                          const TraceOptions& traceOptions) {

  // Possibly terminate at this edge
  /*
  std::cout << "  norm tracevec = " << norm(traceVecInHalfedge) << std::endl;
  if (traceOptions.allowEndOnEdge && norm(traceVecInHalfedge) < traceOptions.allowEndOnEdgeEps) {
    std::cout << "TERMINATING AT EDGE" << std::endl;
    TraceSubResult result;
    result.terminated = true;
    result.endPoint = SurfacePoint(fromHe, tCrossFrom);
    result.incomingDirToPoint = traceVecInHalfedge;
    if (fromHe != fromHe.edge().halfedge()) result.incomingDirToPoint *= -1.;
    return result;
  }
  */

  // Gather some values
  Halfedge faceHe = fromHe.twin(); // the halfedge in hte face we're heading in to
  Face face = faceHe.face();
  std::array<Vector2, 3> vertexCoords = vertexCoordinatesInTriangle(geom, face);

  if (TRACE_PRINT) cout << "  face trace from edge " << fromHe << " vec = " << traceVecInHalfedgeDir << endl;


  // Project the cartesian vector to definitely point in the right direction
  Vector2 traceVecInFaceHalfedgeDir = -traceVecInHalfedgeDir;
  if (TRACE_PRINT) cout << "    vec in face before project " << traceVecInFaceHalfedgeDir << endl;
  traceVecInFaceHalfedgeDir.y = std::fmax(traceVecInFaceHalfedgeDir.y, TRACE_EPS_LOOSE);
  if (TRACE_PRINT) cout << "    vec in face after project " << traceVecInFaceHalfedgeDir << endl;

  // Convert to face coordinates
  Vector2 heDir = geom.halfedgeVectorsInFace[faceHe].normalize();
  Vector2 traceVecInFaceDir = heDir * traceVecInFaceHalfedgeDir;
  if (TRACE_PRINT) cout << "    traceVec in face " << traceVecInFaceDir << endl;

  // Convert to barycentric
  Vector3 vecBaryCanonicalDir = cartesianVectorToBarycentric(vertexCoords, traceVecInFaceDir);
  if (TRACE_PRINT) cout << "    vecBaryCanonical " << vecBaryCanonicalDir << endl;
  Vector3 vecBaryFromEdgeDir = permuteBarycentricFromCanonical(vecBaryCanonicalDir, faceHe);

  if (TRACE_PRINT) cout << "    vec bary before project " << vecBaryFromEdgeDir << endl;
  { // Project to ensure the vector is in the right direction
    vecBaryFromEdgeDir.z = std::fmax(vecBaryFromEdgeDir.z, TRACE_EPS_TIGHT);

    // Manual displacement projection to sum to 0 which perserves above properties
    double diff = -sum(vecBaryFromEdgeDir);
    if (diff > 0) {
      vecBaryFromEdgeDir.z += diff;
    } else {
      vecBaryFromEdgeDir.x += diff / 3.;
      vecBaryFromEdgeDir.y += diff / 3.;
      vecBaryFromEdgeDir.z += diff / 3.;
    }
  }
  if (TRACE_PRINT) cout << "    vec bary after project " << vecBaryFromEdgeDir << endl;

  // Project ensure tCrossFrom is valid
  tCrossFrom = clamp(tCrossFrom, 0., 1.);


  // Assemble data to call the general trace function
  int iHe = halfedgeIndexInTriangle(faceHe);
  Vector3 startPoint{0., 0., 0.};
  startPoint[iHe] = tCrossFrom; // notice: switched from what you'd expect becasue tCrossFrom is defined on twin
  startPoint[(iHe + 1) % 3] = 1.0 - tCrossFrom;
  Vector3 vecBaryCanonicalFixedDir = permuteBarycentricToCanonical(vecBaryFromEdgeDir, faceHe);
  if (TRACE_PRINT) {
    cout << "    iHe = " << iHe << endl;
    cout << "    startPoint = " << startPoint << endl;
    cout << "    canonical bary " << vecBaryCanonicalFixedDir << endl;
  }
  std::array<bool, 3> hittable = {{true, true, true}};
  hittable[iHe] = false;

  return traceInFaceBarycentric(geom, face, startPoint, vecBaryCanonicalFixedDir * traceVecInHalfedgeLen,
                                traceVecInFaceDir, traceVecInHalfedgeLen, hittable, traceOptions);
}


// Trace starting from an edge
inline TraceSubResult traceGeodesic_fromEdge(IntrinsicGeometryInterface& geom, Edge currEdge, double tEdge,
                                             Vector2 currVecDir, double currVecLen, const TraceOptions& traceOptions) {

  if (TRACE_PRINT) cout << "  edge trace " << currEdge << " tEdge = " << tEdge << " edge vec = " << currVecDir << endl;

  // Project to ensure tEdge is valid
  tEdge = clamp(tEdge, 0., 1.);

  // Find coordinates in adjacent face

  // Check which side of the face we're exiting
  Halfedge traceHe;
  Vector2 halfedgeTraceDir;
  if (currVecDir.y >= 0.) {
    traceHe = currEdge.halfedge().twin();
    halfedgeTraceDir = -currVecDir;
    tEdge = 1.0 - tEdge;
  } else {
    traceHe = currEdge.halfedge();

    // Can't go anyywhere if boundary halfedge
    if (!traceHe.twin().isInterior() || (traceOptions.barrierEdges && (*traceOptions.barrierEdges)[traceHe.edge()])) {
      TraceSubResult result;
      result.terminated = true;
      result.endPoint = SurfacePoint(currEdge, tEdge);
      result.incomingDirToPoint = currVecDir;

      return result;
    }

    halfedgeTraceDir = currVecDir;
  }

  return traceInFaceFromEdge(geom, traceHe, tEdge, halfedgeTraceDir, currVecLen, traceOptions);
}

// Trace starting from a face
inline TraceSubResult traceGeodesic_fromFace(IntrinsicGeometryInterface& geom, Face currFace, Vector3 faceBary,
                                             Vector2 currVecDir, double currVecLen, const TraceOptions& traceOptions) {

  // Convert the vector to barycentric
  std::array<Vector2, 3> vertexCoords = vertexCoordinatesInTriangle(geom, currFace);
  Vector3 vecBary = cartesianVectorToBarycentric(vertexCoords, currVecDir) * currVecLen;

  return traceInFaceBarycentric(geom, currFace, faceBary, vecBary, currVecDir, currVecLen, {true, true, true},
                                traceOptions);
}


// Trace starting from a vertex (with a rescaled cartesian vector)
inline TraceSubResult traceGeodesic_fromVertex(IntrinsicGeometryInterface& geom, Vertex currVert, Vector2 currVecDir,
                                               double currVecLen, const TraceOptions& traceOptions) {
  if (TRACE_PRINT) cout << "  vertex trace " << currVert << " edge vec = " << currVecDir << endl;

  // Find the halfedge opening the wedge where tracing will start
  Halfedge wedgeHe;
  Vector2 traceDirRelativeToStart;

  // Normally, one of the interval tests below will return positive and we'll simply launch the trace in to that
  // interval. However, due to numerical misfortune, it is possible that none of the intervals will test positive. In
  // that case, we'll simply launch along whichever halfedge was closest.
  double minCross = std::numeric_limits<double>::infinity();
  Halfedge minCrossHalfedge;
  Vector2 minCrossHalfedgeDir; // the trace vector in this closest halfedge

  Halfedge currHe = currVert.halfedge();
  do {

    // Once we hit the boundary we're done
    // (and traversal below doesn't work on boundary loop, so need to exit specially)
    if (!currHe.isInterior()) {
      break;
    }

    Halfedge nextHe = currHe.next().next().twin();

    // The interval spanned by this edge, which we are currently testing
    Vector2 intervalStart = geom.halfedgeVectorsInVertex[currHe].normalize();
    Vector2 intervalEnd = geom.halfedgeVectorsInVertex[nextHe].normalize();

    if (TRACE_PRINT) {
      cout << "  testing wedge " << intervalStart << " -- " << intervalEnd << endl;
      cout << "    testing wedge (un norm) " << geom.halfedgeVectorsInVertex[currHe] << " -- "
           << geom.halfedgeVectorsInVertex[nextHe] << endl;
      cout << "    corner angle " << geom.cornerAngles[currHe.corner()] << endl;
      cout << "    corner " << currHe.corner() << endl;
      Vector2 relAngle = intervalEnd / intervalStart;
      cout << "    wedge width " << relAngle << " radians: " << relAngle.arg() << endl;
    }


    // Check if our trace vector lies within the interval
    double crossStart = cross(intervalStart, currVecDir);
    double crossEnd = cross(intervalEnd, currVecDir);
    if (crossStart > 0. && crossEnd <= 0.) {
      wedgeHe = currHe;
      traceDirRelativeToStart = currVecDir / intervalStart;
      if (TRACE_PRINT) cout << "    wedge match! relative angle " << traceDirRelativeToStart << endl;
      if (TRACE_PRINT) cout << "    cross start = " << crossStart << " cross end = " << crossEnd << endl;
      break;
    }

    // Keep track of the closest halfedge, as described above
    if (std::fabs(crossStart) < minCross) {
      minCross = std::fabs(crossStart);
      minCrossHalfedge = currHe;
      minCrossHalfedgeDir = Vector2{1, TRACE_EPS_TIGHT};
    }
    if (std::fabs(crossEnd) < minCross) {
      minCross = std::fabs(crossEnd);
      minCrossHalfedge = nextHe;
      minCrossHalfedgeDir = Vector2{1, -TRACE_EPS_TIGHT};
    }

    currHe = nextHe;
  } while (currHe != currVert.halfedge());

  // None of the interval tests passed (probably due to unfortunate numerics), so just trace along the closest
  // halfedge
  if (wedgeHe == Halfedge()) {
    if (TRACE_PRINT) cout << "  no wedge worked. following closest edge with dir " << minCrossHalfedgeDir << endl;
    // Convert to edge coordinates
    currVecDir = convertVecToEdge(minCrossHalfedge, minCrossHalfedgeDir);
    return traceGeodesic_fromEdge(geom, minCrossHalfedge.edge(), convertTToEdge(minCrossHalfedge, 0.), currVecDir,
                                  currVecLen, traceOptions);
  }

  // Compute the actual starting face point, slightly inside and adjacent face
  Face startFace = wedgeHe.face();
  int iHe = halfedgeIndexInTriangle(wedgeHe);

  // Need to convert from "powered" representation to flat vector in face
  double sum = currVert.isBoundary() ? M_PI : 2. * M_PI;
  traceDirRelativeToStart = traceDirRelativeToStart.pow(geom.vertexAngleSums[currVert] / sum);
  traceDirRelativeToStart = traceDirRelativeToStart.normalize();

  // Compute the starting vector
  Vector2 startDirInFace = geom.halfedgeVectorsInFace[wedgeHe].normalize();
  Vector2 traceDirInFace = traceDirRelativeToStart * startDirInFace;
  if (TRACE_PRINT) {
    cout << "  starting vector" << endl;
    cout << "    start wedge vec " << geom.halfedgeVectorsInFace[wedgeHe] << endl;
    cout << "    start wedge vec unit " << geom.halfedgeVectorsInFace[wedgeHe].normalize() << endl;
    cout << "    trace dir in face " << traceDirInFace << endl;
  }


  return traceInFaceTowardsEdge(geom, wedgeHe.next(), traceDirInFace, currVecLen, traceOptions);
}


// Run tracing iteratively in faces, after on of the variants below has gotten it started.
// Will internally add the point path point encoded by prevTraceEnd, don't add beforehand.
void traceGeodesic_iterative(IntrinsicGeometryInterface& geom, TraceGeodesicResult& result, TraceSubResult prevTraceEnd,
                             const TraceOptions& traceOptions) {

  // Now, points are always in faces. Trace until termination.
  size_t iter = 0;
  while (!prevTraceEnd.terminated) {

    // Terminate on iterations
    if (traceOptions.maxIters != INVALID_IND && iter >= traceOptions.maxIters) {

      // Use the last trace as ending data
      result.endPoint = SurfacePoint(prevTraceEnd.crossHe, prevTraceEnd.tCross);
      result.endingDir = prevTraceEnd.traceVectorInHalfedgeDir;
      result.length -= prevTraceEnd.traceVectorInHalfedgeLen;

      if (traceOptions.includePath) {
        result.pathPoints.push_back(result.endPoint);
      }

      return;
    }


    // Construct a point where the previous trace ended
    if (traceOptions.includePath) {
      SurfacePoint currPoint(prevTraceEnd.crossHe.edge(), convertTToEdge(prevTraceEnd.crossHe, prevTraceEnd.tCross));
      result.pathPoints.push_back(currPoint);
    }

    if (TRACE_PRINT) {
      cout << "> tracing from " << prevTraceEnd.crossHe << " t = " << prevTraceEnd.tCross
           << " dir = " << prevTraceEnd.traceVectorInHalfedgeDir << endl;
    }

    // Execute the next step of tracing
    prevTraceEnd =
        traceInFaceFromEdge(geom, prevTraceEnd.crossHe, prevTraceEnd.tCross, prevTraceEnd.traceVectorInHalfedgeDir,
                            prevTraceEnd.traceVectorInHalfedgeLen, traceOptions);
    iter++;
  }

  // Add the final ending point
  if (traceOptions.includePath) {
    result.pathPoints.push_back(prevTraceEnd.endPoint);
  }
  result.endPoint = prevTraceEnd.endPoint;
  result.endingDir = prevTraceEnd.incomingDirToPoint;

  // if we still have remaining length, subtract off from result length
  result.length -= prevTraceEnd.traceVectorInHalfedgeLen;

  // if (std::abs(norm(result.endingDir) - 1.) > .1) throw std::runtime_error("norm problem");

  if (prevTraceEnd.endPoint.type == SurfacePointType::Edge) {
    result.hitBoundary = true;
  }
}


} // namespace

TraceGeodesicResult traceGeodesic(IntrinsicGeometryInterface& geom, SurfacePoint startP, Vector2 traceVec,
                                  const TraceOptions& traceOptions) {
  geom.requireVertexAngleSums();
  geom.requireHalfedgeVectorsInVertex();
  geom.requireHalfedgeVectorsInFace();

  // The output data
  TraceGeodesicResult result;
  result.hasPath = traceOptions.includePath;
  if (traceOptions.includePath) {
    result.pathPoints.push_back(startP);
  }

  if (TRACE_PRINT) cout << "\n>>> Trace query from " << startP << " vec = " << traceVec << endl;

  // Quick out with a zero vector
  if (traceVec.norm2() == 0) {
    geom.unrequireVertexAngleSums();
    geom.unrequireHalfedgeVectorsInVertex();
    geom.unrequireHalfedgeVectorsInFace();

    result.endingDir = Vector2::zero();
    result.length = 0;

    // probably want to ensure we still return a point in a face...
    if (traceOptions.errorOnProblem) {
      throw std::runtime_error("zero vec passed to trace, do something good here");
    }

    return result;
  }

  result.length = norm(traceVec); // store length in result for now

  // Trace the first point, based on what kind of input we got
  TraceSubResult prevTraceEnd;
  switch (startP.type) {
  case SurfacePointType::Vertex: {
    prevTraceEnd = traceGeodesic_fromVertex(geom, startP.vertex, unit(traceVec), norm(traceVec), traceOptions);
    break;
  }
  case SurfacePointType::Edge: {
    prevTraceEnd =
        traceGeodesic_fromEdge(geom, startP.edge, startP.tEdge, unit(traceVec), norm(traceVec), traceOptions);
    break;
  }
  case SurfacePointType::Face: {
    prevTraceEnd =
        traceGeodesic_fromFace(geom, startP.face, startP.faceCoords, unit(traceVec), norm(traceVec), traceOptions);
    break;
  }
  }

  // Keep tracing through triangles until finished
  traceGeodesic_iterative(geom, result, prevTraceEnd, traceOptions);

  geom.unrequireVertexAngleSums();
  geom.unrequireHalfedgeVectorsInVertex();
  geom.unrequireHalfedgeVectorsInFace();

  return result;
}


TraceGeodesicResult traceGeodesic(IntrinsicGeometryInterface& geom, Face startFace, Vector3 startBary,
                                  Vector3 traceBaryVec, const TraceOptions& traceOptions) {


  geom.requireVertexAngleSums();
  geom.requireHalfedgeVectorsInVertex();
  geom.requireHalfedgeVectorsInFace();

  // The output data
  TraceGeodesicResult result;
  result.hasPath = traceOptions.includePath;
  if (traceOptions.includePath) {
    result.pathPoints.push_back(SurfacePoint(startFace, startBary));
  }

  if (TRACE_PRINT) {
    cout << "\n>>> Trace query (barycentric) from " << startFace << " " << startBary << " vec = " << traceBaryVec
         << endl;
  }

  // Early-out if zero
  if (traceBaryVec.norm2() == 0) {
    geom.unrequireVertexAngleSums();
    geom.unrequireHalfedgeVectorsInVertex();
    geom.unrequireHalfedgeVectorsInFace();

    // probably want to ensure we still return a point in a face...
    if (traceOptions.errorOnProblem) {
      throw std::runtime_error("zero vec passed to trace, do something good here");
    }

    result.endingDir = Vector2::zero();
    result.length = 0;
    return result;
  }

  // Make sure the input is sane
  startBary = projectInsideTriangle(startBary);
  traceBaryVec -= Vector3::constant(sum(traceBaryVec) / 3);

  // Construct the cartesian equivalent
  std::array<Vector2, 3> vertexCoords = vertexCoordinatesInTriangle(geom, startFace);
  Vector2 traceVectorCartesian = barycentricDisplacementToCartesian(vertexCoords, traceBaryVec);
  result.length = norm(traceVectorCartesian); // store length in result for now

  // Trace the first point starting inside the face
  TraceSubResult prevTraceEnd =
      traceInFaceBarycentric(geom, startFace, startBary, traceBaryVec, unit(traceVectorCartesian),
                             norm(traceVectorCartesian), {true, true, true}, traceOptions);

  // Keep tracing through triangles until finished
  traceGeodesic_iterative(geom, result, prevTraceEnd, traceOptions);

  geom.unrequireVertexAngleSums();
  geom.unrequireHalfedgeVectorsInVertex();
  geom.unrequireHalfedgeVectorsInFace();

  return result;
}

bool trimTraceResult(TraceGeodesicResult& traceResult, Vertex targetVertex) {

  while (traceResult.pathPoints.size() > 1) {
    SurfacePoint& b = traceResult.pathPoints.back();

    // Remove any edge crossings connected to the target vertex: they're numerical noise because we're already in the
    // 1-ring
    if (b.type == SurfacePointType::Edge &&
        (b.edge.halfedge().vertex() == targetVertex || b.edge.halfedge().twin().vertex() == targetVertex)) {
      traceResult.pathPoints.pop_back();
      traceResult.endingDir = Vector2::undefined();
      continue;
    }

    // Always trim face points
    if (b.type == SurfacePointType::Face) {
      traceResult.pathPoints.pop_back();
      traceResult.endingDir = Vector2::undefined();
      continue;
    }

    // Always trim vertex points
    if (b.type == SurfacePointType::Vertex) {
      traceResult.pathPoints.pop_back();
      traceResult.endingDir = Vector2::undefined();
      continue;
    }

    // we're done here
    break;
  }


  // Check success
  if (traceResult.pathPoints.empty()) return false;

  SurfacePoint& b = traceResult.pathPoints.back();

  switch (b.type) {
  case SurfacePointType::Vertex: {
    if (b.vertex == targetVertex) return true;
    for (Vertex n : b.vertex.adjacentVertices()) {
      if (n == targetVertex) return true;
    }
    break;
  }
  case SurfacePointType::Edge: {
    Halfedge bHe = b.edge.halfedge();
    if (bHe.vertex() == targetVertex) return true;
    if (bHe.twin().vertex() == targetVertex) return true;
    if (bHe.next().next().vertex() == targetVertex) return true;
    if (bHe.twin().next().next().vertex() == targetVertex) return true;
    return false;
    break;
  }
  case SurfacePointType::Face: {
    for (Vertex v : b.face.adjacentVertices()) {
      if (v == targetVertex) return true;
    }
    return false;
    break;
  }
  }

  return false;
}


} // namespace surface
} // namespace geometrycentral
