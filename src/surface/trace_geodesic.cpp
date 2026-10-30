#include "geometrycentral/surface/trace_geodesic.h"

#include "geometrycentral/surface/barycentric_coordinate_helpers.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include <Eigen/Dense>

#include <iomanip>
#include <array>
#include <cmath>
#include <limits>

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


namespace {
// Position of a surface point in face f's layout (f.halfedge()'s tail at the origin, as
// vertexCoordinatesInTriangle); false if the point is not on f.
bool pointInFace(IntrinsicGeometryInterface& geom, Face f, const SurfacePoint& p, Vector2& pos) {
  std::array<Vector2, 3> c = vertexCoordinatesInTriangle(geom, f);
  Halfedge h = f.halfedge();
  switch (p.type) {
  case SurfacePointType::Vertex:
    for (int k = 0; k < 3; k++, h = h.next())
      if (h.tailVertex() == p.vertex) {
        pos = c[k];
        return true;
      }
    return false;
  case SurfacePointType::Edge:
    for (int k = 0; k < 3; k++, h = h.next())
      if (h.edge() == p.edge) {
        Vector2 tail = c[k], tip = c[(k + 1) % 3];
        if (h != p.edge.halfedge()) std::swap(tail, tip);
        pos = tail + p.tEdge * (tip - tail);
        return true;
      }
    return false;
  case SurfacePointType::Face:
    if (p.face != f) return false;
    pos = c[0] * p.faceCoords.x + c[1] * p.faceCoords.y + c[2] * p.faceCoords.z;
    return true;
  }
  return false;
}
} // namespace

namespace {
// Whether surface point p lies on face f (vertex of f, edge of f, or in f): pointInFace without the layout.
bool onFace(const SurfacePoint& p, Face f) {
  Halfedge h = f.halfedge();
  switch (p.type) {
  case SurfacePointType::Vertex:
    for (int k = 0; k < 3; k++, h = h.next())
      if (h.tailVertex() == p.vertex) return true;
    return false;
  case SurfacePointType::Edge:
    for (int k = 0; k < 3; k++, h = h.next())
      if (h.edge() == p.edge) return true;
    return false;
  case SurfacePointType::Face:
    return p.face == f;
  }
  return false;
}

// Whether two consecutive path points share a face, over the faces of a in the order surfacePathLength tries them.
bool shareFace(const SurfacePoint& a, const SurfacePoint& b) {
  switch (a.type) {
  case SurfacePointType::Vertex:
    for (Face f : a.vertex.adjacentFaces())
      if (onFace(b, f)) return true;
    return false;
  case SurfacePointType::Edge:
    if (a.edge.halfedge().isInterior() && onFace(b, a.edge.halfedge().face())) return true;
    if (a.edge.halfedge().twin().isInterior() && onFace(b, a.edge.halfedge().twin().face())) return true;
    return false;
  case SurfacePointType::Face:
    return onFace(b, a.face);
  }
  return false;
}
} // namespace

double surfacePathLength(IntrinsicGeometryInterface& geom, const std::vector<SurfacePoint>& path) {
  VertexPositionGeometry* positions = dynamic_cast<VertexPositionGeometry*>(&geom);
  if (positions) {
    // from the positions (the edge lengths lose the width of a sliver); a step's face is only checked, so no layout
    const VertexData<Vector3>& pos = positions->inputVertexPositions;
    double L = 0.;
    Vector3 prev = path.empty() ? Vector3::zero() : path[0].interpolate(pos);
    for (size_t i = 0; i + 1 < path.size(); i++) {
      if (!shareFace(path[i], path[i + 1])) return std::numeric_limits<double>::quiet_NaN();
      const Vector3 next = path[i + 1].interpolate(pos);
      L += norm(next - prev);
      prev = next;
    }
    return L;
  }
  geom.requireHalfedgeVectorsInFace();
  double L = 0.;
  std::vector<Face> faces;
  for (size_t i = 0; i + 1 < path.size(); i++) {
    const SurfacePoint& p = path[i];
    faces.clear();
    switch (p.type) {
    case SurfacePointType::Vertex:
      for (Face f : p.vertex.adjacentFaces()) faces.push_back(f);
      break;
    case SurfacePointType::Edge:
      if (p.edge.halfedge().isInterior()) faces.push_back(p.edge.halfedge().face());
      if (p.edge.halfedge().twin().isInterior()) faces.push_back(p.edge.halfedge().twin().face());
      break;
    case SurfacePointType::Face:
      faces.push_back(p.face);
      break;
    }
    bool shared = false;
    for (Face f : faces) {
      Vector2 a, b;
      if (pointInFace(geom, f, p, a) && pointInFace(geom, f, path[i + 1], b)) {
        L += norm(b - a);
        shared = true;
        break;
      }
    }
    if (!shared) {
      geom.unrequireHalfedgeVectorsInFace();
      return std::numeric_limits<double>::quiet_NaN();
    }
  }
  geom.unrequireHalfedgeVectorsInFace();
  return L;
}

namespace {
// One unfolding of the line from startVert at tangent-space angle theta (see traceGeodesicByUnfolding).
// Returns true with the path in out on arrival; otherwise aimOffset is the angle (in the start's corner,
// radians) of endVert's unfolded position as seen from the start, where it was laid out nearest the line
// (NaN if it never was).
bool unfoldLine(IntrinsicGeometryInterface& geom, VertexPositionGeometry* positions, Vertex startVert,
                double theta, double length, Vertex endVert, double arriveTol, size_t maxIters,
                std::vector<SurfacePoint>& out, double& aimOffset) {
  out.clear();
  aimOffset = std::numeric_limits<double>::quiet_NaN();
  auto arg = [](Vector2 v) {
    double t = std::atan2(v.y, v.x);
    return t < 0. ? t + 2. * M_PI : t;
  };

  // The corner of startVert holding the direction, and the direction's angle from that corner's first side
  const double sumScale = geom.vertexAngleSums[startVert] / (startVert.isBoundary() ? M_PI : 2. * M_PI);
  Halfedge he;
  double d0 = 0.;
  {
    double bestOut = std::numeric_limits<double>::infinity();
    Halfedge curr = startVert.halfedge();
    do {
      if (!curr.isInterior()) break;
      const double width = geom.cornerAngles[curr.corner()] / sumScale; // in the tangent space's angles
      double off = theta - arg(geom.halfedgeVectorsInVertex[curr]);
      off -= 2. * M_PI * std::floor(off / (2. * M_PI));
      if (off <= width) {
        he = curr;
        d0 = off * sumScale;
        break;
      }
      // closest corner, in case rounding leaves the direction in none
      const double outBy = std::fmin(off - width, 2. * M_PI - off);
      if (outBy < bestOut) {
        bestOut = outBy;
        he = curr;
        d0 = (2. * M_PI - off < off - width) ? 0. : geom.cornerAngles[curr.corner()];
      }
      curr = curr.next().next().twin();
    } while (curr != startVert.halfedge());
  }
  if (he == Halfedge()) return false;
  double c0 = geom.cornerAngles[he.corner()];
  if (positions) { // the corner from the positions: from the lengths it is off by about 1e-8 near 0 and pi
    const Vector3 p0 = positions->inputVertexPositions[startVert];
    const Vector3 e1 = positions->inputVertexPositions[he.tipVertex()] - p0;
    const Vector3 e2 = positions->inputVertexPositions[he.next().tipVertex()] - p0;
    if (norm(e1) > 0. && norm(e2) > 0.) c0 = std::atan2(norm(cross(e1, e2)), dot(e1, e2));
  }
  d0 = clamp(d0, 0., c0);

  // The current face: the line enters it at the corner T (first face: the start) or through an edge, and leaves
  // it through heCross (u -> w, u right of the line, w left). Positions in the segment's frame.
  struct State {
    Halfedge heCross;
    Vertex u, w, T;
    Vector2 pu, pw, pT;
    size_t nOut, steps;
  };
  State st;
  st.heCross = he.next();
  st.u = he.tipVertex();
  st.w = he.next().tipVertex();
  st.T = startVert;
  const double lu = geom.edgeLengths[he.edge()], lw = geom.edgeLengths[he.next().next().edge()];
  st.pu = Vector2{lu * std::cos(d0), -lu * std::sin(d0)};
  st.pw = Vector2{lw * std::cos(c0 - d0), lw * std::sin(c0 - d0)};
  st.pT = Vector2{0., 0.};
  st.steps = 0;
  out.push_back(SurfacePoint(startVert));
  st.nOut = out.size();

  // Where the line passes within rounding of a vertex, the side taken; the other one is tried on failure
  struct Decision {
    State before;
    Halfedge heUR, heRW;
    Vertex r;
    Vector2 pr;
  };
  std::vector<Decision> decisions;
  const double hitTol = 1e-9 * length;
  size_t attempts = 0;
  const size_t maxAttempts = 64;
  double nearest = std::numeric_limits<double>::infinity();

  // the edge each crossing lies on, in the segment's frame (out[k + 1] on sides[k]), and whether the crossed
  // halfedge is its edge's halfedge (the crossing's t from u)
  std::vector<std::pair<Vector2, Vector2>> sides;
  std::vector<char> sideForward;

  // endVert laid out at pB: arrived if the chord from the start to pB has the segment's length (within arriveTol)
  // and crosses every edge crossed so far between its ends, i.e. is a straight line from startVert to endVert
  // in these faces. The crossings are moved onto the chord, so the path ends exactly at endVert. Otherwise
  // pB, if nearest the line so far, aims the next unfolding.
  auto tryArrive = [&](Vector2 pB) {
    const double lB = norm(pB);
    if (std::fabs(lB - length) > arriveTol * length) return false;
    std::vector<double> ts(sides.size());
    for (size_t k = 0; k < sides.size(); k++) {
      const Vector2 pu = sides[k].first, pw = sides[k].second, d = pw - pu;
      const double den = cross(pB, d);
      if (den == 0.) return false;
      const double t = -cross(pB, pu) / den;   // along u -> w
      const double lam = cross(pu, d) / den;   // along the chord
      const double eps = 1e-9;
      if (!(t >= -eps && t <= 1. + eps && lam >= -eps && lam <= 1. + eps)) return false;
      ts[k] = clamp(t, 0., 1.);
    }
    for (size_t k = 0; k < sides.size(); k++) out[k + 1].tEdge = sideForward[k] ? ts[k] : 1. - ts[k];
    out.push_back(SurfacePoint(endVert));
    return true;
  };
  auto noteEnd = [&](Vector2 pB) {
    if (std::fabs(pB.y) < nearest) {
      nearest = std::fabs(pB.y);
      aimOffset = std::atan2(pB.y, pB.x);
    }
  };

  // the line takes the side of r given by left (exit through u -> r) or right (through r -> w)
  auto take = [&](State& s, Halfedge heUR, Halfedge heRW, Vertex r, Vector2 pr, bool left) {
    if (left) {
      s.T = s.w;
      s.pT = s.pw;
      s.heCross = heUR;
      s.w = r;
      s.pw = pr;
    } else {
      s.T = s.u;
      s.pT = s.pu;
      s.heCross = heRW;
      s.u = r;
      s.pu = pr;
    }
  };
  // on failure: back to the last open decision, other side; false if none is left
  auto backtrack = [&]() {
    if (decisions.empty() || ++attempts > maxAttempts) return false;
    Decision d = decisions.back();
    decisions.pop_back();
    st = d.before;
    out.resize(st.nOut);
    sides.resize(st.nOut - 1);
    sideForward.resize(st.nOut - 1);
    take(st, d.heUR, d.heRW, d.r, d.pr, !(d.pr.y > 0.));
    return true;
  };

  while (true) {
    // where the line leaves the current face
    const double den = st.pu.y - st.pw.y;
    const double s = clamp(den != 0. ? st.pu.y / den : 0.5, 0., 1.);
    const double xc = st.pu.x + s * (st.pw.x - st.pu.x);
    bool failed = false;
    if (xc >= length * (1. - arriveTol)) {
      // the line ends in this face (or on its exit edge): endVert should be a corner here
      const Vertex vs[3] = {st.u, st.w, st.T};
      const Vector2 ps[3] = {st.pu, st.pw, st.pT};
      for (int k = 0; k < 3; k++)
        if (vs[k] == endVert) {
          if (tryArrive(ps[k])) return true;
          noteEnd(ps[k]);
        }
      failed = xc > length * (1. + arriveTol);
    }
    Halfedge heE = st.heCross.twin();
    if (!failed && (!heE.isInterior() || ++st.steps > maxIters)) failed = true;
    if (failed) {
      if (!backtrack()) {
        out.clear();
        return false;
      }
      continue;
    }
    const bool fwd = st.heCross == st.heCross.edge().halfedge();
    out.push_back(SurfacePoint(st.heCross.edge(), fwd ? s : 1. - s));
    sides.emplace_back(st.pu, st.pw);
    sideForward.push_back(fwd);

    // the next face: r beyond u -> w, from the vertex positions if the geometry has them (the edge lengths
    // lose the width of a sliver), else from the edge lengths
    Halfedge heUR = heE.next(), heRW = heUR.next();
    Vertex r = heRW.tailVertex();
    Vector2 e = st.pw - st.pu;
    const double L = norm(e);
    const double a = geom.edgeLengths[heUR.edge()], b = geom.edgeLengths[heRW.edge()];
    Vector2 pr;
    if (L > 0.) {
      const Vector2 eh = e / L, nh{eh.y, -eh.x}; // nh: right of u -> w, away from the current face
      double along = (a * a - b * b + L * L) / (2. * L);
      double h = std::sqrt(std::fmax(0., a * a - along * along));
      if (positions) {
        const Vector3 pU = positions->inputVertexPositions[st.u];
        const Vector3 uw = positions->inputVertexPositions[st.w] - pU, ur = positions->inputVertexPositions[r] - pU;
        const double l3 = norm(uw);
        if (l3 > 0.) {
          along = dot(ur, uw) / l3;
          h = norm(cross(uw, ur)) / l3;
        }
      }
      pr = st.pu + along * eh + h * nh;
    } else {
      pr = st.pu + Vector2{a, 0.};
    }
    if (r == endVert) {
      if (tryArrive(pr)) return true;
      noteEnd(pr);
    }
    const bool left = pr.y > 0.;
    if (std::fabs(pr.y) <= hitTol) {
      State before = st;
      before.nOut = out.size(); // the crossing of heCross stays
      decisions.push_back(Decision{before, heUR, heRW, r, pr});
    }
    take(st, heUR, heRW, r, pr, left);
    st.nOut = out.size();
  }
}
} // namespace

std::vector<SurfacePoint> traceGeodesicByUnfolding(IntrinsicGeometryInterface& geom, Vertex startVert,
                                                   Vector2 traceVec, Vertex endVert, double arriveTol,
                                                   size_t maxIters) {
  std::vector<SurfacePoint> out;
  const double length = norm(traceVec);
  if (!(length > 0) || !std::isfinite(length) || startVert == endVert) return out;
  geom.requireEdgeLengths();
  geom.requireCornerAngles();
  geom.requireVertexAngleSums();
  geom.requireHalfedgeVectorsInVertex();
  if (maxIters == INVALID_IND) maxIters = 10 * startVert.getMesh()->nFaces();
  VertexPositionGeometry* positions = dynamic_cast<VertexPositionGeometry*>(&geom);

  // The direction aims the line; where it misses endVert (a direction carried through degenerate
  // faces can be off while the length is right), it is turned toward endVert's unfolded position and
  // the line unfolded again
  const double sumScale = geom.vertexAngleSums[startVert] / (startVert.isBoundary() ? M_PI : 2. * M_PI);
  double theta = std::atan2(traceVec.y, traceVec.x);
  const int maxShots = 8;
  bool ok = false;
  for (int shot = 0; shot < maxShots && !ok; shot++) {
    double aim;
    ok = unfoldLine(geom, positions, startVert, theta, length, endVert, arriveTol, maxIters, out, aim);
    if (!ok) {
      if (!std::isfinite(aim)) break;
      theta += aim / sumScale;
    }
  }
  if (!ok) out.clear();

  geom.unrequireEdgeLengths();
  geom.unrequireCornerAngles();
  geom.unrequireVertexAngleSums();
  geom.unrequireHalfedgeVectorsInVertex();
  return out;
}

namespace {
// An upper bound on the distance from surface point p to vertex v, from the edge lengths: exact on an edge of v,
// by the triangle inequality over the corners in a face of v; infinite if p is not on an edge or face of v.
double endDistanceBound(IntrinsicGeometryInterface& geom, const SurfacePoint& p, Vertex v) {
  const double inf = std::numeric_limits<double>::infinity();
  switch (p.type) {
  case SurfacePointType::Vertex:
    return p.vertex == v ? 0. : inf;
  case SurfacePointType::Edge: {
    Halfedge h = p.edge.halfedge();
    geom.requireEdgeLengths();
    const double l = geom.edgeLengths[p.edge];
    geom.unrequireEdgeLengths();
    if (h.tailVertex() == v) return std::fabs(p.tEdge) * l;
    if (h.tipVertex() == v) return std::fabs(1. - p.tEdge) * l;
    return inf;
  }
  case SurfacePointType::Face: {
    Halfedge h = p.face.halfedge();
    const double b[3] = {p.faceCoords.x, p.faceCoords.y, p.faceCoords.z}; // at h's tail, next's tail, next's tip
    geom.requireEdgeLengths();
    double d = inf;
    for (int k = 0; k < 3; k++, h = h.next()) {
      if (h.tailVertex() != v) continue;
      // v at corner k: the other two corners along h (corner k + 1) and along the previous halfedge (k + 2)
      d = std::fabs(b[(k + 1) % 3]) * geom.edgeLengths[h.edge()] +
          std::fabs(b[(k + 2) % 3]) * geom.edgeLengths[h.next().next().edge()];
    }
    geom.unrequireEdgeLengths();
    return d;
  }
  }
  return inf;
}
} // namespace

std::vector<SurfacePoint> traceGeodesicBetweenVertices(IntrinsicGeometryInterface& geom, Vertex startVert,
                                                       Vector2 traceVec, Vertex endVert,
                                                       const TraceOptions& traceOptions, double lengthTol) {
  TraceOptions options = traceOptions;
  options.includePath = true;
  TraceGeodesicResult result = traceGeodesic(geom, SurfacePoint(startVert), traceVec, options);
  const double len = norm(traceVec);
  // A trace that ends at endVert and crossed only edge interiors on the way went right: then the trimmed trace is
  // kept without measuring it. At a vertex (or an edge's end) the trace picks a new direction, and can still end at
  // endVert by another route (near-coincident vertices), so such a trace is measured.
  bool atEnd = !result.pathPoints.empty() && endDistanceBound(geom, result.pathPoints.back(), endVert) <= lengthTol * len;
  for (size_t i = 1; atEnd && i + 1 < result.pathPoints.size(); i++) {
    const SurfacePoint& p = result.pathPoints[i];
    if (p.type != SurfacePointType::Edge || !(p.tEdge > 0. && p.tEdge < 1.)) atEnd = false;
  }
  const bool trimmed = trimTraceResult(result, endVert);
  if (trimmed) {
    result.pathPoints.push_back(SurfacePoint(endVert));
    if (atEnd) return std::move(result.pathPoints);
    if (std::fabs(surfacePathLength(geom, result.pathPoints) - len) <= lengthTol * len) return std::move(result.pathPoints);
  }
  std::vector<SurfacePoint> unfolded = traceGeodesicByUnfolding(geom, startVert, traceVec, endVert, 1e-6, options.maxIters);
  if (!unfolded.empty()) return unfolded;
  if (trimmed) return std::move(result.pathPoints);
  return traceGeodesic(geom, SurfacePoint(startVert), traceVec, options).pathPoints;
}


} // namespace surface
} // namespace geometrycentral
