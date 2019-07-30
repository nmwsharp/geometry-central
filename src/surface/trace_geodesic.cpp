#include "geometrycentral/surface/trace_geodesic.h"

#include "geometrycentral/surface/vertex_position_geometry.h"

using std::cout;
using std::endl;

namespace geometrycentral {
namespace surface {


// Helper functions which support tracing
namespace {


// === Geometric subroutines for tracing

const double TRACE_EPS_TIGHT = 1e-12;
const double TRACE_EPS_LOOSE = 1e-9;


// Line-line intersection

struct IntersectionResult {
  double tRay;            // where along the input ray the hit happens, in [0,1]
  double tLine;           // where along the line the hit happens, in [0,1]
  double orientationSign; // positive if ray is travelling in correct direction to hit boundary (based on ccw edges and
                          // point in a polygon)
};


inline IntersectionResult rayLineIntersection(Vector2 rayStart, Vector2 ray, Vector2 lineA, Vector2 lineB) {

  Vector2 v1 = rayStart - lineA;
  Vector2 v2 = lineB - lineA;
  Vector2 v3 = ray.rotate90();
  double tRay = cross(v2, v1) / dot(v2, v3);
  double tLine = dot(v1, v3) / dot(v2, v3);
  double orientationSign = -dot(ray, (lineB - lineA).rotate90());

  return IntersectionResult{tRay, tLine, orientationSign};
}


inline std::array<Vector2, 3> vertexCoordinatesInTriangle(IntrinsicGeometryInterface& geom, Face face) {
  return {Vector2{0., 0.}, geom.halfedgeVectorsInFace[face.halfedge()],
          -geom.halfedgeVectorsInFace[face.halfedge().next().next()]};
}

inline int halfedgeIndexInFace(Halfedge he, Face f) {
  if (he == f.halfedge()) return 0;
  if (he == f.halfedge().next()) return 1;
  if (he == f.halfedge().next().next()) return 2;

  return -7777777;
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

  /*
  if (b0 < -TRACE_EPS_LOOSE || b1 < -TRACE_EPS_LOOSE || b2 < -TRACE_EPS_LOOSE || (b0 + b1 + b2) > 1 + TRACE_EPS_LOOSE) {
    cout << result << endl;
    cout << (b0 + b1 + b2) << endl;
    throw std::runtime_error("r2 to bary problem");
  }
  */

  return result;
}

inline Vector2 baryCoordsToFaceCoords(const std::array<Vector2, 3>& vertCoords, Vector3 baryCoord) {
  // Warning: bakes in assumption that vertCoords[0] == (0,0)
  return baryCoord.y * vertCoords[1] + baryCoord.z * vertCoords[2];
}

// Converst tCross from halfedge to edge coordinates, handling sign conventions
inline double convertTToEdge(Halfedge he, double tCross) {
  if (he == he.edge().halfedge()) return tCross;
  return 1.0 - tCross;
}
// Converst vectors in halfedge basis from halfedge to edge coordinates, handling sign conventions
inline Vector2 convertVecToEdge(Halfedge he, Vector2 halfedgeVec) {
  if (he == he.edge().halfedge()) return halfedgeVec;
  return -halfedgeVec;
}

// === Tracing subroutines


inline std::tuple<SurfacePoint, Vector2> traceInFaceBasis(IntrinsicGeometryInterface& geom, Face currFace,
                                                          Vector2 faceCoords, Vector2 currVec,
                                                          const std::array<bool, 3>& edgeIsHittable) {
  // Work in a planar basis for this face
  std::array<Vector2, 3> vertexCoords = vertexCoordinatesInTriangle(geom, currFace);

  // Compute a rough characteristic length. We'll rescale to unit length with this, then transform back to output
  // results. This helps with numerics and saves us from thinking about epsilon sizes.
  double charLength = std::fmax(std::fmax(std::fabs(vertexCoords[1].x), std::fabs(vertexCoords[1].y)),
                                std::fmax(std::fabs(vertexCoords[2].x), std::fabs(vertexCoords[2].y)));
  for (int i = 0; i < 3; i++) {
    vertexCoords[i] /= charLength;
  }
  currVec /= charLength;
  faceCoords /= charLength;


  // Do all three hit tests
  std::array<IntersectionResult, 3> intersections;
  for (int i = 0; i < 3; i++) {
    if (edgeIsHittable[i]) {
      intersections[i] = rayLineIntersection(faceCoords, currVec, vertexCoords[i], vertexCoords[(i + 1) % 3]);
    } else {
      intersections[i] = IntersectionResult{std::numeric_limits<double>::infinity(), -1.0, 0.0};
    }
  }


  // Test to find which intersection happened first
  Halfedge crossingHalfedge;
  double tCross;
  double tRay = std::numeric_limits<double>::infinity();
  Vector2 crossingEdgeVec;
  Halfedge currHe = currFace.halfedge();
  bool haveHit = false;
  for (int i = 0; i < 3; i++) {
    IntersectionResult& isect = intersections[i];

    // Skip not hittable
    if (!edgeIsHittable[i]) {
      currHe = currHe.next();
      continue;
    }

    if (intersections[i].orientationSign > 0 && intersections[i].tRay > -TRACE_EPS_LOOSE &&
        intersections[i].tRay <= tRay) {
      tRay = intersections[i].tRay;
      tCross = intersections[i].tLine;
      crossingHalfedge = currHe;
      crossingEdgeVec = (vertexCoords[(i + 1) % 3] - vertexCoords[i]);
      haveHit = true;
    }

    currHe = currHe.next();
  }

  // If we didn't hit anything, the inputs were bad or an edge case occurred. Stop immediately.
  // Note: this can lead to a repeated point in the path
  if (!haveHit) {
    return std::make_tuple(SurfacePoint(currFace, faceCoordsToBaryCoords(vertexCoords, faceCoords)), Vector2::zero());
  }

  // If the ray would end before exiting the face, end it
  if (tRay >= 1.0) {
    Vector2 endingPos = faceCoords + currVec;
    return std::make_tuple(SurfacePoint(currFace, faceCoordsToBaryCoords(vertexCoords, endingPos)), Vector2::zero());
  }

  // Stay safe kids
  tCross = clamp(tCross, 0.0, 1.0);

  // Quit if we hit a boundary
  if (!crossingHalfedge.twin().isInterior()) {
    return std::make_tuple(SurfacePoint(crossingHalfedge.edge(), convertTToEdge(crossingHalfedge, tCross)),
                           Vector2::zero());
  }

  // Shorten the vector by the amount of distance we traversed in this face
  currVec *= (1.0 - tRay);

  // Transform the trace vector to edge tangent space
  Vector2 halfedgeTraceVec = currVec / crossingEdgeVec.normalize();

  // We know the resulting vector should lie in the -y half halfspace, because we just came from this triangle. Enforce
  // that fact so we can count on it.
  if (halfedgeTraceVec.y > TRACE_EPS_LOOSE) throw std::runtime_error("bad transform over edge");
  halfedgeTraceVec.y = std::min(halfedgeTraceVec.y, -TRACE_EPS_TIGHT);

  // Correct for the canonical orientation of the edge
  currVec = convertVecToEdge(crossingHalfedge, halfedgeTraceVec);

  // Undo scaling for vector
  // (scalings of position don't matter because we return barycentric coords)
  currVec *= charLength;

  // Return a new surface point along the edge
  return std::make_tuple(SurfacePoint(crossingHalfedge.edge(), convertTToEdge(crossingHalfedge, tCross)), currVec);
}


// Trace starting from an edge
inline std::tuple<SurfacePoint, Vector2> traceGeodesic_fromEdge(IntrinsicGeometryInterface& geom, Edge currEdge,
                                                                double tEdge, Vector2 currVec) {

  // Find coordinates in adjacent face

  // Check which side of the face we're exiting
  Halfedge traceHe;
  Vector2 halfedgeTraceVec;
  if (currVec.y >= 0.) {
    traceHe = currEdge.halfedge();
    halfedgeTraceVec = currVec;
  } else {

    traceHe = currEdge.halfedge().twin();
    halfedgeTraceVec = -currVec;
    tEdge = 1.0 - tEdge;

    // Can't go anyywhere if boundary halfedge
    if (!traceHe.isInterior()) {
      return std::make_tuple(SurfacePoint(currEdge, tEdge), Vector2::zero());
    }
  }

  // Find barycentric coordinates in the new face
  Face traceFace = traceHe.face();
  int iHe = halfedgeIndexInFace(traceHe, traceFace);
  std::array<Vector2, 3> vertexCoords = vertexCoordinatesInTriangle(geom, traceFace);
  Vector2 pointPos = (1. - tEdge) * vertexCoords[iHe] + tEdge * vertexCoords[(iHe + 1) % 3];

  // Rotate direction of travel to vertex coords
  Vector2 heDir = geom.halfedgeVectorsInFace[traceHe].normalize();
  Vector2 traceVecInFace = heDir * halfedgeTraceVec;

  // The edge we just came from can't be hit
  std::array<bool, 3> hittable = {{true, true, true}};
  hittable[iHe] = false;

  return traceInFaceBasis(geom, traceFace, pointPos, traceVecInFace, hittable);
}

// Trace starting from a face
inline std::tuple<SurfacePoint, Vector2> traceGeodesic_fromFace(IntrinsicGeometryInterface& geom, Face currFace,
                                                                Vector3 faceBary, Vector2 currVec) {
  Vector2 faceCoords = baryCoordsToFaceCoords(vertexCoordinatesInTriangle(geom, currFace), faceBary);
  return traceInFaceBasis(geom, currFace, faceCoords, currVec, {true, true, true});
}


// Trace starting from a vertex
inline std::tuple<SurfacePoint, Vector2> traceGeodesic_fromVertex(IntrinsicGeometryInterface& geom, Vertex currVert,
                                                                  Vector2 currVec) {
  double traceLen = currVec.norm();

  // Find the halfedge opening the wedge where tracing will start

  Halfedge wedgeHe;
  Vector2 traceVecRelativeToStart;


  // Normally, one of the interval tests below will return positive and we'll simply launch the trace in to that
  // interval. However, due to numerical misfortune, it is possible that none of the intervals will test positive. In
  // that case, we'll simply launch along whichever halfedge was closest.
  double minCross = std::numeric_limits<double>::infinity();
  Halfedge minCrossHalfedge;
  Vector2 minCrossHalfedgeVec; // the trace vector in this closest halfedge

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

    // Check if our trace vector lies within the interval
    double crossStart = cross(intervalStart, currVec);
    double crossEnd = cross(intervalEnd, currVec);
    if (crossStart > 0. && crossEnd <= 0.) {
      wedgeHe = currHe;
      traceVecRelativeToStart = currVec / intervalStart;
      break;
    }

    // Keep track of the closest halfedge, as described above
    if (std::fabs(crossStart) < minCross) {
      minCross = std::fabs(crossStart);
      minCrossHalfedge = currHe;
      minCrossHalfedgeVec = Vector2{1, TRACE_EPS_TIGHT} * traceLen;
    }
    if (std::fabs(crossEnd) < minCross) {
      minCross = std::fabs(crossEnd);
      minCrossHalfedge = nextHe;
      minCrossHalfedgeVec = Vector2{1, -TRACE_EPS_TIGHT} * traceLen;
    }

    currHe = nextHe;
  } while (currHe != currVert.halfedge());

  // None of the interval tests passed (probably due to unfortunate numerics), so just trace along the closest halfedge
  if (wedgeHe == Halfedge()) {
    // Convert to edge coordinates
    currVec = convertVecToEdge(minCrossHalfedge, minCrossHalfedgeVec);
    return traceGeodesic_fromEdge(geom, minCrossHalfedge.edge(),
                                  convertTToEdge(minCrossHalfedge, 1.0 - TRACE_EPS_LOOSE), currVec);
  }

  // Compute the actual starting face point, slightly inside and adjacent face
  Face startFace = wedgeHe.face();
  int iHe = halfedgeIndexInFace(wedgeHe, startFace);

  // Compute the starting vector
  Vector2 startDirInFace = geom.halfedgeVectorsInFace[wedgeHe].normalize();
  Vector2 traceVecInFace = traceVecRelativeToStart * startDirInFace;

  // Gather data to trace in this face
  Vector2 pointInFace = vertexCoordinatesInTriangle(geom, startFace)[iHe];
  std::array<bool, 3> hittable = {{false, false, false}};
  hittable[(iHe + 1) % 3] = true;

  return traceInFaceBasis(geom, startFace, pointInFace, traceVecInFace, hittable);
}

} // namespace


TraceGeodesicResult traceGeodesic(IntrinsicGeometryInterface& geom, SurfacePoint startP, Vector2 traceVec,
                                  bool includePath) {

  geom.requireHalfedgeVectorsInVertex();
  geom.requireHalfedgeVectorsInFace();

  SurfacePoint currPoint = startP;
  Vector2 currVec = traceVec;

  // The output data
  TraceGeodesicResult result;
  if (includePath) {
    result.pathPoints.push_back(currPoint);
  }

  while (currVec.norm2() > 0) {

    SurfacePoint newPoint;
    Vector2 newVec;

    // Dispatch to the appropriate trace subroutine
    switch (currPoint.type) {
    case SurfacePointType::Vertex: {
      std::tie(newPoint, newVec) = traceGeodesic_fromVertex(geom, currPoint.vertex, currVec);
      break;
    }
    case SurfacePointType::Edge: {
      std::tie(newPoint, newVec) = traceGeodesic_fromEdge(geom, currPoint.edge, currPoint.tEdge, currVec);
      break;
    }
    case SurfacePointType::Face: {
      std::tie(newPoint, newVec) = traceGeodesic_fromFace(geom, currPoint.face, currPoint.faceCoords, currVec);
      break;
    }
    }

    // Update the current location and vector
    currPoint = newPoint;
    currVec = newVec;

    // Include this latest point along the path
    if (includePath) {
      result.pathPoints.push_back(currPoint);
    }
  }

  result.endPoint = currPoint;

  geom.unrequireHalfedgeVectorsInVertex();
  geom.unrequireHalfedgeVectorsInFace();

  return result;
}


} // namespace surface
} // namespace geometrycentral
