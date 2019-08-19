#pragma once

#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/surface_point.h"


namespace geometrycentral {
namespace surface {


struct TraceGeodesicResult {
  SurfacePoint endPoint;                // the point the trace ended at
  std::vector<SurfacePoint> pathPoints; // all points along the path, including start and end
  Vector2 endingDir;                    // the incoming direction to the final point, in its tangent space
  bool hitBoundary = false;             // did the trace stop early because we hit a boundary?
  bool hasPath = false;                 // is pathPoints populated?
};

// These trace routines will always yield a path that looks like:
//   - the start point
//   - 0 or more edge crossings
//   - the end point in a face
// the only exception is tracing on a surface with boundary, which may yield an end point on a edge if the trace hit the
// boundary

// Trace from a surface point, and a vector in the canonical tangent space of that point (represented as a vector in
// that tangent space)
TraceGeodesicResult traceGeodesic(IntrinsicGeometryInterface& geom, SurfacePoint startP, Vector2 traceVec,
                                  bool includePath = false);


// Trace from a point in barycentric coordinates inside some face, where the trace vector is a barycentric displacement
// (which must sum to 0)
TraceGeodesicResult traceGeodesic(IntrinsicGeometryInterface& geom, Face startFace, Vector3 startBary,
                                  Vector3 traceBaryVec, bool includePath = false);


// For a trace which was expected to end very near targetVertex, try to clean up the end of the path to end directly at
// targetVertex
void trimTraceResult(TraceGeodesicResult& traceResult, Vertex targetVertex);

} // namespace surface
} // namespace geometrycentral
