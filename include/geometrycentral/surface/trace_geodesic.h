#pragma once

#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/surface_point.h"


namespace geometrycentral {
namespace surface {


struct TraceGeodesicResult {
  SurfacePoint endPoint;                // the point the trace ended at
  std::vector<SurfacePoint> pathPoints; // all points along the path, including start and end
  Vector2 endingDir;                    // the incoming direction to the final point, in its tangent space
  bool hitBoundary = false;             // did the trace stop early because we hit a boundary?
  bool hasPath = false;                 // is pathPoints populated?
  double length;                        // length of traced path (generally equals norm of traceVec,
                                        // unless tracing stopped early)
};

struct TraceOptions {
  bool includePath = false;
  bool errorOnProblem = false;
  EdgeData<bool>* barrierEdges = nullptr; // if set, traces will stop when they hit barrier edges
  size_t maxIters = INVALID_IND;
};
extern const TraceOptions defaultTraceOptions;

// These trace routines will always yield a path that looks like:
//   - the start point
//   - 0 or more edge crossings
//   - the end point in a face (unless allowEndOnEdge set)
// the only exception is tracing on a surface with boundary, which may yield an end point on a edge if the trace hit the
// boundary

// Trace from a surface point, and a vector in the canonical tangent space of that point (represented as a vector in
// that tangent space)
TraceGeodesicResult traceGeodesic(IntrinsicGeometryInterface& geom, SurfacePoint startP, Vector2 traceVec,
                                  const TraceOptions& traceOptions = defaultTraceOptions);


// Trace from a point in barycentric coordinates inside some face, where the trace vector is a barycentric displacement
// (which must sum to 0)
TraceGeodesicResult traceGeodesic(IntrinsicGeometryInterface& geom, Face startFace, Vector3 startBary,
                                  Vector3 traceBaryVec, const TraceOptions& traceOptions = defaultTraceOptions);


// For a trace which was expected to end very near targetVertex, try to clean up the end of the path to end directly at
// targetVertex
// TODO currently DOES NOT fix up traceResult.endingDir, so that field is invalid after calling
// Return value indicates success. If true, the resulting ends in the 1-ring of targetVertex as expected.
bool trimTraceResult(TraceGeodesicResult& traceResult, Vertex targetVertex);

// Trace a geodesic segment known to run from vertex startVert to vertex endVert (an intrinsic edge between two
// vertices of the input mesh): traceVec in startVert's tangent space, its norm the segment's length. The result
// runs from startVert through the edge crossings to endVert. The face-by-face trace above is used and trimmed
// with trimTraceResult. It is kept at once when it ends within lengthTol * norm(traceVec) of endVert (bounded from
// the edge lengths) and crossed only edge interiors on the way; otherwise it is kept if the trimmed path's length
// (surfacePathLength) is within lengthTol relative of norm(traceVec). When neither holds (the trace did not reach
// endVert's 1-ring, or turned off course at a degenerate face or at a vertex it passed through), the segment is
// traced again by unfolding its faces along the line, which has neither failure. If that fails too, the
// face-by-face trace is returned as before: trimmed if it trimmed, else untrimmed.
std::vector<SurfacePoint> traceGeodesicBetweenVertices(IntrinsicGeometryInterface& geom, Vertex startVert,
                                                       Vector2 traceVec, Vertex endVert,
                                                       const TraceOptions& traceOptions = defaultTraceOptions,
                                                       double lengthTol = 1e-9);

// Trace by unfolding: the faces the segment crosses are laid out one by one in the segment's frame (startVert at
// the origin, the segment along +x), each crossed edge chosen by the side of the line its new vertex lies on;
// from the vertex positions if geom is a VertexPositionGeometry (the edge lengths lose the width of a sliver).
// Where the line passes within rounding of a vertex both sides are tried. The line arrives where endVert is laid
// out at a point whose chord from the start has the segment's length (within arriveTol relative) and crosses
// every edge crossed so far; the crossings are then put on that chord, so the path ends exactly at endVert.
// traceVec only aims the line: where the line misses endVert (a direction carried through degenerate faces can
// be off while the length is right), it is turned toward endVert's unfolded position and unfolded again (up to
// 8 times). Returns an empty vector if it never arrives.
std::vector<SurfacePoint> traceGeodesicByUnfolding(IntrinsicGeometryInterface& geom, Vertex startVert,
                                                   Vector2 traceVec, Vertex endVert, double arriveTol = 1e-6,
                                                   size_t maxIters = INVALID_IND);

// Length of a path of surface points (from the vertex positions if geom is a VertexPositionGeometry), each step
// measured in a face the two points share; NaN if two consecutive points share no face.
double surfacePathLength(IntrinsicGeometryInterface& geom, const std::vector<SurfacePoint>& path);

} // namespace surface
} // namespace geometrycentral
