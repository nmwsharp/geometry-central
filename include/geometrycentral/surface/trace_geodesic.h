#pragma once

#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/surface_point.h"


namespace geometrycentral {
namespace surface {


struct TraceGeodesicResult {
  SurfacePoint endPoint;                    // the point the trace ended at
  std::vector<SurfacePoint> pathPoints;     // all points along the path, including start and end
  Vector2 endingDir;                        // the incoming direction to the final point
  bool withPath = false;
};


TraceGeodesicResult traceGeodesic(IntrinsicGeometryInterface& geom, SurfacePoint startP, Vector2 traceVec, bool includePath=false);

// For a trace which was expected to end very near targetVertex, try to clean up the end of the path to end directly at targetVertex
void trimTraceResult(TraceGeodesicResult& traceResult, Vertex targetVertex);

}}
