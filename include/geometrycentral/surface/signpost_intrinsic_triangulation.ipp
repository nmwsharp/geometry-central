#pragma once

#include "geometrycentral/surface/signpost_intrinsic_triangulation.h"

namespace geometrycentral {
namespace surface {

inline double SignpostIntrinsicTriangulation::standardizeAngle(Vertex vert, double angle) const {
  if (vert.isBoundary()) {
    // can't wrap around at vertices
    return angle;
  }
  return std::fmod(angle, vertexAngleSums[vert]);
}

inline Vector2 SignpostIntrinsicTriangulation::halfedgeVector(Halfedge he) const {
  double edgeAngle = signpostAngle[he];
  double scaleFac = 1.0 / vertexAngleScaling(he.vertex());
  Vector2 traceVec = Vector2::fromAngle(edgeAngle * scaleFac) * edgeLengths[he.edge()];
  return traceVec;
}

inline Vector2 SignpostIntrinsicTriangulation::rescaledVertexVector(Vertex v, double angle, double len) const {
  double scaleFac = 1.0 / vertexAngleScaling(v);
  Vector2 traceVec = Vector2::fromAngle(angle * scaleFac) * len;
  return traceVec;
}

inline double SignpostIntrinsicTriangulation::vertexAngleScaling(Vertex v) const {
  return vertexAngleSums[v] / (v.isBoundary() ? M_PI : 2. * M_PI);
}


} // namespace surface
} // namespace geometrycentral
