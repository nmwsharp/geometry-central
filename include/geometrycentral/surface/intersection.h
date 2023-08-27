#pragma once

#include "geometrycentral/surface/embedded_geometry_interface.h"
#include "geometrycentral/utilities/utilities.h"

#include <cmath>
#include <utility>
#include <vector>


namespace geometrycentral {
namespace surface {

struct SurfaceIntersectionResult {
  std::vector<Vector3> points;
  std::vector<std::array<size_t, 2>> edges;
  bool hasIntersections;
};

SurfaceIntersectionResult selfIntersections(EmbeddedGeometryInterface& geometry);
SurfaceIntersectionResult intersections(EmbeddedGeometryInterface& geometry1, EmbeddedGeometryInterface& geometry2,
                                        bool selfCheck = false);

} // namespace surface
} // namespace geometrycentral
