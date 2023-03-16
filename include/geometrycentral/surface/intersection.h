#pragma once

#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/utilities/utilities.h"

#include <cmath>
#include <utility>
#include <vector>


namespace geometrycentral {
namespace surface {

struct SurfaceIntersectionResult {
   std::vector<Vector3> points;
   std::vector<std::array<size_t,2>> edges;
   bool hasIntersections;
};

SurfaceIntersectionResult selfIntersections(VertexPositionGeometry& geometry);
SurfaceIntersectionResult intersections(VertexPositionGeometry& geometry1,
                                        VertexPositionGeometry& geometry2,
                                        bool selfCheck = false );

} // namespace surface
} // namespace geometrycentral
