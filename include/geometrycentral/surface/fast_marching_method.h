#pragma once

#include "geometrycentral/surface/barycentric_vector.h"
#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/utilities/utilities.h"

#include <cmath>
#include <utility>
#include <vector>


namespace geometrycentral {
namespace surface {

VertexData<double> FMMDistance(IntrinsicGeometryInterface& geometry,
                               const std::vector<std::vector<std::pair<SurfacePoint, double>>>& initialDistances,
                               bool sign = false);

VertexData<double> FMMDistance(IntrinsicGeometryInterface& geometry,
                               const std::vector<std::pair<Vertex, double>>& initialDistances, bool sign = false);


} // namespace surface
} // namespace geometrycentral
