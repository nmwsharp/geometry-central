#pragma once

#include "geometrycentral/surface/barycentric_vector.h"

namespace geometrycentral {
namespace surface {

std::vector<std::vector<SurfacePoint>> marchingTriangles(IntrinsicGeometryInterface& geom, const VertexData<double>& u,
                                                         double isoval = 0.);

}
} // namespace geometrycentral