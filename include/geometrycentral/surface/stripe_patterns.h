#pragma once

#include "geometrycentral/surface/intrinsic_geometry_interface.h"


namespace geometrycentral {
namespace surface {

std::tuple<CornerData<double>, FaceData<int>, FaceData<int>>
computeStripePattern(IntrinsicGeometryInterface& geometry, const VertexData<double>& frequencies,
                     const VertexData<Vector2>& directionField);

} // namespace surface
} // namespace geometrycentral
