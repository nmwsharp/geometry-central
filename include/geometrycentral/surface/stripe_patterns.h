#pragma once

#include "geometrycentral/surface/intrinsic_geometry_interface.h"


namespace geometrycentral {
namespace surface {

/**
 * This function is an implementation of "Stripe Patterns on Surfaces" [Knoppel et al. 2015]
 * It takes as input a geometry along with vertex-based frequencies and a line field (2-RoSy) and outputs a
 * 2\pi-periodic function defined on triangle corners such that the 0 (mod 2\pi) isolines of this function are stripes
 * following the direction field spaced according to the target frequencies
 */
std::tuple<CornerData<double>, FaceData<int>, FaceData<int>>
computeStripePattern(IntrinsicGeometryInterface& geometry, const VertexData<double>& frequencies,
                     const VertexData<Vector2>& directionField);

} // namespace surface
} // namespace geometrycentral
