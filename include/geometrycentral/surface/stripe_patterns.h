#pragma once

#include "geometrycentral/surface/embedded_geometry_interface.h"
#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include <vector>

namespace geometrycentral {
namespace surface {

// Implementation of "Stripe Patterns on Surfaces" [Knoppel et al. 2015]
// Based on Keenan Crane's original implementation available here:
// https://www.cs.cmu.edu/~kmcrane/Projects/StripePatterns/code.zip

// Takes as input a geometry along with vertex-based frequencies and a line field (2-RoSy) and outputs a 2\pi-periodic
// function defined on triangle corners such that the 0 (mod 2\pi)
// Isolines of this function are stripes perpendicular to the direction field spaced according to the target frequencies
std::tuple<CornerData<double>, FaceData<int>, FaceData<int>>
computeStripePattern(IntrinsicGeometryInterface& geometry, const VertexData<double>& frequencies,
                     const VertexData<Vector2>& directionField);

// Extracts the zero (mod 2pi) level set of the function values defined on the corners
// Returns a list of vertices and edges suitable for rendering (e.g. with Polyscope)
// (requires access to explicit vertex positions)
std::tuple<std::vector<Vector3>, std::vector<std::array<size_t, 2>>>
extractPolylinesFromStripePattern(EmbeddedGeometryInterface& geometry, const CornerData<double>& values,
                                  const FaceData<int>& stripeIndices, const FaceData<int>& fieldIndices,
                                  const VertexData<Vector2>& directionField, bool connectOnSingularities);

// Runs both of the above functions (per-corner data computation and polyline extraction)
// can optionally connect isolines separated by a singularity using a directionField alignment heuristic
std::tuple<std::vector<Vector3>, std::vector<std::array<size_t, 2>>>
computeStripePatternPolylines(EmbeddedGeometryInterface& geometry, const VertexData<double>& frequencies,
                              const VertexData<Vector2>& directionField, bool connectIsolinesOnSingularities = true);

} // namespace surface
} // namespace geometrycentral
