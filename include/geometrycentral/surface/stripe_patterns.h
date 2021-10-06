#pragma once

#include "geometrycentral/surface/embedded_geometry_interface.h"
#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include <vector>

namespace geometrycentral {
namespace surface {

// Implementation of "Stripe Patterns on Surfaces" [Knoppel et al. 2015]

// Takes as input a geometry along with vertex-based frequencies and a line field (2-RoSy) and outputs a 2\pi-periodic
// function defined on triangle corners such that the 0 (mod 2\pi)
// Isolines of this function are stripes following the direction field spaced according to the target frequencies
std::tuple<CornerData<double>, FaceData<int>, FaceData<int>>
computeStripePattern(IntrinsicGeometryInterface& geometry, const VertexData<double>& frequencies,
                     const VertexData<Vector2>& directionField);

struct Isoline {
  std::vector<std::pair<Halfedge, double>> barycenters;
  bool open;
};

// Extracts the zero (mod 2pi) level set of the function values defined on the corners
// and returns a list of barycentric coordinates and their corresponding halfedges
// WARNING this only works if there is no more than one isoline crossing at any edge of the geometry
std::vector<Isoline> extractIsolinesFromStripePattern(IntrinsicGeometryInterface& geometry,
                                                      const CornerData<double>& stripesValues,
                                                      const FaceData<int>& stripesIndices,
                                                      const FaceData<int>& fieldIndices);

// Same as above, but returns a representation suitable for Polyscope rendering (requires access to explicit vertex
// positions)
std::tuple<std::vector<Vector3>, std::vector<std::array<int, 2>>>
extractPolylinesFromStripePattern(EmbeddedGeometryInterface& geometry, const CornerData<double>& stripesValues,
                                  const FaceData<int>& stripesIndices, const FaceData<int>& fieldIndices);

} // namespace surface
} // namespace geometrycentral
