#pragma once


#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/intrinsic_geometry_interface.h"


namespace geometrycentral {
namespace surface {

// A simplified interface for flip-based intrinsic Delaunay triangulations, which only supports a length-based
// representation. See signpost_intrinsic_triangulation.h for much more advanced functionality

// Modifies both the underlying mesh and edge lengths to make the intrinsic Delaunay
enum class FlipType { Euclidean = 0, Hyperbolic };
size_t flipToDelaunay(SurfaceMesh& mesh, EdgeData<double>& edgeLengths, FlipType flipType = FlipType::Euclidean,
                      double delaunayEPS = 1e-6);

} // namespace surface
} // namespace geometrycentral
