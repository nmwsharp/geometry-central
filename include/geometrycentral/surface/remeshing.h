#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include <queue>

namespace geometrycentral {
namespace surface {

// makes all triangles delaunay
void fixDelaunay(SurfaceMesh& mesh, VertexPositionGeometry& geometry);

// average positions of vertices based on surrounding vertex positions
void smoothByLaplacian(SurfaceMesh& mesh, VertexPositionGeometry& geometry);

// average positions of vertices based on surrounding triangle circumenters
void smoothByCircumcenter(SurfaceMesh& mesh, VertexPositionGeometry& geometry);

// applies splits and collapses to adjust edge lengths based on the curvature
// flatLength: specifies how long the target edge length should be in flat regions
// epsilon: controls how much variation in target length occurs due to curvature
// minLength: specifies the minimum possible length an edge could become
bool adjustEdgeLengths(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geometry, double flatLength, double epsilon,
                       double minLength, bool curvatureAdaptive = true);

} // namespace surface
} // namespace geometrycentral
