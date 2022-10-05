#pragma once

#include "geometrycentral/surface/barycentric_coordinate_helpers.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/mutation_manager.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include <deque>

namespace geometrycentral {
namespace surface {

// makes all triangles delaunay
// returns the number of flips performed
size_t fixDelaunay(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geometry);
size_t fixDelaunay(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geometry, MutationManager& mm);

// average positions of vertices based on surrounding vertex positions
// returns the average amount each vertex was moved by
double smoothByLaplacian(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geometry, double stepSize = 1);
double smoothByLaplacian(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geometry, MutationManager& mm,
                         double stepSize = 1);

// average positions of vertices based on surrounding triangle circumenters
// returns the average amount each vertex was moved by
double smoothByCircumcenter(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geometry, double stepSize = 1);
double smoothByCircumcenter(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geometry, MutationManager& mm,
                            double stepSize = 1);

// applies splits and collapses to adjust edge lengths based on the curvature
// flatLength:          specifies how long the target edge length should be in flat regions
// curvatureAdaptation: controls how much variation in target length occurs due to curvature. Set curvatureAdaptation to
//                      0 if you want lengths to be approximately flatLength everywhere
// minRelativeLength:   specifies the minimum possible length an edge could become. Defined relative to flatLength
bool adjustEdgeLengths(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geometry, double flatLength,
                       double curvatureAdaptation = 0.2, double minRelativeLength = 0.05);
bool adjustEdgeLengths(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geometry, MutationManager& mm,
                       double flatLength, double curvatureAdaptation = 0., double minRelativeLength = 0.05);

void remesh(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geometry, double targetEdgeLength,
            size_t maxIterations = 10, double curvatureAdaptation = 0., double minRelativeLength = 0.05);
void remesh(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geometry, MutationManager& mm, double targetEdgeLength,
            size_t maxIterations = 10, double curvatureAdaptation = 0., double minRelativeLength = 0.05);

} // namespace surface
} // namespace geometrycentral
