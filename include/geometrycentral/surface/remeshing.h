#pragma once

#include "geometrycentral/surface/barycentric_coordinate_helpers.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/mutation_manager.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include <deque>

namespace geometrycentral {
namespace surface {

enum class RemeshBoundaryCondition { Fixed, Tangential, Free };
enum class RemeshSmoothStyle { Circumcentric, Laplacian };

struct RemeshOptions {
  double targetEdgeLength = -1; // the target edge length in flat regions. If `targetEdgeLength` is negative, the target
                                // edge length is set to relative the input mesh's mean edge length
  size_t maxIterations = 10;    // the maximum number of iterations to run for
  double curvatureAdaptation = 0;  // how much target length should vary due to curvature. Set curvatureAdaptation
                                   // to 0 if you want lengths to be approximately targetEdgeLength everywhere
  double minRelativeLength = 0.05; // the minimum possible edge length allowed in the output mesh. Defined relative to
                                   // targetEdgeLength
  RemeshSmoothStyle smoothStyle = RemeshSmoothStyle::Circumcentric; // smoothing function to use
  RemeshBoundaryCondition boundaryCondition =
      RemeshBoundaryCondition::Tangential; // allowed movement of boundary vertices
};
extern const RemeshOptions defaultRemeshOptions;

// Improve mesh using repeated rounds of edge flipping, vertex position smoothing, and edge splits/collapses
void remesh(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, RemeshOptions options = defaultRemeshOptions);
void remesh(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, MutationManager& mm,
            RemeshOptions options = defaultRemeshOptions);

// Try to make all triangles Delaunay
// Returns the number of flips performed
size_t fixDelaunay(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom);
size_t fixDelaunay(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, MutationManager& mm);

// Average positions of vertices based on surrounding vertex positions
// Returns the average amount each vertex was moved by
double smoothByLaplacian(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, double stepSize = 1,
                         RemeshBoundaryCondition bc = RemeshBoundaryCondition::Tangential);
double smoothByLaplacian(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, MutationManager& mm,
                         double stepSize = 1, RemeshBoundaryCondition bc = RemeshBoundaryCondition::Tangential);

// Average positions of vertices based on surrounding triangle circumenters as in [Chen & Holst 2011]
// Returns the average amount each vertex was moved by
double smoothByCircumcenter(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, double stepSize = 1,
                            RemeshBoundaryCondition bc = RemeshBoundaryCondition::Tangential);
double smoothByCircumcenter(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, MutationManager& mm,
                            double stepSize = 1, RemeshBoundaryCondition bc = RemeshBoundaryCondition::Tangential);

// applies splits and collapses to adjust edge lengths based on the curvature
bool adjustEdgeLengths(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom,
                       RemeshOptions options = defaultRemeshOptions);
bool adjustEdgeLengths(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, MutationManager& mm,
                       RemeshOptions options = defaultRemeshOptions);

} // namespace surface
} // namespace geometrycentral
