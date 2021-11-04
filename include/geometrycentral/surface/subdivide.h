#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/mutation_manager.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

namespace geometrycentral {
namespace surface {

// Subdivide to form quad mesh
// Note: these do not take MutationManager arguments since MutationManager only supports triangular meshes at the moment
void linearSubdivide(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geo);
void catmullClarkSubdivide(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geo);

void loopSubdivide(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geo);
void loopSubdivide(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geo, MutationManager& mm);

} // namespace surface
} // namespace geometrycentral
