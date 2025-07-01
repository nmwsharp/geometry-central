#pragma once

#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vector_heat_method.h"

namespace geometrycentral {
namespace surface {

// Find a center of a collection of points at vertices
SurfacePoint findCenter(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom,
                        const std::vector<Vertex>& vertexPts, int p = 2,
                        LogMapStrategy strategy = LogMapStrategy::VectorHeat);
SurfacePoint findCenter(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom, VectorHeatMethodSolver& solver,
                        const std::vector<Vertex>& vertexPts, int p = 2,
                        LogMapStrategy strategy = LogMapStrategy::VectorHeat);

// Find a center of distribution
SurfacePoint findCenter(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom,
                        const VertexData<double>& distribution, int p = 2,
                        LogMapStrategy strategy = LogMapStrategy::VectorHeat);
SurfacePoint findCenter(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom, VectorHeatMethodSolver& solver,
                        const VertexData<double>& distribution, int p = 2,
                        LogMapStrategy strategy = LogMapStrategy::VectorHeat);

} // namespace surface
} // namespace geometrycentral
