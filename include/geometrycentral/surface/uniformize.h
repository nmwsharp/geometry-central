#pragma once


#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"


namespace geometrycentral {
namespace surface {

// Uniformize a topological disk to the flat plane w/ "natural" (u = 0 Dirichlet) boundary conditions
EdgeData<double> uniformizeDisk(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geometry,
                                bool withEdgeFlips = true);


} // namespace surface
} // namespace geometrycentral
