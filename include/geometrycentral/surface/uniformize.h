#pragma once


#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/halfedge_mesh.h"


namespace geometrycentral {
namespace surface {

// Uniformize a topological disk to the flat plane w/ "natural" (u = 0 Dirichlet) boundary conditions
EdgeData<double> uniformizeDisk(IntrinsicGeometryInterface& geometry, bool withEdgeFlips=true);


}}
