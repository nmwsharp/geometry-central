#pragma once


#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"


namespace geometrycentral {
namespace surface {

  // === High level interface
  
  // Parameterize a disk-like surface (without introducing any cuts or cones)
  VertexData<Vector2> parameterizeDisk(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geometry);
  
  // Parameterize a surface with the specified cut (which must render the surface a disk)
  // TODO not implemented
  CornerData<Vector2> parameterize(SurfaceMesh& mesh, IntrinsicGeometryInterface& geometry, const EdgeData<char>& cut);
  
  // Parameterize a surface with spherical topology over the unit sphere
  VertexData<Vector3> parameterizeSphere(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geometry);

}}
