#pragma once


#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/halfedge_mesh.h"


namespace geometrycentral {
namespace surface {

  // === High level interface
  
  // Paramerize a disk-like surface (without introducing any cuts or cones)
  VertexData<Vector2> parameterizeDisk(IntrinsicGeometryInterface& geometry);
  
  // Paramerize a surface with the specified cut (which must render the surface a disk)
  // TODO not implemented
  CornerData<Vector2> parameterize(IntrinsicGeometryInterface& geometry, const EdgeData<char>& cut);


}}
