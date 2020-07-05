#pragma once

#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/surface_mesh.h"

namespace geometrycentral {
namespace surface {

// A simple but effective method to resolve geometrically degenerate faces in a mesh, by adding a tiny epsilon to
// intrinsic edge lengths.

// Offset via a factor relative to the mean edge length in the mesh
void mollifyIntrinsic(SurfaceMesh& mesh, EdgeData<double>& edgeLengths, double relativeFactor = 1e-6);

// Offset via an absolute factor
void mollifyIntrinsicAbsolute(SurfaceMesh& mesh, EdgeData<double>& edgeLengths, double absoluteFactor);

} // namespace surface
} // namespace geometrycentral
