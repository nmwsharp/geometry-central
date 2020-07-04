#pragma once

#include "geometrycentral/numerical/linear_algebra_utilities.h"
#include "geometrycentral/surface/embedded_geometry_interface.h"
#include "geometrycentral/surface/surface_mesh.h"

namespace geometrycentral {
namespace surface {


std::tuple<SparseMatrix<double>, SparseMatrix<double>>
buildTuftedLaplacian(SurfaceMesh& mesh, EmbeddedGeometryInterface& geom, double relativeMollificationFactor = 0.);

// Modifies the input mesh and edge lengths to be the tufted cover!
void buildIntrinsicTuftedCover(SurfaceMesh& mesh, EdgeData<double>& edgeLengths,
                               EmbeddedGeometryInterface* posGeom = nullptr);


} // namespace surface
} // namespace geometrycentral
