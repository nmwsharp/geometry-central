#pragma once

#include "geometrycentral/surface/embedded_geometry_interface.h"

namespace geometrycentral {
namespace surface {

// Helper functions for Bunge et al. "Polygon Laplacian Made Simple" (2020).
// Copyright (C) 2020 Astrid Bunge, Philipp Herholz, Misha Kazhdan, Mario Botsch, MIT license
// (Modified to work in geometry-central. Original code can be found here:
// https://github.com/mbotsch/polygon-laplacian)

Eigen::VectorXd simplePolygonVirtualVertex(const Eigen::MatrixXd& poly);

// Helper functions for de Goes et al. "Discrete Differential Operators on Polygonal Meshes" (2020).
// Use of this source code is governed by a LGPL-3.0 license.
// (Modified to work in geometry-central. Original code can be found here:
// https://github.com/DGtal-team/DGtal/blob/master/src/DGtal/dec/PolygonalCalculus.h)

Eigen::MatrixXd polygonAveragingMatrix(const Face& f);

Eigen::MatrixXd polygonDerivativeMatrix(const Face& f);

// helpers to the helper functions: generic linear algebra stuff, though probably wouldn't find much use elsewhere
// so keeping them here -- also they use Eigen::Vectors here for matrix-multiply compatibility.
Eigen::Matrix3d bracket(const Eigen::Vector3d& n);

Eigen::Vector3d project(const Eigen::Vector3d& u, const Eigen::Vector3d& n);

Eigen::MatrixXd kroneckerWithI2(const Eigen::MatrixXd& M);

} // namespace surface
} // namespace geometrycentral