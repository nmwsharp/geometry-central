#pragma once

#include "geometrycentral/surface/embedded_geometry_interface.h"

namespace geometrycentral {
namespace surface {

// Helper functions for Bunge et al. "Polygon Laplacian Made Simple" (2020).
// Copyright (C) 2020 Astrid Bunge, Philipp Herholz, Misha Kazhdan, Mario Botsch, MIT license
// (Modified to work in geometry-central. Original code can be found here:
// https://github.com/mbotsch/polygon-laplacian)

// FaceData<Eigen::VectorXd> virtualRefinementAreaWeights;
// DependentQuantityD<FaceData<Eigen::VectorXd>> virtualRefinementAreaWeightsQ;
// void computeVirtualRefinementAreaWeights(EmbeddedGeometryInterface& geometry);
FaceData<Eigen::VectorXd> computeVirtualRefinementAreaWeights(EmbeddedGeometryInterface& geometry);
Eigen::MatrixXd simplePolygonMassMatrix(EmbeddedGeometryInterface& geometry,
                                        const FaceData<Eigen::VectorXd>& virtualRefinementAreaWeights, const Face& f);
Eigen::MatrixXd simplePolygonStiffnessMatrix(EmbeddedGeometryInterface& geometry,
                                             const FaceData<Eigen::VectorXd>& virtualRefinementAreaWeights,
                                             const Face& f);
Eigen::MatrixXd polygonPositionMatrix(EmbeddedGeometryInterface& geometry, const Face& f);
Eigen::VectorXd simplePolygonVirtualVertex(const Eigen::MatrixXd& poly);

// Helper functions for de Goes et al. "Discrete Differential Operators on Polygonal Meshes" (2020).
// Use of this source code is governed by a LGPL-3.0 license.
// (Modified to work in geometry-central. Original code can be found here:
// https://github.com/DGtal-team/DGtal/blob/master/src/DGtal/dec/PolygonalCalculus.h)

Eigen::MatrixXd polygonPerFaceLaplacian(EmbeddedGeometryInterface& geometry, const Face& f, const double polygonLambda);
Eigen::MatrixXd polygonPerFaceInnerProductMatrix(EmbeddedGeometryInterface& geometry, const Face& f,
                                                 const double polygonLambda);
Eigen::MatrixXd polygonProjectionMatrix(EmbeddedGeometryInterface& geometry, const Face& f);
Eigen::MatrixXd polygonPerFaceGradientMatrix(EmbeddedGeometryInterface& geometry, const Face& f);
Eigen::MatrixXd polygonCoGradientMatrix(EmbeddedGeometryInterface& geometry, const Face& f);
Eigen::MatrixXd polygonAveragingMatrix(const Face& f);
Eigen::MatrixXd polygonDerivativeMatrix(const Face& f);
Eigen::MatrixXd polygonEdgeVectorMatrix(EmbeddedGeometryInterface& geometry, const Face& f);
Eigen::MatrixXd polygonEdgeMidpointMatrix(EmbeddedGeometryInterface& geometry, const Face& f);
Eigen::MatrixXd polygonFlat(EmbeddedGeometryInterface& geometry, const Face& f);
Eigen::MatrixXd polygonSharp(EmbeddedGeometryInterface& geometry, const Face& f);
Eigen::Vector3d polygonCentroid(EmbeddedGeometryInterface& geometry, const Face& f);
// connections
Eigen::MatrixXd polygonPerFaceConnectionLaplacian(EmbeddedGeometryInterface& geometry, const Face& f,
                                                  const double polygonLambda);
Eigen::MatrixXd polygonBlockConnection(EmbeddedGeometryInterface& geometry, const Face& f);
Eigen::MatrixXd polygonCovariantGradient(EmbeddedGeometryInterface& geometry, const Face& f);
Eigen::MatrixXd polygonCovariantProjection(EmbeddedGeometryInterface& geometry, const Face& f);
// tangent space helpers
Eigen::MatrixXd Tv(EmbeddedGeometryInterface& geometry, const Vertex& v);
Eigen::MatrixXd Tf(EmbeddedGeometryInterface& geometry, const Face& f);
Eigen::Matrix2d Rvf(EmbeddedGeometryInterface& geometry, const Vertex& v, const Face& f);
Eigen::Matrix3d Qvf(EmbeddedGeometryInterface& geometry, const Vertex& v, const Face& f);
// helpers to the helper functions: generic linear algebra stuff, though probably wouldn't find much use elsewhere
// so keeping them here -- also they use Eigen::Vectors here for matrix-multiply compatibility.
Eigen::Matrix3d bracket(const Eigen::Vector3d& n);
Eigen::Vector3d project(const Eigen::Vector3d& u, const Eigen::Vector3d& n);
Eigen::MatrixXd kroneckerWithI2(const Eigen::MatrixXd& M);

} // namespace surface
} // namespace geometrycentral