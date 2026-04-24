#pragma once

#include "geometrycentral/surface/extrinsic_geometry_interface.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/utilities/vector3.h"

#include <Eigen/SparseCore>

namespace geometrycentral {
namespace surface {


class EmbeddedGeometryInterface : public ExtrinsicGeometryInterface {

protected:
  // Constructor is protected, because this class is an interface which is not meant to be instantiated directly.
  // Instantiate it via some realization which encapsulates input data, like VertexPositionGeometry.
  EmbeddedGeometryInterface(SurfaceMesh& mesh_);

public:
  virtual ~EmbeddedGeometryInterface() {}

  // == Quantities

  // Vertex positions
  VertexData<Vector3> vertexPositions;
  void requireVertexPositions();
  void unrequireVertexPositions();

  // Face normal
  FaceData<Vector3> faceNormals;
  void requireFaceNormals();
  void unrequireFaceNormals();

  // Vertex normal
  VertexData<Vector3> vertexNormals;
  void requireVertexNormals();
  void unrequireVertexNormals();

  // Face tangent basis
  FaceData<std::array<Vector3, 2>> faceTangentBasis;
  void requireFaceTangentBasis();
  void unrequireFaceTangentBasis();

  // Vertex tangent basis
  VertexData<std::array<Vector3, 2>> vertexTangentBasis;
  void requireVertexTangentBasis();
  void unrequireVertexTangentBasis();

  // Vertex mean curvature normals
  // These are defined by the property that the mean curvature normals are the laplacian of the vertex positions
  // WARNING: this means that vertexMeanCurvatures != vertexMeanCurvatureNormals.norm()
  VertexData<Vector3> vertexDualMeanCurvatureNormals;
  void requireVertexDualMeanCurvatureNormals();
  void unrequireVertexDualMeanCurvatureNormals();

  // == Polygon Operators

  // = Bunge et al. "Polygon Laplacian Made Simple" (2020), based on virtual refinement (virtual node method).
  // Copyright (C) 2020 Astrid Bunge, Philipp Herholz, Misha Kazhdan, Mario Botsch, MIT license
  // (Modified to work in geometry-central. Original code can be found here:
  // https://github.com/mbotsch/polygon-laplacian)

  // Laplacian
  Eigen::SparseMatrix<double> simplePolygonLaplacian;
  void requireSimplePolygonLaplacian();
  void unrequireSimplePolygonLaplacian();

  // Divergence
  Eigen::SparseMatrix<double> simplePolygonDivergenceMatrix;
  void requireSimplePolygonDivergenceMatrix();
  void unrequireSimplePolygonDivergenceMatrix();

  // Gradient
  Eigen::SparseMatrix<double> simplePolygonGradientMatrix;
  void requireSimplePolygonGradientMatrix();
  void unrequireSimplePolygonGradientMatrix();

  Eigen::SparseMatrix<double> simplePolygonProlongationMatrix;
  void requireSimplePolygonProlongationMatrix();
  void unrequireSimplePolygonProlongationMatrix();

  // connection Laplacian
  Eigen::SparseMatrix<std::complex<double>> simplePolygonVertexConnectionLaplacian;
  void requireSimplePolygonVertexConnectionLaplacian();
  void unrequireSimplePolygonVertexConnectionLaplacian();

  // Vertex Galerkin mass matrix (unlumped)
  Eigen::SparseMatrix<double> simplePolygonVertexGalerkinMassMatrix;
  void requireSimplePolygonVertexGalerkinMassMatrix();
  void unrequireSimplePolygonVertexGalerkinMassMatrix();

  // Vertex mass matrix (lumped)
  Eigen::SparseMatrix<double> simplePolygonVertexLumpedMassMatrix;
  void requireSimplePolygonVertexLumpedMassMatrix();
  void unrequireSimplePolygonVertexLumpedMassMatrix();

protected:
  // == Implmentations of quantities from base classes
  virtual void computeEdgeLengths() override;
  virtual void computeEdgeDihedralAngles() override;

  // == Quantities

  DependentQuantityD<VertexData<Vector3>> vertexPositionsQ;
  virtual void computeVertexPositions() = 0;

  DependentQuantityD<FaceData<Vector3>> faceNormalsQ;
  virtual void computeFaceNormals();

  DependentQuantityD<VertexData<Vector3>> vertexNormalsQ;
  virtual void computeVertexNormals();

  DependentQuantityD<FaceData<std::array<Vector3, 2>>> faceTangentBasisQ;
  virtual void computeFaceTangentBasis();

  DependentQuantityD<VertexData<std::array<Vector3, 2>>> vertexTangentBasisQ;
  virtual void computeVertexTangentBasis();

  DependentQuantityD<VertexData<Vector3>> vertexDualMeanCurvatureNormalsQ;
  virtual void computeVertexDualMeanCurvatureNormals();

  // == Overrides to compute things better using vertex positions
  virtual void computeFaceAreas() override;
  virtual void computeCornerAngles() override;
  virtual void computeHalfedgeCotanWeights() override;
  virtual void computeEdgeCotanWeights() override;

  // == Polygon Operators

  // = Bunge et al. "Polygon Laplacian Made Simple" (2020), based on virtual refinement (virtual node method).
  // Copyright (C) 2020 Astrid Bunge, Philipp Herholz, Misha Kazhdan, Mario Botsch, MIT license
  // (Modified to work in geometry-central. Original code can be found here:
  // https://github.com/mbotsch/polygon-laplacian)

  // Laplacian
  DependentQuantityD<Eigen::SparseMatrix<double>> simplePolygonLaplacianQ;
  virtual void computeSimplePolygonLaplacian();

  // Divergence
  DependentQuantityD<Eigen::SparseMatrix<double>> simplePolygonDivergenceMatrixQ;
  virtual void computeSimplePolygonDivergenceMatrix();

  // Gradient
  DependentQuantityD<Eigen::SparseMatrix<double>> simplePolygonGradientMatrixQ;
  virtual void computeSimplePolygonGradientMatrix();

  // Prolongation
  DependentQuantityD<Eigen::SparseMatrix<double>> simplePolygonProlongationMatrixQ;
  virtual void computeSimplePolygonProlongationMatrix();

  // Connection Laplacian
  DependentQuantityD<Eigen::SparseMatrix<std::complex<double>>> simplePolygonVertexConnectionLaplacianQ;
  virtual void computeSimplePolygonVertexConnectionLaplacian();

  // Vertex mass matrix (unlumped)
  DependentQuantityD<Eigen::SparseMatrix<double>> simplePolygonVertexGalerkinMassMatrixQ;
  virtual void computeSimplePolygonVertexGalerkinMassMatrix();

  // Vertex mass matrix (lumped)
  DependentQuantityD<Eigen::SparseMatrix<double>> simplePolygonVertexLumpedMassMatrixQ;
  virtual void computeSimplePolygonVertexLumpedMassMatrix();

  // helper functions
  Eigen::VectorXd simplePolygonVirtualVertex(const Eigen::MatrixXd& poly);
  Eigen::Vector3d gradientHatFunction(const Eigen::Vector3d& a, const Eigen::Vector3d& b,
                                      const Eigen::Vector3d& c) const;

  // helper functions -- these all depend on quantities in EmbeddedGeometryInterface, which makes them hard to separate
  // their declarations into a separate file.
  FaceData<Eigen::VectorXd> virtualRefinementAreaWeights;
  FaceData<Eigen::Vector3d> virtualRefinementAreaPoints;
  DependentQuantityD<FaceData<Eigen::VectorXd>> virtualRefinementAreaWeightsQ; // affine weights for each virtual node
  DependentQuantityD<FaceData<Eigen::Vector3d>> virtualRefinementAreaPointsQ;
  virtual void computeVirtualRefinementAreaWeights();
  virtual Eigen::MatrixXd simplePolygonMassMatrix(const Face& f);
  virtual Eigen::MatrixXd simplePolygonStiffnessMatrix(const Face& f);
  virtual SparseMatrix<double> simplePolygonGradientMassMatrix();
  virtual Eigen::MatrixXd polygonPositionMatrix(const Face& f);
};


} // namespace surface
} // namespace geometrycentral
