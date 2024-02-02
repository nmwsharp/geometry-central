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

  // Laplacian
  Eigen::SparseMatrix<double> simplePolygonLaplacian;
  void requireSimplePolygonLaplacian();
  void unrequireSimplePolygonLaplacian();

  // Vertex Galerkin mass matrix (unlumped)
  Eigen::SparseMatrix<double> simplePolygonVertexGalerkinMassMatrix;
  void requireSimplePolygonVertexGalerkinMassMatrix();
  void unrequireSimplePolygonVertexGalerkinMassMatrix();

  // Vertex mass matrix (lumped)
  Eigen::SparseMatrix<double> simplePolygonVertexLumpedMassMatrix;
  void requireSimplePolygonVertexLumpedMassMatrix();
  void unrequireSimplePolygonVertexLumpedMassMatrix();

  // = de Goes et al. "Discrete Differential Operators on Polygonal Meshes" (2020), based on the virtual element method.

  // Laplacian
  Eigen::SparseMatrix<double> polygonLaplacian;
  void requirePolygonLaplacian();
  void unrequirePolygonLaplacian();

  // Gradient matrix
  Eigen::SparseMatrix<double> polygonGradientMatrix;
  void requirePolygonGradientMatrix();
  void unrequirePolygonGradientMatrix();

  // Divergence matrix
  Eigen::SparseMatrix<double> polygonDivergenceMatrix;
  void requirePolygonDivergenceMatrix();
  void unrequirePolygonDivergenceMatrix();

  // Vertex mass matrix (lumped)
  Eigen::SparseMatrix<double> polygonVertexLumpedMassMatrix;
  void requirePolygonVertexLumpedMassMatrix();
  void unrequirePolygonVertexLumpedMassMatrix();

  // Vertex connection Laplacian
  Eigen::SparseMatrix<std::complex<double>> polygonVertexConnectionLaplacian;
  void requirePolygonVertexConnectionLaplacian();
  void unrequirePolygonVertexConnectionLaplacian();

  Eigen::SparseMatrix<double> polygonHodge0, polygonHodge0Inverse, polygonHodge1, polygonHodge2, polygonHodge2Inverse,
      polygonD0, polygonD1;
  void requirePolygonDECOperators();
  void unrequirePolygonDECOperators();

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

  // Laplacian
  DependentQuantityD<Eigen::SparseMatrix<double>> simplePolygonLaplacianQ;
  virtual void computeSimplePolygonLaplacian();

  // Vertex mass matrix (unlumped)
  DependentQuantityD<Eigen::SparseMatrix<double>> simplePolygonVertexGalerkinMassMatrixQ;
  virtual void computeSimplePolygonVertexGalerkinMassMatrix();

  // Vertex mass matrix (lumped)
  DependentQuantityD<Eigen::SparseMatrix<double>> simplePolygonVertexLumpedMassMatrixQ;
  virtual void computeSimplePolygonVertexLumpedMassMatrix();

  // helper functions
  FaceData<Eigen::VectorXd> virtualRefinementAreaWeights;
  DependentQuantityD<FaceData<Eigen::VectorXd>> virtualRefinementAreaWeightsQ; // affine weights for each virtual node
  virtual void computeVirtualRefinementAreaWeights();
  virtual Eigen::MatrixXd simplePolygonMassMatrix(const Face& f);
  virtual Eigen::MatrixXd simplePolygonStiffnessMatrix(const Face& f);
  virtual SparseMatrix<double> simplePolygonProlongationMatrix();
  virtual Eigen::MatrixXd polygonPositionMatrix(const Face& f);
  virtual Eigen::VectorXd simplePolygonVirtualVertex(const Eigen::MatrixXd& poly) const;

  // = de Goes et al. "Discrete Differential Operators on Polygonal Meshes" (2020), based on the virtual element method.

  // Laplacian
  DependentQuantityD<Eigen::SparseMatrix<double>> polygonLaplacianQ;
  virtual void computePolygonLaplacian();

  DependentQuantityD<Eigen::SparseMatrix<double>> polygonGradientMatrixQ;
  virtual void computePolygonGradientMatrix();

  DependentQuantityD<Eigen::SparseMatrix<double>> polygonDivergenceMatrixQ;
  virtual void computePolygonDivergenceMatrix();

  // Vertex mass matrix (lumped)
  DependentQuantityD<Eigen::SparseMatrix<double>> polygonVertexLumpedMassMatrixQ;
  virtual void computePolygonVertexLumpedMassMatrix();

  // Vertex connection Laplacian
  DependentQuantityD<Eigen::SparseMatrix<std::complex<double>>> polygonVertexConnectionLaplacianQ;
  virtual void computePolygonVertexConnectionLaplacian();

  // DEC Operators
  std::array<Eigen::SparseMatrix<double>*, 7> polygonDECOperatorArray;
  DependentQuantityD<std::array<Eigen::SparseMatrix<double>*, 7>> polygonDECOperatorsQ;
  virtual void computePolygonDECOperators();

  // helper functions
  const double polygonLambda = 1.0;
  VertexData<Eigen::VectorXd> polygonVertexNormals;
  DependentQuantityD<VertexData<Eigen::VectorXd>> polygonVertexNormalsQ;
  virtual void computePolygonVertexNormals();
  virtual Eigen::MatrixXd polygonPerFaceLaplacian(const Face& f);
  virtual Eigen::MatrixXd polygonPerFaceInnerProductMatrix(const Face& f);
  virtual Eigen::MatrixXd polygonProjectionMatrix(const Face& f);
  virtual Eigen::MatrixXd polygonPerFaceGradientMatrix(const Face& f);
  virtual Eigen::MatrixXd polygonCoGradientMatrix(const Face& f);
  virtual Eigen::MatrixXd polygonAveragingMatrix(const Face& f) const;
  virtual Eigen::MatrixXd polygonDerivativeMatrix(const Face& f) const;
  virtual Eigen::MatrixXd polygonEdgeVectorMatrix(const Face& f);
  virtual Eigen::MatrixXd polygonEdgeMidpointMatrix(const Face& f);
  virtual Eigen::MatrixXd polygonFlat(const Face& f);
  virtual Eigen::MatrixXd polygonSharp(const Face& f);
  virtual Eigen::Vector3d polygonVectorArea(const Face& f);
  virtual double polygonArea(const Face& f);
  virtual Eigen::Vector3d polygonNormal(const Face& f);
  virtual Eigen::Vector3d polygonCentroid(const Face& f);
  // connections
  virtual Eigen::MatrixXd polygonPerFaceConnectionLaplacian(const Face& f);
  virtual Eigen::MatrixXd polygonBlockConnection(const Face& f);
  virtual Eigen::MatrixXd polygonCovariantGradient(const Face& f);
  virtual Eigen::MatrixXd polygonCovariantProjection(const Face& f);
  // tangent space helpers
  virtual Eigen::MatrixXd Tv(const Vertex& v);
  virtual Eigen::MatrixXd Tf(const Face& f);
  virtual Eigen::Matrix2d Rvf(const Vertex& v, const Face& f);
  virtual Eigen::Matrix3d Qvf(const Vertex& v, const Face& f);
  // helpers to the helper functions: generic linear algebra stuff, though probably wouldn't find much use elsewhere
  // so keeping them here -- also they use Eigen::Vectors here for matrix-multiply compatibility.
  virtual Eigen::Matrix3d bracket(const Eigen::Vector3d& n) const;
  virtual Eigen::Vector3d project(const Eigen::Vector3d& u, const Eigen::Vector3d& n) const;
  virtual Eigen::MatrixXd kroneckerWithI2(const Eigen::MatrixXd& M) const;
};


} // namespace surface
} // namespace geometrycentral
