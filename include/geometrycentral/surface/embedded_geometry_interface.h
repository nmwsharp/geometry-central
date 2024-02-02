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

  FaceData<Vector<double>> virtualRefinementAreaWeights;

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
  DependentQuantityD<FaceData<Eigen::VectorXd>> virtualRefinementAreaWeightsQ; // affine weights for each virtual node
  virtual void computeVirtualRefinementAreaWeights();
  virtual Eigen::MatrixXd simplePolygonMassMatrix(const Face& f) const;
  virtual Eigen::MatrixXd simplePolygonStiffnessMatrix(const Face& f) const;
  virtual SparseMatrix<double> simplePolygonProlongationMatrix() const;
  virtual Eigen::MatrixXd polygonPositionMatrix(const Face& f) const;
  virtual Eigen::VectorXd simplePolygonVirtualVertex(const Eigen::MatrixXd& poly) const;
  virtual Vector3 gradientHatFunction(const Vector3& a, const Vector3& b, const Vector3& c) const; // TODO: axe?

  // = de Goes et al. "Discrete Differential Operators on Polygonal Meshes" (2020), based on the virtual element method.

  // Laplacian
  DependentQuantityD<Eigen::SparseMatrix<double>> polygonLaplacianQ;
  virtual void computePolygonLaplacian();

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
  virtual Eigen::MatrixXd polygonPerFaceLaplacian(const Face& f) const;
  virtual Eigen::MatrixXd polygonPerFaceInnerProductMatrix(const Face& f) const;
  virtual Eigen::MatrixXd polygonProjectionMatrix(const Face& f) const;
  virtual Eigen::MatrixXd polygonCoGradientMatrix(const Face& f) const;
  virtual Eigen::MatrixXd polygonGradientMatrix(const Face& f) const;
  virtual Eigen::MatrixXd polygonAveragingMatrix(const Face& f) const;
  virtual Eigen::MatrixXd polygonDerivativeMatrix(const Face& f) const;
  virtual Eigen::MatrixXd polygonEdgeVectorMatrix(const Face& f) const;
  virtual Eigen::MatrixXd polygonEdgeMidpointMatrix(const Face& f) const;
  virtual Eigen::MatrixXd polygonFlat(const Face& f) const;
  virtual Eigen::MatrixXd polygonSharp(const Face& f) const;
  virtual Eigen::Vector3d polygonVectorArea(const Face& f) const;
  virtual double polygonArea(const Face& f) const;
  virtual Eigen::Vector3d polygonNormal(const Face& f) const;
  virtual Eigen::Vector3d polygonCentroid(const Face& f) const;
  // connections
  virtual Eigen::MatrixXd polygonPerFaceConnectionLaplacian(const Face& f) const;
  virtual Eigen::MatrixXd polygonBlockConnection(const Face& f) const;
  virtual Eigen::MatrixXd polygonCovariantGradient(const Face& f) const;
  virtual Eigen::MatrixXd polygonCovariantProjection(const Face& f) const;
  // tangent space helpers
  virtual Eigen::MatrixXd Tv(const Vertex& v) const;
  virtual Eigen::MatrixXd Tf(const Face& f) const;
  virtual Eigen::Matrix2d Rvf(const Vertex& v, const Face& f) const;
  virtual Eigen::Matrix3d Qvf(const Vertex& v, const Face& f) const;
  // helpers to the helper functions: generic linear algebra stuff, though probably wouldn't find much use elsewhere
  // so keeping them here -- also they use Eigen::Vectors here for matrix-multiply compatibility.
  virtual Eigen::Matrix3d bracket(const Eigen::Vector3d& n) const;
  virtual Eigen::Vector3d project(const Eigen::Vector3d& u, const Eigen::Vector3d& n) const;
  virtual Eigen::MatrixXd kroneckerWithI2(const Eigen::MatrixXd& M) const;
};


} // namespace surface
} // namespace geometrycentral
