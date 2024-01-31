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
  Eigen::SparseMatrix<double> virtualRefinementLaplacian;
  void requireVirtualRefinementPolygonLaplacian();
  void unrequireVirtualRefinementPolygonLaplacian();

  // Vertex Galerkin mass matrix (unlumped)
  Eigen::SparseMatrix<double> virtualRefinementVertexGalerkinMassMatrix;
  void requireVirtualRefinementVertexGalerkinMassMatrix();
  void unrequireVirtualRefinementVertexGalerkinMassMatrix();

  // Vertex mass matrix (lumped)
  Eigen::SparseMatrix<double> virtualRefinementVertexLumpedMassMatrix;
  void requireVirtualRefinementVertexLumpedMassMatrix();
  void unrequireVirtualRefinementVertexLumpedMassMatrix();

  // DEC Operators
  Eigen::SparseMatrix<double> virtualRefinementHodge0, virtualRefinementHodge0Inverse, virtualRefinementHodge1,
      virtualRefinementHodge1Inverse, virtualRefinementHodge2, virtualRefinementHodge2Inverse, virtualRefinementD0,
      virtualRefinementD1;
  void requireVirtualRefinementDECOperators();
  void unrequireVirtualRefinementDECOperators();

  // = de Goes et al. "Discrete Differential Operators on Polygonal Meshes" (2020), based on the virtual element method.

  // Laplacian
  Eigen::SparseMatrix<double> virtualElementLaplacian;
  void requireVirtualElementLaplacian();
  void unrequireVirtualElementLaplacian();

  // Vertex mass matrix (unlumped)
  Eigen::SparseMatrix<double> virtualElementVertexGalerkinMassMatrix;
  void requireVirtualElementVertexGalerkinMassMatrix();
  void requireVirtualElementVertexGalerkinMassMatrix();

  // Vertex mass matrix (lumped)
  Eigen::SparseMatrix<double> virtualElementVertexLumpedMassMatrix;
  void requireVirtualElementVertexLumpedMassMatrix();
  void unrequireVirtualElementVertexLumpedMassMatrix();

  // Vertex connection Laplacian
  Eigen::SparseMatrix<std::complex<double>> virtualElementVertexConnectionLaplacian;
  void requireVirtualElementVertexConnectionLaplacian();
  void unrequireVirtualElementVertexConnectionLaplacian();

  Eigen::SparseMatrix<double> virtualElementHodge0, virtualElementHodge0Inverse, virtualElementHodge1,
      virtualElementHodge1Inverse, virtualElementHodge2, virtualElementHodge2Inverse, virtualElementD0,
      virtualElementD1;
  void requireVirtualElementDECOperators();
  void unrequireVirtualElementDECOperators();

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

  // helper functions
  DependentQuantityD<FaceData<Eigen::VectorXd>> virtualRefinementAreaWeightsQ; // affine weights for each virtual node
  virtual void computeVirtualRefinementAreaWeights();
  virtual Eigen::MatrixXd buildPolygonMassMatrix(const Face& f) const;
  virtual Eigen::MatrixXd buildPolygonStiffnessMatrix(const Face& f) const;
  virtual SparseMatrix<double> buildDivergenceMatrix() const;
  virtual SparseMatrix<double> buildGradientMatrix() const;
  virtual SparseMatrix<double> buildGradientMassMatrix() const;
  virtual SparseMatrix<double> buildProlongationMatrix() const;
  virtual Eigen::MatrixXd getPolygonPositionMatrix(const Face& f) const;
  virtual Eigen::VectorXd computeVirtualVertex(const Eigen::MatrixXd& poly) const;
  virtual Vector3 gradientHatFunction(const Vector3& a, const Vector3& b, const Vector3& c) const;

  // Laplacian
  DependentQuantityD<Eigen::SparseMatrix<double>> virtualRefinementLaplacianQ;
  virtual void computeVirtualRefinementLaplacian();

  // Vertex mass matrix (unlumped)
  DependentQuantityD<Eigen::SparseMatrix<double>> virtualRefinementVertexGalerkinMassMatrixQ;
  virtual void computeVirtualRefinementVertexGalerkinMassMatrix();

  // Vertex mass matrix (lumped)
  DependentQuantityD<Eigen::SparseMatrix<double>> virtualRefinementVertexLumpedMassMatrixQ;
  virtual void computeVirtualRefinementVertexLumpedMassMatrix();

  // DEC Operators
  std::array<Eigen::SparseMatrix<double>*, 8> virtualRefinementDECOperatorArray;
  DependentQuantityD<std::array<Eigen::SparseMatrix<double>*, 8>> virtualRefinementDECOperatorsQ;
  virtual void computeVirtualRefinementDECOperators();

  // = de Goes et al. "Discrete Differential Operators on Polygonal Meshes" (2020), based on the virtual element method.

  // Laplacian
  DependentQuantityD<Eigen::SparseMatrix<double>> virtualElementLaplacianQ;
  virtual void computeVirtualElementLaplacian();

  // Vertex mass matrix (unlumped)
  DependentQuantityD<Eigen::SparseMatrix<double>> virtualElementVertexGalerkinMassMatrixQ;
  virtual void computeVirtualElementVertexGalerkinMassMatrix();

  // Vertex mass matrix (lumped)
  DependentQuantityD<Eigen::SparseMatrix<double>> virtualElementVertexLumpedMassMatrixQ;
  virtual void computeVirtualElementVertexLumpedMassMatrix();

  // Vertex connection Laplacian
  DependentQuantityD<Eigen::SparseMatrix<std::complex<double>>> virtualElementVertexConnectionLaplacianQ;
  virtual void computeVirtualElementVertexConnectionLaplacian();

  // DEC Operators
  std::array<Eigen::SparseMatrix<double>*, 8> virtualElementDECOperatorArray;
  DependentQuantityD<std::array<Eigen::SparseMatrix<double>*, 8>> virtualElementDECOperatorsQ;
  virtual void computeVirtualElementDECOperators();
};


} // namespace surface
} // namespace geometrycentral
