#pragma once

#include "geometrycentral/surface/base_geometry_interface.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/utilities/vector2.h"

#include <Eigen/SparseCore>

#include <complex>

namespace geometrycentral {
namespace surface {


class IntrinsicGeometryInterface : public BaseGeometryInterface {

protected:
  // Constructor is protected, because this class is an interface which is not meant to be instantiated directly.
  // Instantiate it via some realization which encapsulates input data, like EdgeLengthGeometry or
  // VertexPositionGeometry.
  IntrinsicGeometryInterface(SurfaceMesh& mesh_);

public:
  virtual ~IntrinsicGeometryInterface() {}

  // == Lengths, areas, and angles

  // Edge lengths
  EdgeData<double> edgeLengths;
  void requireEdgeLengths();
  void unrequireEdgeLengths();

  // Face areas
  FaceData<double> faceAreas;
  void requireFaceAreas();
  void unrequireFaceAreas();

  // Vertex dual areas
  VertexData<double> vertexDualAreas;
  void requireVertexDualAreas();
  void unrequireVertexDualAreas();

  // Corner angles
  CornerData<double> cornerAngles;
  void requireCornerAngles();
  void unrequireCornerAngles();

  // Vertex angle sums
  VertexData<double> vertexAngleSums;
  void requireVertexAngleSums();
  void unrequireVertexAngleSums();

  // Corner scaled angles
  CornerData<double> cornerScaledAngles;
  void requireCornerScaledAngles();
  void unrequireCornerScaledAngles();

  // Vertex gaussian curvature
  VertexData<double> vertexGaussianCurvatures;
  void requireVertexGaussianCurvatures();
  void unrequireVertexGaussianCurvatures();

  // Face gaussian curvature
  FaceData<double> faceGaussianCurvatures;
  void requireFaceGaussianCurvatures();
  void unrequireFaceGaussianCurvatures();

  // Halfedge cotan weight
  HalfedgeData<double> halfedgeCotanWeights;
  void requireHalfedgeCotanWeights();
  void unrequireHalfedgeCotanWeights();

  // Edge cotan weight
  EdgeData<double> edgeCotanWeights;
  void requireEdgeCotanWeights();
  void unrequireEdgeCotanWeights();

  // Shape length scale
  // (computed as sqrt(total_area), so it is a property of the shape, not the mesh)
  double shapeLengthScale = -1;
  void requireShapeLengthScale();
  void unrequireShapeLengthScale();

  // Mesh length scale
  // (computed as mean edge length, so it is a property of the mesh moreso than the shape)
  double meshLengthScale = -1;
  void requireMeshLengthScale();
  void unrequireMeshLengthScale();


  // == Tangent vectors and transport

  // Halfedge vectors in face tangent space
  HalfedgeData<Vector2> halfedgeVectorsInFace;
  void requireHalfedgeVectorsInFace();
  void unrequireHalfedgeVectorsInFace();

  // Face tangent vector transport across halfedges
  HalfedgeData<Vector2> transportVectorsAcrossHalfedge;
  void requireTransportVectorsAcrossHalfedge();
  void unrequireTransportVectorsAcrossHalfedge();

  // Halfedge vectors in vertex tangent space
  HalfedgeData<Vector2> halfedgeVectorsInVertex;
  void requireHalfedgeVectorsInVertex();
  void unrequireHalfedgeVectorsInVertex();

  // Vertex transport across halfedges
  HalfedgeData<Vector2> transportVectorsAlongHalfedge;
  void requireTransportVectorsAlongHalfedge();
  void unrequireTransportVectorsAlongHalfedge();


  // == Operators

  // Cotan laplacian
  Eigen::SparseMatrix<double> cotanLaplacian;
  void requireCotanLaplacian();
  void unrequireCotanLaplacian();

  // Vertex lumped mass matrix
  Eigen::SparseMatrix<double> vertexLumpedMassMatrix;
  void requireVertexLumpedMassMatrix();
  void unrequireVertexLumpedMassMatrix();

  // Vertex Galerkin Mass Matrix
  Eigen::SparseMatrix<double> vertexGalerkinMassMatrix;
  void requireVertexGalerkinMassMatrix();
  void unrequireVertexGalerkinMassMatrix();

  // Vertex connection Laplacian
  Eigen::SparseMatrix<std::complex<double>> vertexConnectionLaplacian;
  void requireVertexConnectionLaplacian();
  void unrequireVertexConnectionLaplacian();

  // Face Galerkin Mass Matrix
  Eigen::SparseMatrix<double> faceGalerkinMassMatrix;
  void requireFaceGalerkinMassMatrix();
  void unrequireFaceGalerkinMassMatrix();

  // Face connection Laplacian
  Eigen::SparseMatrix<std::complex<double>> faceConnectionLaplacian;
  void requireFaceConnectionLaplacian();
  void unrequireFaceConnectionLaplacian();

  // Crouzeix-Raviart Laplace matrix
  Eigen::SparseMatrix<double> crouzeixRaviartLaplacian;
  void requireCrouzeixRaviartLaplacian();
  void unrequireCrouzeixRaviartLaplacian();

  // Crouzeix-Raviart mass matrix
  Eigen::SparseMatrix<double> crouzeixRaviartMassMatrix;
  void requireCrouzeixRaviartMassMatrix();
  void unrequireCrouzeixRaviartMassMatrix();

  // Crouzeix-Raviart connection Laplacian. Corresponds to the complex version, not the 2E x 2E version.
  Eigen::SparseMatrix<std::complex<double>> crouzeixRaviartConnectionLaplacian;
  void requireCrouzeixRaviartConnectionLaplacian();
  void unrequireCrouzeixRaviartConnectionLaplacian();

  // DEC Operators
  Eigen::SparseMatrix<double> hodge0, hodge0Inverse, hodge1, hodge1Inverse, hodge2, hodge2Inverse, d0, d1;
  void requireDECOperators();
  void unrequireDECOperators();

protected:
  // == Lengths, areas, and angles

  // Edge lengths
  // Note that computeEdgeLengths() is pure virtual: some input data class which extends this interface must supply a
  // method for computing edge lengths (EdgeLengthGeometry serves this purpose)
  DependentQuantityD<EdgeData<double>> edgeLengthsQ;
  virtual void computeEdgeLengths() = 0;

  // Face areas
  DependentQuantityD<FaceData<double>> faceAreasQ;
  virtual void computeFaceAreas();

  // Vertex dual area
  DependentQuantityD<VertexData<double>> vertexDualAreasQ;
  virtual void computeVertexDualAreas();

  // Corner angles
  DependentQuantityD<CornerData<double>> cornerAnglesQ;
  virtual void computeCornerAngles();

  // Vertex angle sums
  DependentQuantityD<VertexData<double>> vertexAngleSumsQ;
  virtual void computeVertexAngleSums();

  // Corner scaled angles
  DependentQuantityD<CornerData<double>> cornerScaledAnglesQ;
  virtual void computeCornerScaledAngles();

  // Vertex gaussian curvature
  DependentQuantityD<VertexData<double>> vertexGaussianCurvaturesQ;
  virtual void computeVertexGaussianCurvatures();

  // Face gaussian curvature
  DependentQuantityD<FaceData<double>> faceGaussianCurvaturesQ;
  virtual void computeFaceGaussianCurvatures();

  // Halfedge cotan weight
  DependentQuantityD<HalfedgeData<double>> halfedgeCotanWeightsQ;
  virtual void computeHalfedgeCotanWeights();

  // Edge cotan weight
  DependentQuantityD<EdgeData<double>> edgeCotanWeightsQ;
  virtual void computeEdgeCotanWeights();

  // Shape length scale
  DependentQuantityD<double> shapeLengthScaleQ;
  virtual void computeShapeLengthScale();

  // Mesh length scale
  DependentQuantityD<double> meshLengthScaleQ;
  virtual void computeMeshLengthScale();


  // == Tangent vectors and transport

  // Halfedge vectors in face
  DependentQuantityD<HalfedgeData<Vector2>> halfedgeVectorsInFaceQ;
  virtual void computeHalfedgeVectorsInFace();

  // Face tangent vector transport across halfedges
  DependentQuantityD<HalfedgeData<Vector2>> transportVectorsAcrossHalfedgeQ;
  virtual void computeTransportVectorsAcrossHalfedge();

  // Halfedge vectors in vertex tangent space
  DependentQuantityD<HalfedgeData<Vector2>> halfedgeVectorsInVertexQ;
  virtual void computeHalfedgeVectorsInVertex();

  // Vertex transport across halfedges
  DependentQuantityD<HalfedgeData<Vector2>> transportVectorsAlongHalfedgeQ;
  virtual void computeTransportVectorsAlongHalfedge();


  // == Operators

  // Cotan laplacian
  DependentQuantityD<Eigen::SparseMatrix<double>> cotanLaplacianQ;
  virtual void computeCotanLaplacian();

  // Vertex lumped mass matrix
  DependentQuantityD<Eigen::SparseMatrix<double>> vertexLumpedMassMatrixQ;
  virtual void computeVertexLumpedMassMatrix();

  // Vertex Galerkin Mass Matrix
  DependentQuantityD<Eigen::SparseMatrix<double>> vertexGalerkinMassMatrixQ;
  virtual void computeVertexGalerkinMassMatrix();

  // Vertex connection Laplacian
  DependentQuantityD<Eigen::SparseMatrix<std::complex<double>>> vertexConnectionLaplacianQ;
  virtual void computeVertexConnectionLaplacian();

  // Face Galerkin Mass Matrix
  DependentQuantityD<Eigen::SparseMatrix<double>> faceGalerkinMassMatrixQ;
  virtual void computeFaceGalerkinMassMatrix();

  // Face connection Laplacian
  DependentQuantityD<Eigen::SparseMatrix<std::complex<double>>> faceConnectionLaplacianQ;
  virtual void computeFaceConnectionLaplacian();

  // Crouzeix-Raviart Laplacian
  DependentQuantityD<Eigen::SparseMatrix<double>> crouzeixRaviartLaplacianQ;
  virtual void computeCrouzeixRaviartLaplacian();

  // Crouzeix-Raviart mass matrix
  DependentQuantityD<Eigen::SparseMatrix<double>> crouzeixRaviartMassMatrixQ;
  virtual void computeCrouzeixRaviartMassMatrix();

  // Crouzeix-Raviart connection Laplacian
  DependentQuantityD<Eigen::SparseMatrix<std::complex<double>>> crouzeixRaviartConnectionLaplacianQ;
  virtual void computeCrouzeixRaviartConnectionLaplacian();

  // DEC Operators
  // Note: The DEC operators deviate from the convention of one member per quantity. This extra array allows the
  // DependentQuantityD<> helper type to still manage and clear out these members.
  std::array<Eigen::SparseMatrix<double>*, 8> DECOperatorArray;
  DependentQuantityD<std::array<Eigen::SparseMatrix<double>*, 8>> DECOperatorsQ;
  virtual void computeDECOperators();
};

} // namespace surface
} // namespace geometrycentral
