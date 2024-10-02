#include "geometrycentral/surface/embedded_geometry_interface.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "load_test_meshes.h"

#include "gtest/gtest.h"

#include <iostream>
#include <string>

using namespace geometrycentral;
using namespace geometrycentral::surface;

class PolygonMeshSuite : public MeshAssetSuite {};

/* Test that the lumped mass matrices correspond with their unlumped versions. */
TEST_F(PolygonMeshSuite, MassLumpingTest) {

  double epsilon = 1e-8;
  for (MeshAsset& a : allMeshes()) {
    a.printThyName();
    SurfaceMesh& mesh = *a.mesh;
    VertexPositionGeometry& geometry = *a.geometry;

    geometry.requireSimplePolygonVertexLumpedMassMatrix();
    geometry.requireSimplePolygonVertexGalerkinMassMatrix();

    SparseMatrix<double> Mt = geometry.simplePolygonVertexGalerkinMassMatrix.transpose();
    SparseMatrix<double>& M = geometry.simplePolygonVertexLumpedMassMatrix;
    for (size_t i = 0; i < mesh.nVertices(); i++) {
      double rowSum = 0.;
      for (SparseMatrix<double>::InnerIterator it(Mt, i); it; ++it) {
        rowSum += it.value();
      }
      EXPECT_LT(std::abs(rowSum - M.coeffRef(i, i)), epsilon);
    }

    geometry.unrequireSimplePolygonVertexLumpedMassMatrix();
    geometry.unrequireSimplePolygonVertexGalerkinMassMatrix();
  }
}

/* Check that polygon Laplacians and mass matrices coincide with simplicial versions on a triangle mesh. */
TEST_F(PolygonMeshSuite, TriangularTest) {

  double epsilon = 1e-8;
  for (MeshAsset& a : triangularMeshes()) {
    a.printThyName();
    SurfaceMesh& mesh = *a.mesh;
    VertexPositionGeometry& geometry = *a.geometry;

    geometry.requireVertexGalerkinMassMatrix();
    geometry.requireVertexLumpedMassMatrix();
    geometry.requireCotanLaplacian();
    geometry.requireSimplePolygonLaplacian();
    geometry.requireSimplePolygonVertexLumpedMassMatrix();
    geometry.requireSimplePolygonVertexGalerkinMassMatrix();
    geometry.requirePolygonLaplacian();
    geometry.requirePolygonVertexLumpedMassMatrix();

    double L = geometry.cotanLaplacian.norm();
    double Mg = geometry.vertexGalerkinMassMatrix.norm();
    double Ml = geometry.vertexLumpedMassMatrix.norm();

    EXPECT_LT((geometry.simplePolygonLaplacian - geometry.cotanLaplacian).norm() / L, epsilon);
    EXPECT_LT((geometry.polygonLaplacian - geometry.cotanLaplacian).norm() / L, epsilon);
    EXPECT_LT((geometry.simplePolygonVertexGalerkinMassMatrix - geometry.vertexGalerkinMassMatrix).norm() / Mg,
              epsilon);
    EXPECT_LT((geometry.simplePolygonVertexLumpedMassMatrix - geometry.vertexLumpedMassMatrix).norm() / Ml, epsilon);
    EXPECT_LT((geometry.polygonVertexLumpedMassMatrix - geometry.vertexLumpedMassMatrix).norm() / Ml, epsilon);

    geometry.unrequireVertexGalerkinMassMatrix();
    geometry.unrequireVertexLumpedMassMatrix();
    geometry.unrequireCotanLaplacian();
    geometry.unrequireSimplePolygonLaplacian();
    geometry.unrequireSimplePolygonVertexLumpedMassMatrix();
    geometry.unrequireSimplePolygonVertexGalerkinMassMatrix();
    geometry.unrequirePolygonLaplacian();
    geometry.unrequirePolygonVertexLumpedMassMatrix();
  }
}

/* Check that Laplacian can be assembled as expected from DEC operators. */
TEST_F(PolygonMeshSuite, DECTest) {

  double epsilon = 1e-8;
  for (MeshAsset& a : allMeshes()) {
    a.printThyName();
    SurfaceMesh& mesh = *a.mesh;
    VertexPositionGeometry& geometry = *a.geometry;

    geometry.requirePolygonDECOperators();
    geometry.requirePolygonLaplacian();

    SparseMatrix<double>& L = geometry.polygonLaplacian;
    SparseMatrix<double>& h0 = geometry.polygonHodge0;
    SparseMatrix<double>& h0Inv = geometry.polygonHodge0Inverse;
    SparseMatrix<double>& h1 = geometry.polygonHodge1;
    SparseMatrix<double>& h2 = geometry.polygonHodge2;
    SparseMatrix<double>& h2Inv = geometry.polygonHodge2Inverse;
    SparseMatrix<double>& d0 = geometry.polygonD0;
    SparseMatrix<double>& d1 = geometry.polygonD1;
    EXPECT_LT((L - d0.transpose() * h1 * d0).norm(), epsilon);

    if (mesh.isTriangular()) {
      geometry.requireDECOperators();
      EXPECT_LT((h0 - geometry.hodge0).norm(), epsilon);
      EXPECT_LT((h2 - geometry.hodge2).norm(), epsilon);
      geometry.unrequireDECOperators();
    }

    geometry.unrequirePolygonDECOperators();
    geometry.unrequirePolygonLaplacian();
  }
}