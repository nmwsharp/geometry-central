#include "geometrycentral/surface/integer_coordinates_intrinsic_triangulation.h"
#include "geometrycentral/surface/intrinsic_triangulation.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/signpost_intrinsic_triangulation.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "load_test_meshes.h"

#include "gtest/gtest.h"

#include <iostream>
#include <string>
#include <unordered_set>


using namespace geometrycentral;
using namespace geometrycentral::surface;

class IntrinsicTriangulationSuite : public MeshAssetSuite {};

TEST_F(IntrinsicTriangulationSuite, SignpostFlip) {
  for (const MeshAsset& a : {getAsset("fox.ply", true)}) {
    a.printThyName();
    ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
    VertexPositionGeometry& origGeometry = *a.geometry;

    SignpostIntrinsicTriangulation tri(mesh, origGeometry);

    tri.flipToDelaunay();
  }
}

TEST_F(IntrinsicTriangulationSuite, IntegerFlip) {
  for (const MeshAsset& a : {getAsset("fox.ply", true)}) {
    a.printThyName();
    ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
    VertexPositionGeometry& origGeometry = *a.geometry;

    IntegerCoordinatesIntrinsicTriangulation tri(mesh, origGeometry);

    tri.flipToDelaunay();
  }
}


TEST_F(IntrinsicTriangulationSuite, SignpostTrace) {
  for (const MeshAsset& a : {getAsset("fox.ply", true)}) {
    a.printThyName();
    ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
    VertexPositionGeometry& origGeometry = *a.geometry;

    SignpostIntrinsicTriangulation tri(mesh, origGeometry);

    tri.flipToDelaunay();

    EdgeData<std::vector<SurfacePoint>> out = tri.traceAllIntrinsicEdgesAlongInput();
    for(Edge e : tri.inputMesh.edges()) {
      EXPECT_GE(out[e].size(), 2);
    }
  }
}

// TODO integer tracing
//TEST_F(IntrinsicTriangulationSuite, IntegerTrace) {
  //for (const MeshAsset& a : {getAsset("fox.ply", true)}) {
    //a.printThyName();
    //ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
    //VertexPositionGeometry& origGeometry = *a.geometry;

    //IntegerCoordinatesIntrinsicTriangulation tri(mesh, origGeometry);

    //tri.flipToDelaunay();

    //EdgeData<std::vector<SurfacePoint>> out = tri.traceAllIntrinsicEdgesAlongInput();
    //for(Edge e : tri.inputMesh.edges()) {
      //EXPECT_GE(out[e].size(), 2);
    //}
  //}
//}


TEST_F(IntrinsicTriangulationSuite, SignpostRefine) {
  for (const MeshAsset& a : {getAsset("fox.ply", true)}) {
    a.printThyName();
    ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
    VertexPositionGeometry& origGeometry = *a.geometry;

    SignpostIntrinsicTriangulation tri(mesh, origGeometry);

    tri.delaunayRefine();

    // (technically on some meshes no insertions may be needed, but for this test lets choose meshes that do need it)
    EXPECT_GT(tri.mesh.nVertices(), tri.inputMesh.nVertices());
  }
}

TEST_F(IntrinsicTriangulationSuite, IntegerRefine) {
  for (const MeshAsset& a : {getAsset("fox.ply", true)}) {
    a.printThyName();
    ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
    VertexPositionGeometry& origGeometry = *a.geometry;

    IntegerCoordinatesIntrinsicTriangulation tri(mesh, origGeometry);

    tri.delaunayRefine();

    // (technically on some meshes no insertions may be needed, but for this test lets choose meshes that do need it)
    EXPECT_GT(tri.mesh.nVertices(), tri.inputMesh.nVertices());
  }
}

TEST_F(IntrinsicTriangulationSuite, SignpostCommonSubdivision) {
  for (const MeshAsset& a : {getAsset("fox.ply", true)}) {
    a.printThyName();
    ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
    VertexPositionGeometry& origGeometry = *a.geometry;

    SignpostIntrinsicTriangulation tri(mesh, origGeometry);

    tri.delaunayRefine();
    CommonSubdivision& cs = tri.getCommonSubdivision();
    cs.constructMesh();
    
    EXPECT_GT(cs.mesh->nVertices(), tri.inputMesh.nVertices());
    EXPECT_GT(cs.mesh->nVertices(), tri.mesh.nVertices());
  }
}

TEST_F(IntrinsicTriangulationSuite, IntegerCommonSubdivision) {
  for (const MeshAsset& a : {getAsset("fox.ply", true)}) {
    a.printThyName();
    ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
    VertexPositionGeometry& origGeometry = *a.geometry;

    IntegerCoordinatesIntrinsicTriangulation tri(mesh, origGeometry);

    tri.delaunayRefine();
    CommonSubdivision& cs = tri.getCommonSubdivision();
    cs.constructMesh();
    
    EXPECT_GT(cs.mesh->nVertices(), tri.inputMesh.nVertices());
    EXPECT_GT(cs.mesh->nVertices(), tri.mesh.nVertices());
  }
}

// TODO test signpost and integer against each other to verify they give same results
