#include "geometrycentral/surface/integer_coordinates_intrinsic_triangulation.h"
#include "geometrycentral/surface/intrinsic_triangulation.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/signpost_intrinsic_triangulation.h"
#include "geometrycentral/surface/transfer_functions.h"
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
  for (const MeshAsset& a : {getAsset("fox.ply", true), getAsset("cat_head.obj", true)}) {
    a.printThyName();
    ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
    VertexPositionGeometry& origGeometry = *a.geometry;

    SignpostIntrinsicTriangulation tri(mesh, origGeometry);

    tri.flipToDelaunay();

    EXPECT_TRUE(tri.isDelaunay());
  }
}

TEST_F(IntrinsicTriangulationSuite, IntegerFlip) {
  for (const MeshAsset& a : {getAsset("fox.ply", true), getAsset("cat_head.obj", true)}) {
    a.printThyName();
    ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
    VertexPositionGeometry& origGeometry = *a.geometry;

    IntegerCoordinatesIntrinsicTriangulation tri(mesh, origGeometry);

    tri.flipToDelaunay();

    EXPECT_TRUE(tri.isDelaunay());
  }
}

TEST_F(IntrinsicTriangulationSuite, DelaunayTriangulationsAgree) {
  for (const MeshAsset& a : {getAsset("fox.ply", true), getAsset("cat_head.obj", true)}) {
    a.printThyName();
    ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
    VertexPositionGeometry& origGeometry = *a.geometry;

    IntegerCoordinatesIntrinsicTriangulation tri_int(mesh, origGeometry);
    SignpostIntrinsicTriangulation tri_sign(mesh, origGeometry);

    tri_int.flipToDelaunay();
    tri_sign.flipToDelaunay();

    for (size_t iE = 0; iE < tri_int.intrinsicMesh->nEdges(); iE++) {
      double l_int = tri_int.edgeLengths[tri_int.intrinsicMesh->edge(iE)];
      double l_sign = tri_sign.edgeLengths[tri_sign.intrinsicMesh->edge(iE)];
      EXPECT_NEAR(l_int, l_sign, 1e-5);
    }
  }
}


TEST_F(IntrinsicTriangulationSuite, SignpostTrace) {
  for (const MeshAsset& a : {getAsset("fox.ply", true), getAsset("cat_head.obj", true)}) {
    a.printThyName();
    ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
    VertexPositionGeometry& origGeometry = *a.geometry;

    SignpostIntrinsicTriangulation tri(mesh, origGeometry);

    tri.flipToDelaunay();

    EdgeData<std::vector<SurfacePoint>> out = tri.traceAllIntrinsicEdgesAlongInput();
    for (Edge e : tri.mesh.edges()) {
      EXPECT_GE(out[e].size(), 2);
    }
  }
}

TEST_F(IntrinsicTriangulationSuite, IntegerTrace) {
  for (const MeshAsset& a : {getAsset("fox.ply", true), getAsset("cat_head.obj", true)}) {
    a.printThyName();
    ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
    VertexPositionGeometry& origGeometry = *a.geometry;

    IntegerCoordinatesIntrinsicTriangulation tri(mesh, origGeometry);

    tri.flipToDelaunay();

    EdgeData<std::vector<SurfacePoint>> out = tri.traceAllIntrinsicEdgesAlongInput();
    for (Edge e : tri.mesh.edges()) {
      EXPECT_EQ(out[e].size(), std::max(0, tri.normalCoordinates[e]) + 2);
    }
  }
}

TEST_F(IntrinsicTriangulationSuite, SignpostEquivalentPoint) {
  for (const MeshAsset& a : {getAsset("fox.ply", true), getAsset("cat_head.obj", true)}) {
    a.printThyName();
    ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
    VertexPositionGeometry& origGeometry = *a.geometry;

    SignpostIntrinsicTriangulation tri(mesh, origGeometry);

    tri.flipToDelaunay();

    auto roundtripTest = [&](SurfacePoint origPoint) {
      Vector3 origPos = origPoint.interpolate(origGeometry.vertexPositions);

      // Map to the intrinsic surface
      SurfacePoint intPoint = tri.equivalentPointOnIntrinsic(origPoint);

      // Map back to the input surface
      SurfacePoint returnPoint = tri.equivalentPointOnInput(intPoint);

      // Check that the result is close to where we started
      Vector3 returnPos = returnPoint.interpolate(origGeometry.vertexPositions);

      EXPECT_LT((origPos - returnPos).norm(), 1e-5);
    };


    // Pick a bunch of face points and map them back and forth; verify we get get very similar locations
    std::mt19937 mt(42);
    auto randBary = [&]() {
      std::uniform_real_distribution<double> dist(0.0, 1.0);
      double r1 = dist(mt);
      double r2 = dist(mt);
      Vector3 bary{1 - std::sqrt(r1), std::sqrt(r1) * (1 - r2), std::sqrt(r1) * r2};
      return bary;
    };
    for (Face f : mesh.faces()) {
      SurfacePoint origPoint(f, randBary());
      roundtripTest(origPoint);
    }

    // Test a few vertex points & edge points
    roundtripTest(SurfacePoint(mesh.vertex(12)));
    roundtripTest(SurfacePoint(mesh.vertex(42)));
    roundtripTest(SurfacePoint(mesh.vertex(55)));
    roundtripTest(SurfacePoint(mesh.edge(55), 0.4));
    roundtripTest(SurfacePoint(mesh.edge(55), 0.0));
    roundtripTest(SurfacePoint(mesh.edge(55), 1.0));
    roundtripTest(SurfacePoint(mesh.edge(11), 1.0));
  }
}

TEST_F(IntrinsicTriangulationSuite, IntegerEdgeTraceAgreesWithBulk) {
  for (const MeshAsset& a : {getAsset("fox.ply", true), getAsset("cat_head.obj", true)}) {
    a.printThyName();
    ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
    VertexPositionGeometry& origGeometry = *a.geometry;

    IntegerCoordinatesIntrinsicTriangulation tri(mesh, origGeometry);

    tri.flipToDelaunay();

    // Traced individually
    EdgeData<std::vector<SurfacePoint>> out1(tri.inputMesh);
    for (Edge e : tri.inputMesh.edges()) {
      out1[e] = tri.traceInputHalfedgeAlongIntrinsic(e.halfedge());
    }

    // Traced via common subdivision
    EdgeData<std::vector<SurfacePoint>> out2 = tri.traceAllInputEdgesAlongIntrinsic();

    for (Edge e : tri.inputMesh.edges()) {
      EXPECT_EQ(out1[e].size(), out2[e].size());
      for (size_t iP = 0; iP < out1[e].size(); iP++) {
        EXPECT_EQ(out1[e][iP].type, out2[e][iP].type);
        switch (out1[e][iP].type) {
        case SurfacePointType::Vertex:
          EXPECT_EQ(out1[e][iP].vertex, out2[e][iP].vertex);
          break;
        case SurfacePointType::Edge:
          EXPECT_EQ(out1[e][iP].edge, out2[e][iP].edge);
          EXPECT_NEAR(out1[e][iP].tEdge, out2[e][iP].tEdge, 1e-5);
          break;
        case SurfacePointType::Face:
          EXPECT_EQ(out1[e][iP].face, out2[e][iP].face);
          EXPECT_NEAR((out1[e][iP].faceCoords - out2[e][iP].faceCoords).norm(), 0, 1e-5);
          break;
        }
      }
    }
  }
}

TEST_F(IntrinsicTriangulationSuite, TraceInputEdgeAlongIntrinsicSignpostVsInteger) {
  for (const MeshAsset& a : {getAsset("fox.ply", true), getAsset("cat_head.obj", true)}) {
    a.printThyName();
    ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
    VertexPositionGeometry& origGeometry = *a.geometry;

    SignpostIntrinsicTriangulation tri_signpost(mesh, origGeometry);
    IntegerCoordinatesIntrinsicTriangulation tri_integer(mesh, origGeometry);

    tri_signpost.flipToDelaunay();
    tri_integer.flipToDelaunay();

    for (Halfedge he : mesh.halfedges()) {
      std::vector<SurfacePoint> trace_signpost = tri_signpost.traceInputHalfedgeAlongIntrinsic(he);
      std::vector<SurfacePoint> trace_integer = tri_integer.traceInputHalfedgeAlongIntrinsic(he);

      EXPECT_EQ(trace_signpost.size(), trace_integer.size());
      for (size_t iP = 0; iP < trace_signpost.size(); iP++) {
        SurfacePoint& p_signpost = trace_signpost[iP];
        SurfacePoint& p_integer = trace_integer[iP];
        EXPECT_EQ(p_signpost.type, p_integer.type);
        switch (p_signpost.type) {
        case SurfacePointType::Vertex:
          EXPECT_EQ(p_signpost.vertex, p_integer.vertex);
          break;
        case SurfacePointType::Edge:
          EXPECT_EQ(p_signpost.edge, p_integer.edge);
          EXPECT_NEAR(p_signpost.tEdge, p_integer.tEdge, 1e-5);
          break;
        case SurfacePointType::Face:
          EXPECT_EQ(p_signpost.face, p_integer.face);
          EXPECT_NEAR((p_signpost.faceCoords - p_integer.faceCoords).norm(), 0, 1e-5);
          break;
        }
      }
    }
  }
}


TEST_F(IntrinsicTriangulationSuite, SignpostRefine) {
  for (const MeshAsset& a : {getAsset("fox.ply", true), getAsset("cat_head.obj", true)}) {
    a.printThyName();
    ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
    VertexPositionGeometry& origGeometry = *a.geometry;

    SignpostIntrinsicTriangulation tri(mesh, origGeometry);

    tri.delaunayRefine();
    EXPECT_TRUE(tri.isDelaunay());

    // (technically on some meshes no insertions may be needed, but for this test lets choose meshes that do need it)
    EXPECT_GT(tri.mesh.nVertices(), tri.inputMesh.nVertices());

    // (technically we should check the minimum angle away from needle-like vertices)
    EXPECT_GE(tri.minAngleDegrees(), 25);
  }
}

TEST_F(IntrinsicTriangulationSuite, IntegerRefine) {
  for (const MeshAsset& a : {getAsset("fox.ply", true), getAsset("cat_head.obj", true)}) {
    a.printThyName();
    ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
    VertexPositionGeometry& origGeometry = *a.geometry;

    IntegerCoordinatesIntrinsicTriangulation tri(mesh, origGeometry);

    tri.delaunayRefine();
    EXPECT_TRUE(tri.isDelaunay());

    // (technically on some meshes no insertions may be needed, but for this test lets choose meshes that do need it)
    EXPECT_GT(tri.mesh.nVertices(), tri.inputMesh.nVertices());

    // (technically we should check the minimum angle away from needle-like vertices)
    EXPECT_GE(tri.minAngleDegrees(), 25);
  }
}

TEST_F(IntrinsicTriangulationSuite, SignpostCommonSubdivision) {
  for (const MeshAsset& a : {getAsset("fox.ply", true), getAsset("cat_head.obj", true)}) {
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
  for (const MeshAsset& a : {getAsset("fox.ply", true), getAsset("cat_head.obj", true)}) {
    a.printThyName();
    ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
    VertexPositionGeometry& origGeometry = *a.geometry;

    IntegerCoordinatesIntrinsicTriangulation tri(mesh, origGeometry);

    tri.delaunayRefine();
    CommonSubdivision& cs = tri.getCommonSubdivision();
    bool triangulate = false;
    cs.constructMesh(triangulate);

    EXPECT_GT(cs.mesh->nVertices(), tri.inputMesh.nVertices());
    EXPECT_GT(cs.mesh->nVertices(), tri.mesh.nVertices());

    // Check element counts against a few ways of computing them
    size_t nV, nE, nF;
    std::tie(nV, nE, nF) = cs.elementCounts();

    EXPECT_EQ(cs.mesh->nVertices(), nV);
    EXPECT_EQ(cs.mesh->nEdges(), nE);
    EXPECT_EQ(cs.mesh->nFaces(), nF);

    size_t nV_normal = tri.intrinsicMesh->nVertices();
    for (Edge e : tri.intrinsicMesh->edges()) {
      nV_normal += fmax(0, tri.normalCoordinates[e]);
    }

    EXPECT_EQ(cs.mesh->nVertices(), nV_normal);
  }
}

// TODO test signpost and integer against each other to verify they give same results

TEST_F(IntrinsicTriangulationSuite, CommonSubdivisionCompareIntegerSignpost) {
  for (const MeshAsset& a : {getAsset("fox.ply", true), getAsset("cat_head.obj", true)}) {

    ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
    VertexPositionGeometry& origGeometry = *a.geometry;

    // Construct a common subdivision with both data structures

    IntegerCoordinatesIntrinsicTriangulation tri_int(mesh, origGeometry);
    SignpostIntrinsicTriangulation tri_sign(mesh, origGeometry);

    // Refinement gives a slightly different mesh, which is a little annoying but probably
    // not worth fixing.
    //tri_int.delaunayRefine();
    //tri_sign.delaunayRefine();
    tri_int.flipToDelaunay();
    tri_sign.flipToDelaunay();
    EXPECT_EQ(tri_int.intrinsicMesh->nVertices(), tri_sign.intrinsicMesh->nVertices());

    CommonSubdivision& cs_int = tri_int.getCommonSubdivision();
    CommonSubdivision& cs_sign = tri_sign.getCommonSubdivision();
    cs_int.constructMesh();
    cs_sign.constructMesh();

    // Check that the element counts of the commont subdivision are the same
    EXPECT_EQ(cs_int.mesh->nVertices(), cs_sign.mesh->nVertices());
    EXPECT_EQ(cs_int.mesh->nEdges(), cs_sign.mesh->nEdges());
    EXPECT_EQ(cs_int.mesh->nFaces(), cs_sign.mesh->nFaces());

    // Check that all faces have very similar area
    // (technically we could generate two correct subdivisions with the faces not
    // in correspondence, but in the current implementation they will be, so we'll 
    // just test that)
    origGeometry.requireEdgeLengths();
    EdgeData<double> len_int = cs_int.interpolateEdgeLengthsA(origGeometry.edgeLengths);
    EdgeData<double> len_sign = cs_sign.interpolateEdgeLengthsA(origGeometry.edgeLengths);
    for(size_t iF = 0; iF < cs_int.mesh->nFaces(); iF++) {
      EXPECT_NEAR(len_int[iF], len_sign[iF], 1e-5);
    }
  }
}



// TODO: also test with signposts?
TEST_F(IntrinsicTriangulationSuite, FunctionTransfer) {
  for (const MeshAsset& a : {getAsset("fox.ply", true)}) {
    a.printThyName();
    ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
    VertexPositionGeometry& origGeometry = *a.geometry;

    IntegerCoordinatesIntrinsicTriangulation tri(mesh, origGeometry);

    tri.delaunayRefine();
    CommonSubdivision& cs = tri.getCommonSubdivision();

    AttributeTransfer transfer(cs, origGeometry);
    VertexData<double> data_B(*tri.intrinsicMesh, Vector<double>::Random(tri.intrinsicMesh->nVertices()));
    VertexData<double> data_A_Pointwise = transfer.transferBtoA(data_B, TransferMethod::Pointwise);
    VertexData<double> data_A_L2 = transfer.transferBtoA(data_B, TransferMethod::L2);
    Vector<double> truth = transfer.P_B * data_B.toVector();
    Vector<double> pointwiseA = transfer.P_A * data_A_Pointwise.toVector();
    Vector<double> L2A = transfer.P_A * data_A_L2.toVector();

    double pointwiseErr = (pointwiseA - truth).dot(transfer.M_CS_Galerkin * (pointwiseA - truth));
    double L2Err = (L2A - truth).dot(transfer.M_CS_Galerkin * (L2A - truth));

    EXPECT_LE(L2Err, pointwiseErr);

    SparseMatrix<double> lhs, rhs;
    std::tie(lhs, rhs) = transfer.constructBtoAMatrices();
    Vector<double> residual = lhs * data_A_L2.toVector() - rhs * data_B.toVector();
    EXPECT_LE(residual.norm(), 1e-6);
  }
}

TEST_F(IntrinsicTriangulationSuite, CommonSubdivisionGeometry) {
  for (const MeshAsset& a : {getAsset("fox.ply", true)}) {
    a.printThyName();
    ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
    VertexPositionGeometry& origGeometry = *a.geometry;

    IntegerCoordinatesIntrinsicTriangulation tri(mesh, origGeometry);

    tri.delaunayRefine();
    CommonSubdivision& cs = tri.getCommonSubdivision();
    cs.constructMesh();

    // == Edge lengths
    // Lengths from extrinsic vertex positions
    const VertexData<Vector3>& posCS = cs.interpolateAcrossA(origGeometry.vertexPositions);
    VertexPositionGeometry csGeo(*cs.mesh, posCS);
    csGeo.requireEdgeLengths();
    EdgeData<double> lengthsFromPosA = csGeo.edgeLengths;
    csGeo.unrequireEdgeLengths();

    // Lengths from extrinsic edge lengths
    origGeometry.requireEdgeLengths();
    const EdgeData<double>& lengthsA = origGeometry.edgeLengths;
    EdgeData<double> lengthsFromLenA = cs.interpolateEdgeLengthsA(lengthsA);

    // Lengths from intrinsic edge lengths
    EdgeData<double> lengthsFromLenB = cs.interpolateEdgeLengthsB(tri.edgeLengths);

    EXPECT_NEAR((lengthsFromPosA.toVector() - lengthsFromLenA.toVector()).norm(), 0, 1e-5);
    EXPECT_NEAR((lengthsFromPosA.toVector() - lengthsFromLenB.toVector()).norm(), 0, 1e-5);
    EXPECT_NEAR((lengthsFromLenA.toVector() - lengthsFromLenB.toVector()).norm(), 0, 1e-5);
  }
}
