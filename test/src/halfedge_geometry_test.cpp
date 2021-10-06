
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"

#include "geometrycentral/surface/base_geometry_interface.h"
#include "geometrycentral/surface/edge_length_geometry.h"
#include "geometrycentral/surface/embedded_geometry_interface.h"
#include "geometrycentral/surface/extrinsic_geometry_interface.h"
#include "geometrycentral/surface/geometry.h"
#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/rich_surface_mesh_data.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "load_test_meshes.h"

#include "gtest/gtest.h"

#include <iostream>
#include <string>
#include <unordered_set>


using namespace geometrycentral;
using namespace geometrycentral::surface;
using std::cout;
using std::endl;

class HalfedgeGeometrySuite : public MeshAssetSuite {};


// ============================================================
// =============== Misc test
// ============================================================

// TODO needs to move to utilities tests
TEST_F(HalfedgeGeometrySuite, VectorPod) {
  // Hijacking
  EXPECT_TRUE(std::is_pod<Vector2>::value);
  EXPECT_TRUE(std::is_pod<Vector3>::value);
}

// ============================================================
// =============== Types
// ============================================================

// Make sure we can construct a pointer to all the geometries
TEST_F(HalfedgeGeometrySuite, GeometryPointers) {

  { std::unique_ptr<BaseGeometryInterface> deleteGeom(getAsset("bob_small.ply", true).geometry.release()); }
  { std::unique_ptr<IntrinsicGeometryInterface> deleteGeom(getAsset("bob_small.ply", true).geometry.release()); }
  { std::unique_ptr<ExtrinsicGeometryInterface> deleteGeom(getAsset("bob_small.ply", true).geometry.release()); }
  { std::unique_ptr<EmbeddedGeometryInterface> deleteGeom(getAsset("bob_small.ply", true).geometry.release()); }
}

// ============================================================
// =============== Quantity management tests
// ============================================================

TEST_F(HalfedgeGeometrySuite, RefreshMutationTest) {
  // for (auto& asset : {getAsset("lego.ply", false), getAsset("lego.ply", true)}) {
  for (auto& asset : {getAsset("lego.ply", true)}) {
    ManifoldSurfaceMesh& mesh = *asset.manifoldMesh;
    VertexPositionGeometry& origGeometry = *asset.geometry;

    // Initial element counts
    size_t nVertexOrig = mesh.nVertices();
    size_t nHalfedgeOrig = mesh.nHalfedges();
    size_t nCornerOrig = mesh.nCorners();
    size_t nEdgeOrig = mesh.nEdges();
    size_t nFaceOrig = mesh.nFaces();

    // Require some quantities
    origGeometry.requireVertexGaussianCurvatures();
    origGeometry.requireVertexPositions();
    origGeometry.requireFaceAreas();

    // Check that quantities are valid
    double EPS = 1e-4;
    for (Vertex v : mesh.vertices()) {
      ASSERT_LT(std::abs(origGeometry.vertexGaussianCurvatures[v]), EPS);
      ASSERT_TRUE(isfinite(origGeometry.vertexPositions[v]));
    }
    for (Face f : mesh.faces()) {
      ASSERT_TRUE(std::isfinite(origGeometry.faceAreas[f]));
    }

    double totalAreaBefore = 0.;
    for (Face f : mesh.faces()) {
      totalAreaBefore += origGeometry.faceAreas[f];
    }

    // Split some edges
    std::vector<Edge> origEdges;
    for (Edge e : mesh.edges()) {
      origEdges.push_back(e);
    }
    for (Edge e : origEdges) {
      Vertex vA = e.halfedge().vertex();
      Vertex vB = e.halfedge().twin().vertex();
      Halfedge he = mesh.splitEdgeTriangular(e);
      origGeometry.vertexPositions[he.vertex()] =
          0.5 * (origGeometry.vertexPositions[vA] + origGeometry.vertexPositions[vB]);
    }
    mesh.validateConnectivity();
    for (Face f : mesh.faces()) {
      mesh.triangulate(f);
    }
    mesh.validateConnectivity();

    // Be sure the mesh actually got bigger
    EXPECT_LT(nVertexOrig, mesh.nVertices());
    EXPECT_LT(nHalfedgeOrig, mesh.nHalfedges());
    EXPECT_LT(nCornerOrig, mesh.nCorners());
    EXPECT_LT(nEdgeOrig, mesh.nEdges());
    EXPECT_LT(nFaceOrig, mesh.nFaces());

    // All important refresh
    origGeometry.refreshQuantities();

    // Check that our quantities are still valid after
    for (Vertex v : mesh.vertices()) {
      ASSERT_LT(std::abs(origGeometry.vertexGaussianCurvatures[v]), EPS);
      ASSERT_TRUE(isfinite(origGeometry.vertexPositions[v]));
    }

    for (Face f : mesh.faces()) {
      ASSERT_TRUE(std::isfinite(origGeometry.faceAreas[f]));
    }
    double totalAreaAfter = 0.;
    for (Face f : mesh.faces()) {
      totalAreaAfter += origGeometry.faceAreas[f];
    }
    ASSERT_NEAR(totalAreaBefore, totalAreaAfter, EPS);
  }
}


TEST_F(HalfedgeGeometrySuite, Purge) {
  for (auto& asset : {getAsset("bob_small.ply", false), getAsset("bob_small.ply", true)}) {
    SurfaceMesh& mesh = *asset.mesh;
    BaseGeometryInterface& geometry = *asset.geometry;

    // Make sure the size is zero when empty
    EXPECT_EQ(geometry.vertexIndices.size(), 0);

    // Get them indices
    geometry.requireVertexIndices();
    EXPECT_EQ(geometry.vertexIndices.size(), mesh.nVertices());

    // Unrequire (but should not get rid of yet)
    geometry.unrequireVertexIndices();
    EXPECT_EQ(geometry.vertexIndices.size(), mesh.nVertices());

    // Purge actually deletes
    geometry.purgeQuantities();
    EXPECT_EQ(geometry.vertexIndices.size(), 0);
  }
}


// The DEC operators use a special array to ensure they all get deleted, make sure it works
TEST_F(HalfedgeGeometrySuite, PurgeTestDEC) {
  auto asset = getAsset("bob_small.ply", false);
  SurfaceMesh& mesh = *asset.mesh;
  IntrinsicGeometryInterface& geometry = *asset.geometry;

  // Make sure the size is zero when empty
  EXPECT_EQ(geometry.d0.nonZeros(), 0);

  // Populate
  geometry.requireDECOperators();
  EXPECT_GT(geometry.d0.nonZeros(), 0);

  // Unrequire (but should not get rid of yet)
  geometry.unrequireDECOperators();
  EXPECT_GT(geometry.d0.nonZeros(), 0);

  // Purge actually deletes
  geometry.purgeQuantities();
  EXPECT_EQ(geometry.d0.nonZeros(), 0);
}


// Copying
TEST_F(HalfedgeGeometrySuite, CopyTest) {
  for (auto& asset : {getAsset("bob_small.ply", false), getAsset("bob_small.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;

    VertexPositionGeometry& geometry = *asset.geometry;

    /* This SHOULD NOT compile, there is no implicit copy constructor
    auto testF = [](VertexPositionGeometry g) {
      g.requireVertexIndices();
    };
    testF(geometry);
    */

    // Copy vertex position
    std::unique_ptr<VertexPositionGeometry> copy1 = geometry.copy();
    copy1->requireFaceAreas();

    // Construct edge length geometry
    copy1->requireEdgeLengths();
    EdgeLengthGeometry eGeom(mesh, copy1->edgeLengths);

    // Copy vertex position
    std::unique_ptr<EdgeLengthGeometry> copy2 = eGeom.copy();
    copy2->requireFaceAreas();
  }
}


// ============================================================
// =============== Constructor tests
// ============================================================

TEST_F(HalfedgeGeometrySuite, PositionMatrixConstructionTest) {
  for (auto& asset : {getAsset("bob_small.ply", false)}) {
    SurfaceMesh& mesh = *asset.mesh;
    VertexPositionGeometry& vGeom = *asset.geometry;

    // Build out a matrix
    DenseMatrix<double> vPos(mesh.nVertices(), 3);
    size_t iV = 0;
    for (Vertex v : mesh.vertices()) {
      Vector3 p = vGeom.vertexPositions[v];
      vPos(iV, 0) = p.x;
      vPos(iV, 1) = p.y;
      vPos(iV, 2) = p.z;
      iV++;
    }

    // Construct a geometry
    VertexPositionGeometry newGeom(mesh, vPos);

    vGeom.requireMeshLengthScale();
    newGeom.requireMeshLengthScale();

    EXPECT_EQ(vGeom.meshLengthScale, newGeom.meshLengthScale);
  }
}

TEST_F(HalfedgeGeometrySuite, MatrixFactoryTest) {
  for (auto& asset : {getAsset("bob_small.ply", false)}) {
    SurfaceMesh& mesh = *asset.mesh;
    VertexPositionGeometry& vGeom = *asset.geometry;

    // Get face mat
    DenseMatrix<size_t> F = mesh.getFaceVertexMatrix<size_t>();

    // Build out a matrix
    DenseMatrix<double> vPos(mesh.nVertices(), 3);
    size_t iV = 0;
    for (Vertex v : mesh.vertices()) {
      Vector3 p = vGeom.vertexPositions[v];
      vPos(iV, 0) = p.x;
      vPos(iV, 1) = p.y;
      vPos(iV, 2) = p.z;
      iV++;
    }

    vGeom.requireMeshLengthScale();

    // Invoke factory

    { // general
      std::unique_ptr<SurfaceMesh> newMesh;
      std::unique_ptr<VertexPositionGeometry> newGeom;
      std::tie(newMesh, newGeom) = makeSurfaceMeshAndGeometry(vPos, F);

      newGeom->requireMeshLengthScale();
      EXPECT_EQ(vGeom.meshLengthScale, newGeom->meshLengthScale);
    }

    { // general
      std::unique_ptr<ManifoldSurfaceMesh> newMesh;
      std::unique_ptr<VertexPositionGeometry> newGeom;
      std::tie(newMesh, newGeom) = makeManifoldSurfaceMeshAndGeometry(vPos, F);

      newGeom->requireMeshLengthScale();
      EXPECT_EQ(vGeom.meshLengthScale, newGeom->meshLengthScale);
    }
  }
}


// ============================================================
// =============== Quantity tests
// ============================================================

// Simple tests which ensure that the quantity can be computed and is a reasonable range.

// == Basic indices

TEST_F(HalfedgeGeometrySuite, VertexIndices) {
  for (auto& asset : {getAsset("bob_small.ply", false), getAsset("bob_small.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    BaseGeometryInterface& geometry = *asset.geometry;

    geometry.requireVertexIndices();
    for (Vertex v : mesh.vertices()) {
      EXPECT_GE(geometry.vertexIndices[v], 0);
      EXPECT_LT(geometry.vertexIndices[v], mesh.nVertices());
    }
  }
}

TEST_F(HalfedgeGeometrySuite, HalfedgeIndices) {
  for (auto& asset : {getAsset("bob_small.ply", false), getAsset("bob_small.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    BaseGeometryInterface& geometry = *asset.geometry;

    geometry.requireHalfedgeIndices();
    for (Halfedge e : mesh.halfedges()) {
      EXPECT_GE(geometry.halfedgeIndices[e], 0);
      EXPECT_LT(geometry.halfedgeIndices[e], mesh.nHalfedges());
    }
  }
}

TEST_F(HalfedgeGeometrySuite, CornerIndices) {
  for (auto& asset : {getAsset("bob_small.ply", false), getAsset("bob_small.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    BaseGeometryInterface& geometry = *asset.geometry;

    geometry.requireCornerIndices();
    for (Corner e : mesh.corners()) {
      EXPECT_GE(geometry.cornerIndices[e], 0);
      EXPECT_LT(geometry.cornerIndices[e], mesh.nCorners());
    }
  }
}

TEST_F(HalfedgeGeometrySuite, EdgeIndices) {
  for (auto& asset : {getAsset("bob_small.ply", false), getAsset("bob_small.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    BaseGeometryInterface& geometry = *asset.geometry;

    geometry.requireEdgeIndices();
    for (Edge e : mesh.edges()) {
      EXPECT_GE(geometry.edgeIndices[e], 0);
      EXPECT_LT(geometry.edgeIndices[e], mesh.nEdges());
    }
  }
}

TEST_F(HalfedgeGeometrySuite, FaceIndices) {
  for (auto& asset : {getAsset("bob_small.ply", false), getAsset("bob_small.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    BaseGeometryInterface& geometry = *asset.geometry;

    geometry.requireFaceIndices();
    for (Face e : mesh.faces()) {
      EXPECT_GE(geometry.faceIndices[e], 0);
      EXPECT_LT(geometry.faceIndices[e], mesh.nFaces());
    }
  }
}

TEST_F(HalfedgeGeometrySuite, BoundaryLoopIndices) {
  for (auto& asset : {getAsset("bob_small.ply", false), getAsset("bob_small.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    BaseGeometryInterface& geometry = *asset.geometry;

    geometry.requireBoundaryLoopIndices();
    for (BoundaryLoop e : mesh.boundaryLoops()) {
      EXPECT_GE(geometry.boundaryLoopIndices[e], 0);
      EXPECT_LT(geometry.boundaryLoopIndices[e], mesh.nBoundaryLoops());
    }
  }
}

// == Intrinsic geometry

TEST_F(HalfedgeGeometrySuite, EdgeLengths) {
  for (auto& asset : {getAsset("bob_small.ply", false), getAsset("bob_small.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    IntrinsicGeometryInterface& geometry = *asset.geometry;

    geometry.requireEdgeLengths();
    for (Edge e : mesh.edges()) {
      EXPECT_GT(geometry.edgeLengths[e], 0);
      EXPECT_TRUE(std::isfinite(geometry.edgeLengths[e]));
    }
  }
}


// Test immediate edge length computation from position
TEST_F(HalfedgeGeometrySuite, EdgeLengthImmediate_Position) {
  for (auto& asset : {getAsset("bob_small.ply", false), getAsset("bob_small.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    IntrinsicGeometryInterface& origGeometry = *asset.geometry;

    // Construct position geometry
    VertexPositionGeometry& geometry = *asset.geometry;

    EdgeData<double> edgeLengthImmediate(mesh);
    for (Edge e : mesh.edges()) {
      edgeLengthImmediate[e] = geometry.edgeLength(e);
    }

    geometry.requireEdgeLengths();
    for (Edge e : mesh.edges()) {
      EXPECT_NEAR(geometry.edgeLengths[e], edgeLengthImmediate[e], 1e-6);
    }
  }
}

// Ensure overrides for computing edge length give the same result
TEST_F(HalfedgeGeometrySuite, EdgeLengthOverrides) {
  for (auto& asset : {getAsset("bob_small.ply", false), getAsset("bob_small.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    VertexPositionGeometry& origGeometry = *asset.geometry;

    // Construct edge length geometry
    origGeometry.requireEdgeLengths();
    EdgeLengthGeometry geometry(mesh, origGeometry.edgeLengths);

    geometry.requireEdgeLengths();
    origGeometry.requireEdgeLengths();
    for (Edge e : mesh.edges()) {
      EXPECT_NEAR(geometry.edgeLengths[e], origGeometry.edgeLengths[e], 1e-6);
    }
  }
}


TEST_F(HalfedgeGeometrySuite, FaceAreas) {
  for (auto& asset : {getAsset("bob_small.ply", false), getAsset("bob_small.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    IntrinsicGeometryInterface& geometry = *asset.geometry;

    geometry.requireFaceAreas();
    for (Face e : mesh.faces()) {
      EXPECT_GE(geometry.faceAreas[e], 0);
      EXPECT_TRUE(std::isfinite(geometry.faceAreas[e]));
    }
  }
}

// Ensure overrides for computing face area give the same result
TEST_F(HalfedgeGeometrySuite, FaceAreaOverrides) {
  for (auto& asset : {getAsset("lego.ply", false), getAsset("lego.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    VertexPositionGeometry& origGeometry = *asset.geometry;

    // Construct edge length geometry
    origGeometry.requireEdgeLengths();
    EdgeLengthGeometry geometry(mesh, origGeometry.edgeLengths);

    geometry.requireFaceAreas();
    origGeometry.requireFaceAreas();
    bool allExactSame = true;
    for (Face f : mesh.faces()) {
      EXPECT_NEAR(geometry.faceAreas[f], origGeometry.faceAreas[f], 1e-6);

      allExactSame = allExactSame && (geometry.faceAreas[f] == origGeometry.faceAreas[f]);
    }

    // Ensure that not all of the values are exactly the same to the bit-- this tell us that the type system is behaving
    // nicely and the override version is actually being invoked.
    EXPECT_FALSE(allExactSame);
  }
}

// Test the immediate face area computation from the length geometry
TEST_F(HalfedgeGeometrySuite, FaceAreasImmediate_Length) {
  for (auto& asset : {getAsset("bob_small.ply", false), getAsset("bob_small.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    IntrinsicGeometryInterface& origGeometry = *asset.geometry;

    // Construct edge length geometry
    origGeometry.requireEdgeLengths();
    EdgeLengthGeometry geometry(mesh, origGeometry.edgeLengths);


    FaceData<double> faceAreaImmediate(mesh);
    for (Face e : mesh.faces()) {
      faceAreaImmediate[e] = geometry.faceArea(e);
    }

    geometry.requireFaceAreas();
    for (Face e : mesh.faces()) {
      EXPECT_NEAR(geometry.faceAreas[e], faceAreaImmediate[e], 1e-6);
    }
  }
}

// Test the immediate face area computation from the position geometry
TEST_F(HalfedgeGeometrySuite, FaceAreasImmediate_Position) {
  for (auto& asset : {getAsset("bob_small.ply", false), getAsset("bob_small.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    IntrinsicGeometryInterface& origGeometry = *asset.geometry;

    // Construct position geometry
    VertexPositionGeometry& geometry = *asset.geometry;

    FaceData<double> faceAreaImmediate(mesh);
    for (Face e : mesh.faces()) {
      faceAreaImmediate[e] = geometry.faceArea(e);
    }

    geometry.requireFaceAreas();
    for (Face e : mesh.faces()) {
      EXPECT_NEAR(geometry.faceAreas[e], faceAreaImmediate[e], 1e-6);
    }
  }
}

TEST_F(HalfedgeGeometrySuite, VertexDualAreas) {
  for (auto& asset : {getAsset("bob_small.ply", false), getAsset("bob_small.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    IntrinsicGeometryInterface& geometry = *asset.geometry;

    geometry.requireVertexDualAreas();
    for (Vertex e : mesh.vertices()) {
      EXPECT_GE(geometry.vertexDualAreas[e], 0);
      EXPECT_TRUE(std::isfinite(geometry.vertexDualAreas[e]));
    }
  }
}

TEST_F(HalfedgeGeometrySuite, CornerAngles) {
  for (auto& asset : {getAsset("bob_small.ply", false), getAsset("bob_small.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    IntrinsicGeometryInterface& geometry = *asset.geometry;

    geometry.requireCornerAngles();
    for (Corner e : mesh.corners()) {
      EXPECT_GE(geometry.cornerAngles[e], 0);
      EXPECT_TRUE(std::isfinite(geometry.cornerAngles[e]));
    }
  }
}

// Test immediate corner angles from length geometry
TEST_F(HalfedgeGeometrySuite, CornerAnglesImmediate_Length) {

  for (auto& asset : {getAsset("bob_small.ply", false), getAsset("bob_small.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    IntrinsicGeometryInterface& origGeometry = *asset.geometry;

    // Construct edge length geometry
    origGeometry.requireEdgeLengths();
    EdgeLengthGeometry geometry(mesh, origGeometry.edgeLengths);

    CornerData<double> cornerAnglesImmediate(mesh);
    for (Corner e : mesh.corners()) {
      cornerAnglesImmediate[e] = geometry.cornerAngle(e);
    }

    geometry.requireCornerAngles();
    for (Corner e : mesh.corners()) {
      EXPECT_NEAR(geometry.cornerAngles[e], cornerAnglesImmediate[e], 1e-6);
    }
  }
}

// Test immediate corner angles from position geometry
TEST_F(HalfedgeGeometrySuite, CornerAnglesImmediate_Position) {
  for (auto& asset : {getAsset("bob_small.ply", false), getAsset("bob_small.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    IntrinsicGeometryInterface& origGeometry = *asset.geometry;

    // Construct position geometry
    VertexPositionGeometry& geometry = *asset.geometry;

    CornerData<double> cornerAnglesImmediate(mesh);
    for (Corner e : mesh.corners()) {
      cornerAnglesImmediate[e] = geometry.cornerAngle(e);
    }

    geometry.requireCornerAngles();
    for (Corner e : mesh.corners()) {
      EXPECT_NEAR(geometry.cornerAngles[e], cornerAnglesImmediate[e], 1e-6);
    }
  }
}

// Ensure overrides for computing corner angles give the same result
TEST_F(HalfedgeGeometrySuite, CornerAngleOverrides) {
  for (auto& asset : {getAsset("lego.ply", false), getAsset("lego.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    VertexPositionGeometry& origGeometry = *asset.geometry;

    // Construct edge length geometry
    origGeometry.requireEdgeLengths();
    EdgeLengthGeometry geometry(mesh, origGeometry.edgeLengths);

    geometry.requireCornerAngles();
    origGeometry.requireCornerAngles();
    bool allExactSame = true;
    for (Corner c : mesh.corners()) {
      EXPECT_NEAR(geometry.cornerAngles[c], origGeometry.cornerAngles[c], 1e-6);

      allExactSame = allExactSame && (geometry.cornerAngles[c] == origGeometry.cornerAngles[c]);
    }

    // Ensure that not all of the values are exactly the same to the bit-- this tell us that the type system is behaving
    // nicely and the override version is actually being invoked.
    EXPECT_FALSE(allExactSame);
  }
}


TEST_F(HalfedgeGeometrySuite, VertexAngleSums) {
  for (auto& asset : {getAsset("bob_small.ply", false), getAsset("bob_small.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    IntrinsicGeometryInterface& geometry = *asset.geometry;

    geometry.requireVertexAngleSums();
    for (Vertex e : mesh.vertices()) {
      EXPECT_GE(geometry.vertexAngleSums[e], 0);
      EXPECT_TRUE(std::isfinite(geometry.vertexAngleSums[e]));
    }
  }
}


TEST_F(HalfedgeGeometrySuite, CornerScaledAngles) {
  for (auto& asset : {getAsset("bob_small.ply", false), getAsset("bob_small.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    IntrinsicGeometryInterface& geometry = *asset.geometry;

    geometry.requireCornerScaledAngles();
    for (Corner e : mesh.corners()) {
      EXPECT_GE(geometry.cornerScaledAngles[e], 0);
      EXPECT_TRUE(std::isfinite(geometry.cornerScaledAngles[e]));
    }
  }
}


TEST_F(HalfedgeGeometrySuite, VertexGaussianCurvatures) {
  for (auto& asset : {getAsset("bob_small.ply", false), getAsset("bob_small.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    IntrinsicGeometryInterface& geometry = *asset.geometry;

    geometry.requireVertexGaussianCurvatures();
    for (Vertex e : mesh.vertices()) {
      EXPECT_TRUE(std::isfinite(geometry.vertexGaussianCurvatures[e]));
      // EXPECT_LT(std::abs(geometry.vertexGaussianCurvatures[e]), 100.);
    }
  }
}

TEST_F(HalfedgeGeometrySuite, FaceGaussianCurvatures) {
  for (auto& asset : {getAsset("bob_small.ply", false), getAsset("bob_small.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    IntrinsicGeometryInterface& geometry = *asset.geometry;

    geometry.requireFaceGaussianCurvatures();
    for (Face e : mesh.faces()) {
      EXPECT_TRUE(std::isfinite(geometry.faceGaussianCurvatures[e]));
    }
  }
}

TEST_F(HalfedgeGeometrySuite, HalfedgeCotanWeights) {
  for (auto& asset : {getAsset("lego.ply", false), getAsset("lego.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    IntrinsicGeometryInterface& geometry = *asset.geometry;

    geometry.requireHalfedgeCotanWeights();
    for (Halfedge e : mesh.halfedges()) {
      EXPECT_TRUE(std::isfinite(geometry.halfedgeCotanWeights[e]));
    }
  }
}

// Test immediate halfedge cotan computation from length
TEST_F(HalfedgeGeometrySuite, HalfedgeCotanWeightsImmediate_Length) {
  for (auto& asset : {getAsset("lego.ply", false), getAsset("lego.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    IntrinsicGeometryInterface& origGeometry = *asset.geometry;

    // Construct edge length geometry
    origGeometry.requireEdgeLengths();
    EdgeLengthGeometry geometry(mesh, origGeometry.edgeLengths);

    HalfedgeData<double> halfedgeCotanWeightsImmediate(mesh);
    for (Halfedge e : mesh.halfedges()) {
      halfedgeCotanWeightsImmediate[e] = geometry.halfedgeCotanWeight(e);
    }

    geometry.requireHalfedgeCotanWeights();
    for (Halfedge e : mesh.halfedges()) {
      EXPECT_NEAR(geometry.halfedgeCotanWeights[e], halfedgeCotanWeightsImmediate[e], 1e-6);
    }
  }
}

// Test immediate halfedge cotan computation from position
TEST_F(HalfedgeGeometrySuite, HalfedgeCotanWeightsImmediate_Position) {
  for (auto& asset : {getAsset("lego.ply", false), getAsset("lego.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    IntrinsicGeometryInterface& origGeometry = *asset.geometry;

    // Construct position geometry
    VertexPositionGeometry& geometry = *asset.geometry;

    HalfedgeData<double> halfedgeCotanWeightsImmediate(mesh);
    for (Halfedge e : mesh.halfedges()) {
      halfedgeCotanWeightsImmediate[e] = geometry.halfedgeCotanWeight(e);
    }

    geometry.requireHalfedgeCotanWeights();
    for (Halfedge e : mesh.halfedges()) {
      EXPECT_NEAR(geometry.halfedgeCotanWeights[e], halfedgeCotanWeightsImmediate[e], 1e-6);
    }
  }
}

// Ensure overrides for computing halfedge cotan weights give the same result
TEST_F(HalfedgeGeometrySuite, HalfedgeCotanWeightOverrides) {
  for (auto& asset : {getAsset("lego.ply", false), getAsset("lego.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    VertexPositionGeometry& origGeometry = *asset.geometry;

    // Construct edge length geometry
    origGeometry.requireEdgeLengths();
    EdgeLengthGeometry geometry(mesh, origGeometry.edgeLengths);

    geometry.requireHalfedgeCotanWeights();
    origGeometry.requireHalfedgeCotanWeights();
    bool allExactSame = true;
    for (Halfedge he : mesh.halfedges()) {
      EXPECT_NEAR(geometry.halfedgeCotanWeights[he], origGeometry.halfedgeCotanWeights[he], 1e-6);

      allExactSame = allExactSame && (geometry.halfedgeCotanWeights[he] == origGeometry.halfedgeCotanWeights[he]);
    }

    // Ensure that not all of the values are exactly the same to the bit-- this tell us that the type system is behaving
    // nicely and the override version is actually being invoked.
    EXPECT_FALSE(allExactSame);
  }
}

TEST_F(HalfedgeGeometrySuite, EdgeCotanWeights) {
  for (auto& asset : {getAsset("lego.ply", false), getAsset("lego.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    IntrinsicGeometryInterface& geometry = *asset.geometry;

    geometry.requireEdgeCotanWeights();
    for (Edge e : mesh.edges()) {
      EXPECT_TRUE(std::isfinite(geometry.edgeCotanWeights[e]));
    }
  }
}

// Test immediate cotan computation from length
TEST_F(HalfedgeGeometrySuite, EdgeCotanWeightsImmediate_Length) {
  for (auto& asset : {getAsset("lego.ply", false), getAsset("lego.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    IntrinsicGeometryInterface& origGeometry = *asset.geometry;

    // Construct edge length geometry
    origGeometry.requireEdgeLengths();
    EdgeLengthGeometry geometry(mesh, origGeometry.edgeLengths);

    EdgeData<double> edgeCotanWeightsImmediate(mesh);
    for (Edge e : mesh.edges()) {
      edgeCotanWeightsImmediate[e] = geometry.edgeCotanWeight(e);
    }

    geometry.requireEdgeCotanWeights();
    for (Edge e : mesh.edges()) {
      EXPECT_NEAR(geometry.edgeCotanWeights[e], edgeCotanWeightsImmediate[e], 1e-6);
    }
  }
}

// Test immediate cotan computation from position
TEST_F(HalfedgeGeometrySuite, EdgeCotanWeightsImmediate_Position) {
  for (auto& asset : {getAsset("lego.ply", false), getAsset("lego.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    IntrinsicGeometryInterface& origGeometry = *asset.geometry;

    // Construct position geometry
    VertexPositionGeometry& geometry = *asset.geometry;

    EdgeData<double> edgeCotanWeightsImmediate(mesh);
    for (Edge e : mesh.edges()) {
      edgeCotanWeightsImmediate[e] = geometry.edgeCotanWeight(e);
    }

    geometry.requireEdgeCotanWeights();
    for (Edge e : mesh.edges()) {
      EXPECT_NEAR(geometry.edgeCotanWeights[e], edgeCotanWeightsImmediate[e], 1e-6);
    }
  }
}

// Ensure overrides for computing edge cotan weights give the same result
TEST_F(HalfedgeGeometrySuite, EdgeCotanWeightOverrides) {
  for (auto& asset : {getAsset("lego.ply", false), getAsset("lego.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    VertexPositionGeometry& origGeometry = *asset.geometry;

    // Construct edge length geometry
    origGeometry.requireEdgeLengths();
    EdgeLengthGeometry geometry(mesh, origGeometry.edgeLengths);

    geometry.requireEdgeCotanWeights();
    origGeometry.requireEdgeCotanWeights();
    bool allExactSame = true;
    for (Edge e : mesh.edges()) {
      EXPECT_NEAR(geometry.edgeCotanWeights[e], origGeometry.edgeCotanWeights[e], 1e-6);

      allExactSame = allExactSame && (geometry.edgeCotanWeights[e] == origGeometry.edgeCotanWeights[e]);
    }

    // Ensure that not all of the values are exactly the same to the bit-- this tell us that the type system is behaving
    // nicely and the override version is actually being invoked.
    EXPECT_FALSE(allExactSame);
  }
}

TEST_F(HalfedgeGeometrySuite, HalfedgeVectorsInFace) {
  for (auto& asset : {getAsset("lego.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    IntrinsicGeometryInterface& geometry = *asset.geometry;

    geometry.requireHalfedgeVectorsInFace();
    for (Halfedge he : mesh.halfedges()) {
      if (he.isInterior()) {
        EXPECT_TRUE(isfinite(geometry.halfedgeVectorsInFace[he]));
      } else {
        EXPECT_FALSE(isfinite(geometry.halfedgeVectorsInFace[he]));
      }
    }
  }
}

TEST_F(HalfedgeGeometrySuite, TransportVectorsAcrossHalfedge) {
  for (auto& asset : {getAsset("lego.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    IntrinsicGeometryInterface& geometry = *asset.geometry;

    geometry.requireTransportVectorsAcrossHalfedge();
    for (Halfedge he : mesh.halfedges()) {
      if (he.edge().isBoundary()) {
        EXPECT_FALSE(isfinite(geometry.transportVectorsAcrossHalfedge[he]));
      } else {
        EXPECT_TRUE(isfinite(geometry.transportVectorsAcrossHalfedge[he]));
      }
    }
  }
}

TEST_F(HalfedgeGeometrySuite, HalfedgeVectorsInVertex) {
  for (auto& asset : {getAsset("lego.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    IntrinsicGeometryInterface& geometry = *asset.geometry;

    geometry.requireHalfedgeVectorsInVertex();
    for (Halfedge he : mesh.halfedges()) {
      EXPECT_TRUE(isfinite(geometry.halfedgeVectorsInVertex[he]));
    }
  }
}

TEST_F(HalfedgeGeometrySuite, TransportVectorsAlongHalfedge) {
  for (auto& asset : {getAsset("lego.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    IntrinsicGeometryInterface& geometry = *asset.geometry;

    geometry.requireTransportVectorsAlongHalfedge();
    for (Halfedge he : mesh.halfedges()) {
      EXPECT_TRUE(isfinite(geometry.transportVectorsAlongHalfedge[he]));
    }
  }
}


TEST_F(HalfedgeGeometrySuite, CotanLaplacian) {
  for (auto& asset : {getAsset("lego.ply", false), getAsset("lego.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    IntrinsicGeometryInterface& geometry = *asset.geometry;

    geometry.requireCotanLaplacian();

    EXPECT_EQ(geometry.cotanLaplacian.rows(), (long int)mesh.nVertices());
    EXPECT_EQ(geometry.cotanLaplacian.cols(), (long int)mesh.nVertices());
    EXPECT_EQ(geometry.cotanLaplacian.nonZeros(), (long int)(mesh.nVertices() + 2 * mesh.nEdges()));

    EXPECT_NEAR(geometry.cotanLaplacian.sum(), 0., 1e-6);
  }
}

TEST_F(HalfedgeGeometrySuite, VertexLumpedMassMatrix) {
  for (auto& asset : {getAsset("lego.ply", false), getAsset("lego.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    IntrinsicGeometryInterface& geometry = *asset.geometry;

    geometry.requireVertexLumpedMassMatrix();

    EXPECT_EQ(geometry.vertexLumpedMassMatrix.rows(), (long int)mesh.nVertices());
    EXPECT_EQ(geometry.vertexLumpedMassMatrix.cols(), (long int)mesh.nVertices());
    EXPECT_EQ(geometry.vertexLumpedMassMatrix.nonZeros(), (long int)(mesh.nVertices()));
  }
}

TEST_F(HalfedgeGeometrySuite, VertexGalerkinMassMatrix) {
  for (auto& asset : {getAsset("lego.ply", false), getAsset("lego.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    IntrinsicGeometryInterface& geometry = *asset.geometry;

    geometry.requireVertexGalerkinMassMatrix();

    EXPECT_EQ(geometry.vertexGalerkinMassMatrix.rows(), (long int)mesh.nVertices());
    EXPECT_EQ(geometry.vertexGalerkinMassMatrix.cols(), (long int)mesh.nVertices());
    EXPECT_EQ(geometry.vertexGalerkinMassMatrix.nonZeros(), (long int)(mesh.nVertices() + 2 * mesh.nEdges()));
  }
}

TEST_F(HalfedgeGeometrySuite, VertexConnectionLaplacian) {
  // for (auto& asset : {getAsset("lego.ply", false), getAsset("lego.ply", true)}) { // TODO nonmanifold not supported
  for (auto& asset : {getAsset("lego.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    IntrinsicGeometryInterface& geometry = *asset.geometry;

    geometry.requireVertexConnectionLaplacian();

    EXPECT_EQ(geometry.vertexConnectionLaplacian.rows(), (long int)mesh.nVertices());
    EXPECT_EQ(geometry.vertexConnectionLaplacian.cols(), (long int)mesh.nVertices());
    EXPECT_EQ(geometry.vertexConnectionLaplacian.nonZeros(), (long int)(mesh.nVertices() + 2. * mesh.nEdges()));
  }
}

TEST_F(HalfedgeGeometrySuite, DECOperators) {
  for (auto& asset : {getAsset("lego.ply", false), getAsset("lego.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    IntrinsicGeometryInterface& geometry = *asset.geometry;

    geometry.requireDECOperators();

    // Eigen::SparseMatrix<double> hodge0, hodge0Inverse, hodge1, hodge1Inverse, hodge2, hodge2Inverse, d0, d1;

    // == Check dimensions
    auto dimensionCheck = [](Eigen::SparseMatrix<double>& m, size_t nRow, size_t nCol) {
      EXPECT_EQ(m.rows(), (long int)nRow);
      EXPECT_EQ(m.cols(), (long int)nCol);
    };

    dimensionCheck(geometry.hodge0, mesh.nVertices(), mesh.nVertices());
    dimensionCheck(geometry.hodge0Inverse, mesh.nVertices(), mesh.nVertices());
    dimensionCheck(geometry.hodge1, mesh.nEdges(), mesh.nEdges());
    dimensionCheck(geometry.hodge1Inverse, mesh.nEdges(), mesh.nEdges());
    dimensionCheck(geometry.hodge2, mesh.nFaces(), mesh.nFaces());
    dimensionCheck(geometry.hodge2Inverse, mesh.nFaces(), mesh.nFaces());
    dimensionCheck(geometry.d0, mesh.nEdges(), mesh.nVertices());
    dimensionCheck(geometry.d1, mesh.nFaces(), mesh.nEdges());
  }
}

TEST_F(HalfedgeGeometrySuite, EdgeDihedralAngles) {
  for (auto& asset : {getAsset("lego.ply", false), getAsset("lego.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    ExtrinsicGeometryInterface& geometry = *asset.geometry;

    geometry.requireEdgeDihedralAngles();
    for (Edge e : mesh.edges()) {
      EXPECT_TRUE(std::isfinite(geometry.edgeDihedralAngles[e]));
    }
  }
}

TEST_F(HalfedgeGeometrySuite, VertexPrincipalCurvatureDirections) {
  for (auto& asset : {getAsset("lego.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    ExtrinsicGeometryInterface& geometry = *asset.geometry;

    geometry.requireVertexPrincipalCurvatureDirections();
    for (Vertex v : mesh.vertices()) {
      EXPECT_TRUE(isfinite(geometry.vertexPrincipalCurvatureDirections[v]));
    }
  }
}

TEST_F(HalfedgeGeometrySuite, FaceNormal) {
  for (auto& asset : {getAsset("lego.ply", false), getAsset("lego.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    EmbeddedGeometryInterface& geometry = *asset.geometry;

    geometry.requireFaceNormals();
    for (Face f : mesh.faces()) {
      EXPECT_TRUE(isfinite(geometry.faceNormals[f]));
      EXPECT_NEAR(norm(geometry.faceNormals[f]), 1., 1e-6);
    }
  }
}

// Test immediate face normal computation from position
TEST_F(HalfedgeGeometrySuite, FaceNormalImmediate_Position) {
  for (auto& asset : {getAsset("lego.ply", false), getAsset("lego.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    IntrinsicGeometryInterface& origGeometry = *asset.geometry;

    // Construct position geometry
    VertexPositionGeometry& geometry = *asset.geometry;

    FaceData<Vector3> faceNormalImmediate(mesh);
    for (Face f : mesh.faces()) {
      faceNormalImmediate[f] = geometry.faceNormal(f);
    }

    geometry.requireFaceNormals();
    for (Face f : mesh.faces()) {
      EXPECT_LT(norm(geometry.faceNormals[f] - faceNormalImmediate[f]), 1e-6);
    }
  }
}

TEST_F(HalfedgeGeometrySuite, VertexNormal) {
  for (auto& asset : {getAsset("lego.ply", false), getAsset("lego.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    EmbeddedGeometryInterface& geometry = *asset.geometry;

    geometry.requireVertexNormals();
    for (Vertex v : mesh.vertices()) {
      EXPECT_TRUE(isfinite(geometry.vertexNormals[v]));
      EXPECT_NEAR(norm(geometry.vertexNormals[v]), 1., 1e-6);
    }
  }
}

TEST_F(HalfedgeGeometrySuite, FaceTangentBasis) {
  for (auto& asset : {getAsset("lego.ply", false), getAsset("lego.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    EmbeddedGeometryInterface& geometry = *asset.geometry;

    geometry.requireFaceTangentBasis();
    geometry.requireFaceNormals();
    for (Face f : mesh.faces()) {
      EXPECT_NEAR(norm(geometry.faceTangentBasis[f][0]), 1., 1e-6);
      EXPECT_NEAR(norm(geometry.faceTangentBasis[f][1]), 1., 1e-6);
      EXPECT_NEAR(dot(geometry.faceTangentBasis[f][0], geometry.faceTangentBasis[f][1]), 0., 1e-6);
      EXPECT_NEAR(dot(geometry.faceTangentBasis[f][0], geometry.faceNormals[f]), 0., 1e-6);
    }
  }
}

TEST_F(HalfedgeGeometrySuite, VertexTangentBasis) {
  for (auto& asset : {getAsset("lego.ply", false), getAsset("lego.ply", true)}) {
    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    EmbeddedGeometryInterface& geometry = *asset.geometry;

    geometry.requireVertexTangentBasis();
    geometry.requireVertexNormals();
    for (Vertex v : mesh.vertices()) {
      EXPECT_NEAR(norm(geometry.vertexTangentBasis[v][0]), 1., 1e-6);
      EXPECT_NEAR(norm(geometry.vertexTangentBasis[v][1]), 1., 1e-6);
      EXPECT_NEAR(dot(geometry.vertexTangentBasis[v][0], geometry.vertexTangentBasis[v][1]), 0., 1e-6);
      EXPECT_NEAR(dot(geometry.vertexTangentBasis[v][0], geometry.vertexNormals[v]), 0., 1e-6);
    }
  }
}

// ============================================================
// =============== Geometry tests
// ============================================================

// More interesting geometry tests which check invariants, etc


// Check that the vertex curvatures return the value expected by Gauss-Bonnet
TEST_F(HalfedgeGeometrySuite, VertexGaussianCurvaturesSum) {
  for (auto& asset : closedMeshes()) {
    if (!asset.isTriangular || !asset.isSubclassManifoldSurfaceMesh) continue;

    asset.printThyName();
    ManifoldSurfaceMesh& mesh = *asset.manifoldMesh;
    IntrinsicGeometryInterface& geometry = *asset.geometry;

    geometry.requireVertexGaussianCurvatures();
    double curvatureSum = 0.;
    for (Vertex e : mesh.vertices()) {
      double s = geometry.vertexGaussianCurvatures[e];
      curvatureSum += s;
    }

    double gaussBonnetCurvature = 2. * PI * mesh.eulerCharacteristic();

    EXPECT_NEAR(curvatureSum, gaussBonnetCurvature, 1e-6);
  }
}

// Check that the face curvatures return the value expected by Gauss-Bonnet
TEST_F(HalfedgeGeometrySuite, FaceGaussianCurvaturesSum) {
  for (auto& asset : closedMeshes()) {
    if (!asset.isTriangular || !asset.isSubclassManifoldSurfaceMesh) continue;

    asset.printThyName();
    ManifoldSurfaceMesh& mesh = *asset.manifoldMesh;

    IntrinsicGeometryInterface& geometry = *asset.geometry;

    geometry.requireFaceGaussianCurvatures();
    double curvatureSum = 0.;
    for (Face f : mesh.faces()) {
      curvatureSum += geometry.faceGaussianCurvatures[f];
    }

    double gaussBonnetCurvature = 2. * PI * mesh.eulerCharacteristic();

    EXPECT_NEAR(curvatureSum, gaussBonnetCurvature, 1e-2);
  }
}

// Test that a bunch of quantities which should sum to the surface area, do
TEST_F(HalfedgeGeometrySuite, SurfaceAreaEquivalence) {
  for (auto& asset : triangularMeshes(false, true)) {

    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    IntrinsicGeometryInterface& geometry = *asset.geometry;

    double tol = 1e-6;

    // Face area
    double surfaceArea_faces = 0.;
    geometry.requireFaceAreas();
    for (Face f : mesh.faces()) {
      surfaceArea_faces += geometry.faceAreas[f];
    }

    // Vertex dual area
    double surfaceArea_vertices = 0.;
    geometry.requireVertexDualAreas();
    for (Vertex v : mesh.vertices()) {
      surfaceArea_vertices += geometry.vertexDualAreas[v];
    }
    EXPECT_NEAR(surfaceArea_faces, surfaceArea_vertices, tol);

    // Lumped mass matrix
    geometry.requireVertexLumpedMassMatrix();
    double surfaceArea_vertexLumpedMass = geometry.vertexLumpedMassMatrix.sum();
    EXPECT_NEAR(surfaceArea_vertices, surfaceArea_vertexLumpedMass, tol);

    // Galerkin mass matrix
    geometry.requireVertexGalerkinMassMatrix();
    double surfaceArea_vertexGalerkinMass = geometry.vertexGalerkinMassMatrix.sum();
    EXPECT_NEAR(surfaceArea_vertexLumpedMass, surfaceArea_vertexGalerkinMass, tol);

    // hodge0
    geometry.requireDECOperators();
    double surfaceArea_hodge0 = geometry.hodge0.sum();
    EXPECT_NEAR(surfaceArea_vertexGalerkinMass, surfaceArea_hodge0, tol);

    // hodge2
    geometry.requireDECOperators();
    double surfaceArea_hodge2 = geometry.hodge2Inverse.sum();
    EXPECT_NEAR(surfaceArea_hodge0, surfaceArea_hodge2, tol);
  }
}


// Build the Laplacian two different ways and ensure that they match up
TEST_F(HalfedgeGeometrySuite, CotanLaplacianEquivalence) {
  for (auto& asset : allMeshes()) {
    if (!asset.isTriangular) continue;

    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    IntrinsicGeometryInterface& geometry = *asset.geometry;

    geometry.requireDECOperators();
    geometry.requireCotanLaplacian();

    Eigen::SparseMatrix<double> weak1FormLaplace = geometry.d0.transpose() * geometry.hodge1 * geometry.d0;

    double difference = (geometry.cotanLaplacian - weak1FormLaplace).cwiseAbs().sum();

    EXPECT_LT(difference, 1e-6);
  }
}

// Make sure that principal curature direction is near-zero on flat and spherical meshes
TEST_F(HalfedgeGeometrySuite, VertexPrincipalCurvatureDirectionsUmbilic) {

  { // flat mesh (with boundary)
    for (auto& asset : {getAsset("lego.ply", true)}) {
      SurfaceMesh& mesh = *asset.mesh;
      ExtrinsicGeometryInterface& geometry = *asset.geometry;

      geometry.requireVertexPrincipalCurvatureDirections();
      for (Vertex v : mesh.vertices()) {
        EXPECT_LT(norm(geometry.vertexPrincipalCurvatureDirections[v]), 1e-5);
      }
    }
  }

  { // sphere mesh
    for (auto& asset : {getAsset("sphere_small.ply", true)}) {
      SurfaceMesh& mesh = *asset.mesh;
      ExtrinsicGeometryInterface& geometry = *asset.geometry;

      geometry.requireVertexPrincipalCurvatureDirections();
      for (Vertex v : mesh.vertices()) {
        EXPECT_LT(norm(geometry.vertexPrincipalCurvatureDirections[v]), 1e-2);
      }
    }
  }
}

// Verify that the normal and tangent directions for each face form an orthonormal frame
TEST_F(HalfedgeGeometrySuite, FaceTangentOrthonormal) {
  for (auto& asset : allMeshes()) {
    if (!asset.isTriangular) continue;

    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    EmbeddedGeometryInterface& geometry = *asset.geometry;

    geometry.requireFaceNormals();
    geometry.requireFaceTangentBasis();

    double tol = 1e-6;

    for (Face f : mesh.faces()) {

      Vector3 normal = geometry.faceNormals[f];
      Vector3 basisX = geometry.faceTangentBasis[f][0];
      Vector3 basisY = geometry.faceTangentBasis[f][1];

      // Check unit
      EXPECT_NEAR(norm(normal), 1., tol);
      EXPECT_NEAR(norm(basisX), 1., tol);
      EXPECT_NEAR(norm(basisY), 1., tol);

      // Check orthogonal
      EXPECT_NEAR(dot(normal, basisX), 0., tol);
      EXPECT_NEAR(dot(normal, basisY), 0., tol);
      EXPECT_NEAR(dot(basisX, basisY), 0., tol);

      // Check handed-ness
      Vector3 crossV = cross(basisX, basisY);
      EXPECT_NEAR(norm(normal - crossV), 0., tol);
    }
  }
}

// Verify that the normal and tangent directions for each vertex form an orthonormal frame
TEST_F(HalfedgeGeometrySuite, VertexTangentOrthonormal) {
  for (auto& asset : allMeshes()) {
    if (!asset.isTriangular) continue;

    asset.printThyName();
    SurfaceMesh& mesh = *asset.mesh;
    EmbeddedGeometryInterface& geometry = *asset.geometry;

    geometry.requireVertexNormals();
    geometry.requireVertexTangentBasis();

    double tol = 1e-6;

    for (Vertex v : mesh.vertices()) {

      Vector3 normal = geometry.vertexNormals[v];
      Vector3 basisX = geometry.vertexTangentBasis[v][0];
      Vector3 basisY = geometry.vertexTangentBasis[v][1];

      // Check unit
      EXPECT_NEAR(norm(normal), 1., tol);
      EXPECT_NEAR(norm(basisX), 1., tol);
      EXPECT_NEAR(norm(basisY), 1., tol);

      // Check orthogonal
      EXPECT_NEAR(dot(normal, basisX), 0., tol);
      EXPECT_NEAR(dot(normal, basisY), 0., tol);
      EXPECT_NEAR(dot(basisX, basisY), 0., tol);

      // Check handed-ness
      Vector3 crossV = cross(basisX, basisY);
      EXPECT_NEAR(norm(normal - crossV), 0., tol);
    }
  }
}

// Ensure that a convex shape has all-positive diheral angles
TEST_F(HalfedgeGeometrySuite, ConvexDiheralAngles) {
  for (auto& asset : {getAsset("tet.obj", false), getAsset("tet.obj", true)}) {
    SurfaceMesh& mesh = *asset.mesh;
    ExtrinsicGeometryInterface& geometry = *asset.geometry;

    geometry.requireEdgeDihedralAngles();
    for (Edge e : mesh.edges()) {
      EXPECT_GT(geometry.edgeDihedralAngles[e], 0.);
    }
  }
}


// ============================================================
// =============== Surface point tests
// ============================================================

TEST_F(HalfedgeGeometrySuite, SurfacePointVertexInSomeFaceTest) {
  for (auto& asset : {getAsset("bob_small.ply", false), getAsset("bob_small.ply", true)}) {
    SurfaceMesh& mesh = *asset.mesh;
    VertexPositionGeometry& geom = *asset.geometry;

    double EPS = 1e-4;

    // Test vertex points
    for (Vertex v : mesh.vertices()) {
      SurfacePoint p(v);
      Vector3 posOrig = p.interpolate(geom.vertexPositions);
      p.validate();

      SurfacePoint pEquiv = p.inSomeFace();
      pEquiv.validate();
      Vector3 posEquiv = pEquiv.interpolate(geom.vertexPositions);

      double dist = (posOrig - posEquiv).norm();

      EXPECT_LT(dist, EPS);
    }
  }
}

TEST_F(HalfedgeGeometrySuite, SurfacePointEdgeInSomeFaceTest) {

  for (auto& asset : {getAsset("bob_small.ply", false), getAsset("bob_small.ply", true)}) {
    SurfaceMesh& mesh = *asset.mesh;
    VertexPositionGeometry& geom = *asset.geometry;

    double EPS = 1e-4;
    std::mt19937 mt(42);

    // Test edge points
    auto unitRand = [&]() {
      std::uniform_real_distribution<double> dist(0.0, 1.0);
      return dist(mt);
    };

    for (Edge e : mesh.edges()) {
      SurfacePoint p(e, unitRand());
      Vector3 posOrig = p.interpolate(geom.vertexPositions);
      p.validate();

      SurfacePoint pEquiv = p.inSomeFace();
      Vector3 posEquiv = pEquiv.interpolate(geom.vertexPositions);

      double dist = (posOrig - posEquiv).norm();
      EXPECT_LT(dist, EPS);
    }
  }
}

TEST_F(HalfedgeGeometrySuite, SurfacePointFaceInSomeFaceTest) {

  for (auto& asset : {getAsset("bob_small.ply", false), getAsset("bob_small.ply", true)}) {
    SurfaceMesh& mesh = *asset.mesh;
    VertexPositionGeometry& geom = *asset.geometry;

    double EPS = 1e-4;
    std::mt19937 mt(42);

    // Test face points
    auto unitRand = [&]() {
      std::uniform_real_distribution<double> dist(0.0, 1.0);
      return dist(mt);
    };
    auto unitBary = [&]() {
      std::uniform_real_distribution<double> dist(0.0, 1.0);
      Vector3 res{unitRand(), unitRand(), unitRand()};
      res /= (res.x + res.y + res.z);
      return res;
    };

    for (Face f : mesh.faces()) {
      SurfacePoint p(f, unitBary());
      Vector3 posOrig = p.interpolate(geom.vertexPositions);
      p.validate();

      SurfacePoint pEquiv = p.inSomeFace();
      Vector3 posEquiv = pEquiv.interpolate(geom.vertexPositions);

      double dist = (posOrig - posEquiv).norm();
      EXPECT_LT(dist, EPS);
    }
  }
}
