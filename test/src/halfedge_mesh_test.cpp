
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/rich_surface_mesh_data.h"

#include "load_test_meshes.h"

#include "gtest/gtest.h"

#include <iostream>
#include <string>
#include <unordered_set>


using namespace geometrycentral;
using namespace geometrycentral::surface;
using std::cout;
using std::endl;


// helpers
namespace {
std::mt19937 mt(42);

template <typename T>
void fillRandom(T& vals) {
  auto unitRand = [&]() {
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(mt);
  };
  for (size_t i = 0; i < vals.size(); i++) {
    vals[i] = unitRand();
  }
}
} // namespace

class HalfedgeMeshSuite : public MeshAssetSuite {};

// ============================================================
// =============== Basic validation tests
// ============================================================

TEST_F(HalfedgeMeshSuite, ValidateClosedMeshTest) {
  for (MeshAsset& a : closedMeshes()) {
    a.printThyName();
    a.mesh->validateConnectivity();
  }
}

TEST_F(HalfedgeMeshSuite, ValidateBoundaryMeshTest) {
  for (MeshAsset& a : boundaryMeshes()) {
    a.printThyName();
    a.mesh->validateConnectivity();
  }
}

// ============================================================
// =============== Constructor tests
// ============================================================

TEST_F(HalfedgeMeshSuite, MatrixConstructorTest) {

  std::unique_ptr<SurfaceMesh> mesh = getAsset("spot.ply", false).mesh;

  { // size_t
    DenseMatrix<size_t> F = mesh->getFaceVertexMatrix<size_t>();
    SurfaceMesh newM(F);
    newM.validateConnectivity();
    EXPECT_EQ(newM.nVertices(), mesh->nVertices());
  }

  { // unsigned int
    DenseMatrix<unsigned int> F = mesh->getFaceVertexMatrix<unsigned int>();
    SurfaceMesh newM(F);
    newM.validateConnectivity();
    EXPECT_EQ(newM.nVertices(), mesh->nVertices());
  }

  { // int
    DenseMatrix<int> F = mesh->getFaceVertexMatrix<int>();
    SurfaceMesh newM(F);
    newM.validateConnectivity();
    EXPECT_EQ(newM.nVertices(), mesh->nVertices());
  }


  { // fixed-size int
    DenseMatrix<int> Fbig = mesh->getFaceVertexMatrix<int>();
    Eigen::Matrix<int, Eigen::Dynamic, 3> F = Fbig;
    SurfaceMesh newM(F);
    newM.validateConnectivity();
    EXPECT_EQ(newM.nVertices(), mesh->nVertices());
  }
}

TEST_F(HalfedgeMeshSuite, MatrixConstructorManifoldTest) {

  std::unique_ptr<SurfaceMesh> mesh = getAsset("spot.ply", false).mesh;

  { // size_t
    DenseMatrix<size_t> F = mesh->getFaceVertexMatrix<size_t>();
    ManifoldSurfaceMesh newM(F);
    newM.validateConnectivity();
    EXPECT_EQ(newM.nVertices(), mesh->nVertices());
  }

  { // unsigned int
    DenseMatrix<unsigned int> F = mesh->getFaceVertexMatrix<unsigned int>();
    ManifoldSurfaceMesh newM(F);
    newM.validateConnectivity();
    EXPECT_EQ(newM.nVertices(), mesh->nVertices());
  }

  { // int
    DenseMatrix<int> F = mesh->getFaceVertexMatrix<int>();
    ManifoldSurfaceMesh newM(F);
    newM.validateConnectivity();
    EXPECT_EQ(newM.nVertices(), mesh->nVertices());
  }


  { // fixed-size int
    DenseMatrix<int> Fbig = mesh->getFaceVertexMatrix<int>();
    Eigen::Matrix<int, Eigen::Dynamic, 3> F = Fbig;
    ManifoldSurfaceMesh newM(F);
    newM.validateConnectivity();
    EXPECT_EQ(newM.nVertices(), mesh->nVertices());
  }
}

// ============================================================
// =============== Range iterator tests
// ============================================================

TEST_F(HalfedgeMeshSuite, IterateVerticesTest) {
  for (MeshAsset& a : allMeshes()) {
    a.printThyName();

    std::unordered_set<Vertex> seenElements;
    for (Vertex e : a.mesh->vertices()) {
      seenElements.insert(e);
    }

    EXPECT_EQ(seenElements.size(), a.mesh->nVertices());
  }
}

TEST_F(HalfedgeMeshSuite, IterateHalfedgesTest) {
  for (MeshAsset& a : allMeshes()) {
    a.printThyName();

    std::unordered_set<Halfedge> seenElements;
    for (Halfedge e : a.mesh->halfedges()) {
      seenElements.insert(e);
    }

    EXPECT_EQ(seenElements.size(), a.mesh->nHalfedges());
  }
}

TEST_F(HalfedgeMeshSuite, IterateInteriorHalfedgesTest) {
  for (MeshAsset& a : allMeshes()) {
    a.printThyName();

    std::unordered_set<Halfedge> seenElements;
    for (Halfedge e : a.mesh->interiorHalfedges()) {
      seenElements.insert(e);
    }

    EXPECT_EQ(seenElements.size(), a.mesh->nInteriorHalfedges());
  }
}

TEST_F(HalfedgeMeshSuite, IterateExteriorHalfedgesTest) {
  for (MeshAsset& a : allMeshes()) {
    a.printThyName();

    std::unordered_set<Halfedge> seenElements;
    for (Halfedge e : a.mesh->exteriorHalfedges()) {
      seenElements.insert(e);
    }

    EXPECT_EQ(seenElements.size(), a.mesh->nExteriorHalfedges());
  }
}

TEST_F(HalfedgeMeshSuite, IterateCornersTest) {
  for (MeshAsset& a : allMeshes()) {
    a.printThyName();

    std::unordered_set<Corner> seenElements;
    for (Corner e : a.mesh->corners()) {
      seenElements.insert(e);
    }

    EXPECT_EQ(seenElements.size(), a.mesh->nCorners());
  }
}

TEST_F(HalfedgeMeshSuite, IterateEdgesTest) {
  for (MeshAsset& a : allMeshes()) {
    a.printThyName();

    std::unordered_set<Edge> seenElements;
    for (Edge e : a.mesh->edges()) {
      seenElements.insert(e);
    }

    EXPECT_EQ(seenElements.size(), a.mesh->nEdges());
  }
}

TEST_F(HalfedgeMeshSuite, IterateFacesTest) {
  for (MeshAsset& a : allMeshes()) {
    a.printThyName();

    std::unordered_set<Face> seenElements;
    for (Face e : a.mesh->faces()) {
      seenElements.insert(e);
    }

    EXPECT_EQ(seenElements.size(), a.mesh->nFaces());
  }
}

TEST_F(HalfedgeMeshSuite, IterateBoundaryLoopsTest) {
  for (MeshAsset& a : allMeshes()) {
    a.printThyName();

    std::unordered_set<BoundaryLoop> seenElements;
    for (BoundaryLoop e : a.mesh->boundaryLoops()) {
      seenElements.insert(e);
    }

    EXPECT_EQ(seenElements.size(), a.mesh->nBoundaryLoops());
  }
}


// ============================================================
// =============== Utility and status functions
// ============================================================

TEST_F(HalfedgeMeshSuite, HasBoundaryTest) {
  for (MeshAsset& a : closedMeshes()) {
    a.printThyName();
    EXPECT_FALSE(a.mesh->hasBoundary());
  }
  for (MeshAsset& a : boundaryMeshes()) {
    a.printThyName();
    EXPECT_TRUE(a.mesh->hasBoundary());
  }
}


TEST_F(HalfedgeMeshSuite, IsTriangularTest) {
  EXPECT_EQ(getAsset("tet.obj", false).mesh->isTriangular(), true);
  EXPECT_EQ(getAsset("spot.ply", false).mesh->isTriangular(), true);
  EXPECT_EQ(getAsset("dodecahedron_poly.obj", false).mesh->isTriangular(), false);
  EXPECT_EQ(getAsset("platonic_shelf.obj", false).mesh->isTriangular(), false);
  EXPECT_EQ(getAsset("bob_small.ply", false).mesh->isTriangular(), true);
  EXPECT_EQ(getAsset("lego.ply", false).mesh->isTriangular(), true);
}

TEST_F(HalfedgeMeshSuite, EulerCharacteristicTest) {
  EXPECT_EQ(getAsset("tet.obj", true).manifoldMesh->eulerCharacteristic(), 2);
  EXPECT_EQ(getAsset("spot.ply", true).manifoldMesh->eulerCharacteristic(), 2);
  EXPECT_EQ(getAsset("dodecahedron_poly.obj", true).manifoldMesh->eulerCharacteristic(), 2);
  EXPECT_EQ(getAsset("bob_small.ply", true).manifoldMesh->eulerCharacteristic(), 0);
}

TEST_F(HalfedgeMeshSuite, GenusTest) {
  EXPECT_EQ(getAsset("tet.obj", true).manifoldMesh->genus(), 0);
  EXPECT_EQ(getAsset("spot.ply", true).manifoldMesh->genus(), 0);
  EXPECT_EQ(getAsset("dodecahedron_poly.obj", true).manifoldMesh->genus(), 0);
  EXPECT_EQ(getAsset("bob_small.ply", true).manifoldMesh->genus(), 1);
}

TEST_F(HalfedgeMeshSuite, ConnectedComponentsTest) {
  EXPECT_EQ(getAsset("tet.obj", false).mesh->nConnectedComponents(), 1);
  EXPECT_EQ(getAsset("spot.ply", false).mesh->nConnectedComponents(), 1);
  EXPECT_EQ(getAsset("dodecahedron_poly.obj", false).mesh->nConnectedComponents(), 1);
  EXPECT_EQ(getAsset("platonic_shelf.obj", false).mesh->nConnectedComponents(), 5);
  EXPECT_EQ(getAsset("bob_small.ply", false).mesh->nConnectedComponents(), 1);
  EXPECT_EQ(getAsset("lego.ply", false).mesh->nConnectedComponents(), 1);
}


TEST_F(HalfedgeMeshSuite, PrintElementTest) {
  MeshAsset a = getAsset("lego.ply", true);
  SurfaceMesh& mesh = *a.mesh;
  std::cout << mesh.halfedge(0) << std::endl;
  std::cout << mesh.corner(0) << std::endl;
  std::cout << mesh.vertex(0) << std::endl;
  std::cout << mesh.edge(0) << std::endl;
  std::cout << mesh.face(0) << std::endl;
  std::cout << mesh.boundaryLoop(0) << std::endl;
}


TEST_F(HalfedgeMeshSuite, PrintElementStringTest) {
  MeshAsset a = getAsset("lego.ply", true);
  SurfaceMesh& mesh = *a.mesh;
  std::cout << std::to_string(mesh.halfedge(0)) << std::endl;
  std::cout << std::to_string(mesh.corner(0)) << std::endl;
  std::cout << std::to_string(mesh.vertex(0)) << std::endl;
  std::cout << std::to_string(mesh.edge(0)) << std::endl;
  std::cout << std::to_string(mesh.face(0)) << std::endl;
  std::cout << std::to_string(mesh.boundaryLoop(0)) << std::endl;
}


// ============================================================
// =============== Containers
// ============================================================


// Make sure that nothing explodes if we delete the mesh before the container
TEST_F(HalfedgeMeshSuite, ContainerMeshDestructTest) {
  std::unique_ptr<SurfaceMesh> mesh = getAsset("spot.ply", false).mesh;

  {
    VertexData<double> testD(*mesh);
    for (Vertex v : mesh->vertices()) {
      testD[v] = 42.0;
    }

    mesh.reset(); // delete mesh

  } // scope block triggers testD delete

  ASSERT_EQ(2 + 2, 4); // debugging is easier if failure isn't last
}


TEST_F(HalfedgeMeshSuite, ContainerAccessTest) {
  std::unique_ptr<SurfaceMesh> mesh = getAsset("spot.ply", false).mesh;

  {
    VertexData<double> testDA(*mesh);
    VertexData<double> testDB(*mesh);
    for (Vertex v : mesh->vertices()) {
      testDA[v] = 3.0;
      testDB[v] = 2.0;
    }

    testDA.raw() += testDB.raw();

    for (Vertex v : mesh->vertices()) {
      ASSERT_EQ(testDA[v], 5.);
    }

  } // scope block triggers testD delete
}


// == Container arithmetic
//
// Note: these tests cover many, but not all of the arithmetic operations. They are generated by macro, so
// (theoretically) if one is correct they all should be.
// TODO add tests on non-compressed meshes to make sure there's nothing funky going on there.

TEST_F(HalfedgeMeshSuite, ContainerArithmeticPlus) {
  std::unique_ptr<SurfaceMesh> mesh = getAsset("spot.ply", false).mesh;

  {
    { // basic
      VertexData<double> A(*mesh, 1.);
      VertexData<double> B(*mesh, 2.);
      VertexData<double> C = A + B;
      for (Vertex v : mesh->vertices()) ASSERT_EQ(C[v], 3.);
    }

    { // different type
      VertexData<double> A(*mesh, 1.);
      VertexData<float> B(*mesh, 2.);
      VertexData<double> C = A + B;
      for (Vertex v : mesh->vertices()) ASSERT_EQ(C[v], 3.);
    }

    { // +=
      VertexData<double> A(*mesh, 1.);
      VertexData<float> B(*mesh, 2.);
      A += B;
      for (Vertex v : mesh->vertices()) ASSERT_EQ(A[v], 3.);
    }

    { // scalar
      VertexData<double> A(*mesh, 1.);
      double b = 2;
      float c = 3;
      A = (A + b);
      A = (b + A);
      A += b;
      A = (A + c);
      A = (c + A);
      A += c;
      for (Vertex v : mesh->vertices()) ASSERT_EQ(A[v], 16.);
    }
  }
}

TEST_F(HalfedgeMeshSuite, ContainerArithmeticMinus) {
  std::unique_ptr<SurfaceMesh> mesh = getAsset("spot.ply", false).mesh;

  { // operator-

    { // basic
      VertexData<double> A(*mesh, 1.);
      VertexData<double> B(*mesh, 2.);
      VertexData<double> C = A - B;
      for (Vertex v : mesh->vertices()) ASSERT_EQ(C[v], -1.);
    }

    { // different type
      VertexData<double> A(*mesh, 1.);
      VertexData<float> B(*mesh, 2.);
      VertexData<double> C = A - B;
      for (Vertex v : mesh->vertices()) ASSERT_EQ(C[v], -1.);
    }

    { // -=
      VertexData<double> A(*mesh, 1.);
      VertexData<float> B(*mesh, 2.);
      A -= B;
      for (Vertex v : mesh->vertices()) ASSERT_EQ(A[v], -1.);
    }

    { // scalar
      VertexData<double> A(*mesh, 1.);
      double b = 2;
      float c = 3;
      A = (A - b); // -1
      A = (b - A); // 3
      A -= b;      // 1
      A = (A - c); // -2
      A = (c - A); // 5
      A -= c;      // 2
      for (Vertex v : mesh->vertices()) ASSERT_EQ(A[v], 2.);
    }
  }
}

TEST_F(HalfedgeMeshSuite, ContainerArithmeticMult) {
  std::unique_ptr<SurfaceMesh> mesh = getAsset("spot.ply", false).mesh;

  { // operator*

    { // basic
      VertexData<double> A(*mesh, 3.);
      VertexData<double> B(*mesh, 2.);
      VertexData<double> C = A * B;
      for (Vertex v : mesh->vertices()) ASSERT_EQ(C[v], 6.);
    }

    { // different type
      VertexData<double> A(*mesh, 2.);
      VertexData<float> B(*mesh, 3.);
      VertexData<double> C = A * B;
      for (Vertex v : mesh->vertices()) ASSERT_EQ(C[v], 6.);
    }

    { // very different type
      VertexData<float> A(*mesh, 2.);
      VertexData<Vector3> B(*mesh, Vector3{1., 2., 3.});
      VertexData<Vector3> C = A * B;
      Vector3 r{2., 4., 6.};
      for (Vertex v : mesh->vertices()) ASSERT_EQ(C[v], r);
    }

    { // *=
      VertexData<double> A(*mesh, 2.);
      VertexData<float> B(*mesh, 3.);
      A *= B;
      for (Vertex v : mesh->vertices()) ASSERT_EQ(A[v], 6.);
    }

    { // scalar
      VertexData<double> A(*mesh, 1.);
      double b = 2;
      float c = 3;
      A = (A * b); // 2
      A = (b * A); // 4
      A *= b;      // 8
      A = (A * c); // 24
      A = (c * A); // 72
      A *= c;      // 216
      for (Vertex v : mesh->vertices()) ASSERT_EQ(A[v], 216.);
    }
  }
}

TEST_F(HalfedgeMeshSuite, ContainerArithmeticDiv) {
  std::unique_ptr<SurfaceMesh> mesh = getAsset("spot.ply", false).mesh;

  { // operator /

    { // basic
      VertexData<double> A(*mesh, 3.);
      VertexData<double> B(*mesh, 2.);
      VertexData<double> C = A / B;
      for (Vertex v : mesh->vertices()) ASSERT_EQ(C[v], 1.5);
    }

    { // different type
      VertexData<double> A(*mesh, 3.);
      VertexData<float> B(*mesh, 2.);
      VertexData<double> C = A / B;
      for (Vertex v : mesh->vertices()) ASSERT_EQ(C[v], 1.5);
    }

    { // /=
      VertexData<double> A(*mesh, 3.);
      VertexData<float> B(*mesh, 2.);
      A /= B;
      for (Vertex v : mesh->vertices()) ASSERT_EQ(A[v], 1.5);
    }

    { // scalar
      VertexData<double> A(*mesh, 1.);
      double b = 1;
      float c = 1;
      A = (A / b);
      A = (b / A);
      A /= b;
      A = (A / c);
      A = (c / A);
      A /= c;
      for (Vertex v : mesh->vertices()) ASSERT_EQ(A[v], 1.); // something
    }
  }
}

TEST_F(HalfedgeMeshSuite, ContainerArithmeticMod) {
  std::unique_ptr<SurfaceMesh> mesh = getAsset("spot.ply", false).mesh;

  { // operator %

    { // basic
      VertexData<int> A(*mesh, 3);
      VertexData<int> B(*mesh, 2);
      VertexData<int> C = A % B;
      for (Vertex v : mesh->vertices()) ASSERT_EQ(C[v], 1);
    }

    { // different type
      VertexData<int> A(*mesh, 3);
      VertexData<int> B(*mesh, 2);
      VertexData<int> C = A % B;
      for (Vertex v : mesh->vertices()) ASSERT_EQ(C[v], 1);
    }

    { // %=
      VertexData<int> A(*mesh, 3);
      VertexData<int> B(*mesh, 2);
      A %= B;
      for (Vertex v : mesh->vertices()) ASSERT_EQ(A[v], 1);
    }

    { // scalar
      // (have to worry about arithmetic exceptions here)
      VertexData<int> A(*mesh, 9);
      int b = 7;
      A = (A % b);
      A = VertexData<int>(*mesh, 9);
      A = (b % A);
      A = VertexData<int>(*mesh, 9);
      A %= b;
      A = VertexData<int>(*mesh, 9);
      A = (A % b);
      A = VertexData<int>(*mesh, 9);
      A = (b % A);
      A = VertexData<int>(*mesh, 9);
      A %= b;
      A = VertexData<int>(*mesh, 9);
      // for (Vertex v : mesh->vertices()) ASSERT_EQ(A[v], 1.); // = something?
    }
  }
}

TEST_F(HalfedgeMeshSuite, ContainerArithmeticUnary) {
  std::unique_ptr<SurfaceMesh> mesh = getAsset("spot.ply", false).mesh;

  { // operator +
    VertexData<double> A(*mesh, 3.);
    VertexData<double> B = +A;
    for (Vertex v : mesh->vertices()) ASSERT_EQ(B[v], 3);
  }

  { // operator -
    VertexData<double> A(*mesh, 3.);
    VertexData<double> B = -A;
    for (Vertex v : mesh->vertices()) ASSERT_EQ(B[v], -3);
  }

  { // operator !
    VertexData<bool> A(*mesh, false);
    VertexData<bool> B = !A;
    for (Vertex v : mesh->vertices()) ASSERT_EQ(B[v], true);
  }
}


// ============================================================
// =============== Navigators
// ============================================================

TEST_F(HalfedgeMeshSuite, PrevTest) {
  for (MeshAsset& a : allMeshes()) {
    a.printThyName();

    for (Halfedge he : a.mesh->halfedges()) {
      Halfedge next = he.next();
      EXPECT_EQ(next.prevOrbitFace(), he);
      EXPECT_EQ(next.prevOrbitVertex(), he); // doesn't necessarily work on nonmanifold
    }
  }
}

TEST_F(HalfedgeMeshSuite, VertexAdjacentNavigator) {
  for (MeshAsset& a : allMeshes()) {
    a.printThyName();

    size_t degreeTot = 0;
    for (Vertex v : a.mesh->vertices()) {
      for (Halfedge he : v.outgoingHalfedges()) {
        degreeTot++;
      }
    }

    EXPECT_EQ(degreeTot, a.mesh->nHalfedges());
  }
}


TEST_F(HalfedgeMeshSuite, VertexCornerNavigatorInterior) {
  for (MeshAsset& a : allMeshes()) {
    a.printThyName();

    for (Vertex v : a.mesh->vertices()) {
      for (Corner c : v.adjacentCorners()) {
        EXPECT_FALSE(c.face().isBoundaryLoop());
      }
    }
  }
}


TEST_F(HalfedgeMeshSuite, VertexEdgeNavigator) {
  for (MeshAsset& a : polygonalComplexMeshes()) {
    a.printThyName();

    EdgeData<size_t> edgeCount(*a.mesh, 0);

    for (Vertex v : a.mesh->vertices()) {

      // make sure each edge is returned just once
      std::unordered_set<Edge> seen;

      for (Edge e : v.adjacentEdges()) {
        EXPECT_TRUE(seen.find(e) == seen.end());
        seen.insert(e);
        edgeCount[e]++;
      }
    }

    for (Edge e : a.mesh->edges()) {
      EXPECT_EQ(edgeCount[e], 2);
    }
  }
}

TEST_F(HalfedgeMeshSuite, VertexFaceNavigator) {
  for (MeshAsset& a : polygonalComplexMeshes()) {
    a.printThyName();

    FaceData<size_t> faceCount(*a.mesh, 0);

    for (Vertex v : a.mesh->vertices()) {

      // make sure each face is returned just once
      std::unordered_set<Face> seen;

      for (Face f : v.adjacentFaces()) {
        EXPECT_TRUE(seen.find(f) == seen.end());
        seen.insert(f);
        faceCount[f]++;
      }
    }

    for (Face f : a.mesh->faces()) {
      EXPECT_EQ(faceCount[f], f.degree());
    }
  }
}

TEST_F(HalfedgeMeshSuite, FaceFaceNavigator) {
  for (MeshAsset& a : polygonalComplexMeshes()) {
    a.printThyName();

    FaceData<size_t> faceCount(*a.mesh, 0);

    for (Face f : a.mesh->faces()) {

      // make sure each face is returned just once
      std::unordered_set<Face> seen;

      for (Face fn : f.adjacentFaces()) {
        EXPECT_TRUE(seen.find(fn) == seen.end());
        seen.insert(fn);
        faceCount[fn]++;
      }
    }

    for (Face f : a.mesh->faces()) {

      size_t expectedCount = 0;
      for (Edge e : f.adjacentEdges()) {
        expectedCount += e.degree() - 1;
      }

      EXPECT_EQ(faceCount[f], expectedCount);
    }
  }
}

TEST_F(HalfedgeMeshSuite, EdgeVertexNavigator) {
  // test firstVertex()/secondVertex() and e.adjacentVertices()

  std::unique_ptr<SurfaceMesh> mesh = getAsset("spot.ply", false).mesh;

  for (Edge e : mesh->edges()) {
    EXPECT_FALSE(e.firstVertex() == e.secondVertex());

    Vertex firstV;
    for (Vertex v : e.adjacentVertices()) {
      if (firstV == Vertex()) {
        firstV = v;
      } else {
        EXPECT_FALSE(v == firstV);
      }

      EXPECT_TRUE(v == e.firstVertex() || v == e.secondVertex());
    }
  }
}

TEST_F(HalfedgeMeshSuite, EdgeDiamondNavigator) {

  std::unique_ptr<SurfaceMesh> mesh = getAsset("spot.ply", false).mesh;

  for (Edge e : mesh->edges()) {

    // make sure each halfedge only shows up once
    std::unordered_set<Halfedge> badHalfedges;
    badHalfedges.insert(e.halfedge());
    badHalfedges.insert(e.halfedge().twin());
    for (Halfedge he : e.diamondBoundary()) {
      EXPECT_EQ(badHalfedges.find(he), badHalfedges.end());
      badHalfedges.insert(he);

      EXPECT_TRUE(he.face() == e.halfedge().face() || he.face() == e.halfedge().twin().face());
    }
  }
}


// ============================================================
// =============== Utilities
// ============================================================

TEST_F(HalfedgeMeshSuite, IsManifoldOrientedTest) {

  {
    auto asset = getAsset("lego.ply", false);
    SurfaceMesh& mesh = *asset.mesh;
    EXPECT_TRUE(mesh.isEdgeManifold());
    EXPECT_TRUE(mesh.isManifold());
    EXPECT_TRUE(mesh.isOriented());
  }

  {
    auto asset = getAsset("hourglass_ico.obj", false);
    SurfaceMesh& mesh = *asset.mesh;
    EXPECT_TRUE(mesh.isEdgeManifold());
    EXPECT_FALSE(mesh.isManifold());
    EXPECT_TRUE(mesh.isOriented());
  }

  {
    auto asset = getAsset("triple_vierbein.obj", false);
    SurfaceMesh& mesh = *asset.mesh;
    EXPECT_FALSE(mesh.isEdgeManifold());
    EXPECT_FALSE(mesh.isManifold());
    EXPECT_FALSE(mesh.isOriented());
  }

  {
    auto asset = getAsset("moebius.obj", false);
    SurfaceMesh& mesh = *asset.mesh;
    EXPECT_TRUE(mesh.isEdgeManifold());
    EXPECT_TRUE(mesh.isManifold());
    EXPECT_FALSE(mesh.isOriented());
  }
}

// ============================================================
// =============== Rich mesh
// ============================================================

TEST_F(HalfedgeMeshSuite, RichMeshDataSaveLoadProperties) {

  for (auto& asset : {getAsset("lego.ply", false), getAsset("lego.ply", true)}) {

    SurfaceMesh& mesh = *asset.mesh;
    VertexPositionGeometry& geom = *asset.geometry;

    // Initialize some values
    VertexData<double> vertexValues(mesh);
    fillRandom(vertexValues);
    EdgeData<double> edgeValues(mesh);
    fillRandom(edgeValues);
    HalfedgeData<double> halfedgeValues(mesh);
    fillRandom(halfedgeValues);
    CornerData<double> cornerValues(mesh);
    fillRandom(cornerValues);
    FaceData<double> faceValues(mesh);
    fillRandom(faceValues);
    BoundaryLoopData<double> blValues(mesh);
    fillRandom(blValues);

    RichSurfaceMeshData richData(mesh);

    richData.addVertexProperty("name_1", vertexValues);
    richData.addEdgeProperty("name_2", edgeValues);
    richData.addHalfedgeProperty("name_3", halfedgeValues);
    richData.addCornerProperty("name_4", cornerValues);
    richData.addFaceProperty("name_5", faceValues);
    richData.addBoundaryLoopProperty("name_6", blValues);

    // Write the data to file
    richData.write("test_archive.ply");

    // Read the data back from file
    RichSurfaceMeshData richDataIn(mesh, "test_archive.ply");

    // Get the properties
    VertexData<double> vertexValuesIn = richDataIn.getVertexProperty<double>("name_1");
    EdgeData<double> edgeValuesIn = richDataIn.getEdgeProperty<double>("name_2");
    HalfedgeData<double> halfedgeValuesIn = richDataIn.getHalfedgeProperty<double>("name_3");
    CornerData<double> cornerValuesIn = richDataIn.getCornerProperty<double>("name_4");
    FaceData<double> faceValuesIn = richDataIn.getFaceProperty<double>("name_5");
    BoundaryLoopData<double> blValuesIn = richDataIn.getBoundaryLoopProperty<double>("name_6");

    for (Vertex v : mesh.vertices()) EXPECT_EQ(vertexValues[v], vertexValuesIn[v]);
    for (Edge e : mesh.edges()) EXPECT_EQ(edgeValues[e], edgeValuesIn[e]);
    for (Halfedge he : mesh.halfedges()) EXPECT_EQ(halfedgeValues[he], halfedgeValuesIn[he]);
    for (Corner c : mesh.corners()) EXPECT_EQ(cornerValues[c], cornerValuesIn[c]);
    for (Face f : mesh.faces()) EXPECT_EQ(faceValues[f], faceValuesIn[f]);
    for (BoundaryLoop bl : mesh.boundaryLoops()) EXPECT_EQ(blValues[bl], blValuesIn[bl]);
  }
}

// TODO test these after a deletion

TEST_F(HalfedgeMeshSuite, RichMeshDataSaveLoadMeshGeneral) {

  for (auto& asset : {getAsset("lego.ply", false), getAsset("lego.ply", true)}) {

    SurfaceMesh& mesh = *asset.mesh;
    VertexPositionGeometry& geom = *asset.geometry;

    HalfedgeData<double> halfedgeValues(mesh);
    fillRandom(halfedgeValues);

    // Write the data to file
    RichSurfaceMeshData richData(mesh);
    richData.addMeshConnectivity();
    richData.addGeometry(geom);
    richData.addHalfedgeProperty("he_vals", halfedgeValues);
    richData.write("test_archive.ply");

    // Read the data back from file
    std::unique_ptr<SurfaceMesh> meshIn;
    std::unique_ptr<RichSurfaceMeshData> richDataIn;
    std::unique_ptr<VertexPositionGeometry> geomIn;
    std::tie(meshIn, richDataIn) = RichSurfaceMeshData::readMeshAndData("test_archive.ply");

    // Do some basic sanity checks on the mesh
    meshIn->validateConnectivity();
    ASSERT_EQ(mesh.nVertices(), meshIn->nVertices());
    ASSERT_EQ(mesh.nHalfedges(), meshIn->nHalfedges());
    ASSERT_EQ(mesh.nEdges(), meshIn->nEdges());
    ASSERT_EQ(mesh.nFaces(), meshIn->nFaces());
    ASSERT_EQ(mesh.nBoundaryLoops(), meshIn->nBoundaryLoops());

    // Check contained data and properties
    geomIn = richDataIn->getGeometry();
    HalfedgeData<double> halfedgeValuesIn = richDataIn->getHalfedgeProperty<double>("he_vals");
    for (size_t i = 0; i < mesh.nHalfedges(); i++) EXPECT_EQ(halfedgeValues[i], halfedgeValuesIn[i]);
  }
}

TEST_F(HalfedgeMeshSuite, RichMeshDataSaveLoadMeshManifold) {

  for (auto& asset : {getAsset("lego.ply", false), getAsset("lego.ply", true)}) {

    SurfaceMesh& mesh = *asset.mesh;
    VertexPositionGeometry& geom = *asset.geometry;

    HalfedgeData<double> halfedgeValues(mesh);
    fillRandom(halfedgeValues);

    // Write the data to file
    RichSurfaceMeshData richData(mesh);
    richData.addMeshConnectivity();
    richData.addGeometry(geom);
    richData.addHalfedgeProperty("he_vals", halfedgeValues);
    richData.write("test_archive.ply");

    // Read the data back from file
    std::unique_ptr<SurfaceMesh> meshIn;
    std::unique_ptr<RichSurfaceMeshData> richDataIn;
    std::unique_ptr<VertexPositionGeometry> geomIn;
    std::tie(meshIn, richDataIn) = RichSurfaceMeshData::readMeshAndData("test_archive.ply");

    // Do some basic sanity checks on the mesh
    meshIn->validateConnectivity();
    ASSERT_EQ(mesh.nVertices(), meshIn->nVertices());
    ASSERT_EQ(mesh.nHalfedges(), meshIn->nHalfedges());
    ASSERT_EQ(mesh.nEdges(), meshIn->nEdges());
    ASSERT_EQ(mesh.nFaces(), meshIn->nFaces());
    ASSERT_EQ(mesh.nBoundaryLoops(), meshIn->nBoundaryLoops());

    // Check contained data and properties
    geomIn = richDataIn->getGeometry();
    HalfedgeData<double> halfedgeValuesIn = richDataIn->getHalfedgeProperty<double>("he_vals");
    for (size_t i = 0; i < mesh.nHalfedges(); i++) EXPECT_EQ(halfedgeValues[i], halfedgeValuesIn[i]);
  }
}
