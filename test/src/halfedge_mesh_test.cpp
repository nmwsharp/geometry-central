
#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/meshio.h"

#include "load_test_meshes.h"

#include "gtest/gtest.h"

#include <iostream>
#include <string>
#include <unordered_set>


using namespace geometrycentral;
using namespace geometrycentral::surface;
using std::cout;
using std::endl;

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
  EXPECT_EQ(getAsset("tet.obj").mesh->isTriangular(), true);
  EXPECT_EQ(getAsset("spot.ply").mesh->isTriangular(), true);
  EXPECT_EQ(getAsset("dodecahedron_poly.obj").mesh->isTriangular(), false);
  EXPECT_EQ(getAsset("platonic_shelf.obj").mesh->isTriangular(), false);
  EXPECT_EQ(getAsset("bob_small.ply").mesh->isTriangular(), true);
  EXPECT_EQ(getAsset("lego.ply").mesh->isTriangular(), true);
}

TEST_F(HalfedgeMeshSuite, EulerCharacteristicTest) {
  EXPECT_EQ(getAsset("tet.obj").mesh->eulerCharacteristic(), 2);
  EXPECT_EQ(getAsset("spot.ply").mesh->eulerCharacteristic(), 2);
  EXPECT_EQ(getAsset("dodecahedron_poly.obj").mesh->eulerCharacteristic(), 2);
  EXPECT_EQ(getAsset("bob_small.ply").mesh->eulerCharacteristic(), 0);
}

TEST_F(HalfedgeMeshSuite, GenusTest) {
  EXPECT_EQ(getAsset("tet.obj").mesh->genus(), 0);
  EXPECT_EQ(getAsset("spot.ply").mesh->genus(), 0);
  EXPECT_EQ(getAsset("dodecahedron_poly.obj").mesh->genus(), 0);
  EXPECT_EQ(getAsset("bob_small.ply").mesh->genus(), 1);
}

TEST_F(HalfedgeMeshSuite, ConnectedComponentsTest) {
  EXPECT_EQ(getAsset("tet.obj").mesh->nConnectedComponents(), 1);
  EXPECT_EQ(getAsset("spot.ply").mesh->nConnectedComponents(), 1);
  EXPECT_EQ(getAsset("dodecahedron_poly.obj").mesh->nConnectedComponents(), 1);
  EXPECT_EQ(getAsset("platonic_shelf.obj").mesh->nConnectedComponents(), 5);
  EXPECT_EQ(getAsset("bob_small.ply").mesh->nConnectedComponents(), 1);
  EXPECT_EQ(getAsset("lego.ply").mesh->nConnectedComponents(), 1);
}


// ============================================================
// =============== Containers
// ============================================================


// Make sure that nothing explodes if we delete the mesh before the container
TEST_F(HalfedgeMeshSuite, ContainerMeshDestructTest) {
  std::unique_ptr<HalfedgeMesh> mesh = getAsset("spot.ply").mesh;

  {
    VertexData<double> testD(*mesh);
    for (Vertex v : mesh->vertices()) {
      testD[v] = 42.0;
    }

    mesh.reset(); // delete mesh

  } // scope block triggers testD delete

  ASSERT_EQ(2 + 2, 4); // debugging is easier if failure isn't last
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
