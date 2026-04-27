#include "geometrycentral/surface/flip_geodesics.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "load_test_meshes.h"

#include "gtest/gtest.h"

#include <string>

using namespace geometrycentral;
using namespace geometrycentral::surface;

class FlipGeodesicSuite : public MeshAssetSuite {};


// ============================================================
// =============== Construction
// ============================================================

TEST_F(FlipGeodesicSuite, ConstructFromDijkstra) {
  for (const MeshAsset& a : {getAsset("fox.ply", true), getAsset("cat_head.obj", true)}) {
    a.printThyName();
    ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
    VertexPositionGeometry& geom = *a.geometry;

    std::unique_ptr<FlipEdgeNetwork> network =
        FlipEdgeNetwork::constructFromDijkstraPath(mesh, geom, mesh.vertex(0), mesh.vertex(mesh.nVertices() / 2));

    ASSERT_NE(network, nullptr);
    ASSERT_EQ(network->paths.size(), 1);
    EXPECT_GT(network->paths[0]->size(), 0);
  }
}

TEST_F(FlipGeodesicSuite, ConstructFromPiecewiseDijkstra) {
  for (const MeshAsset& a : {getAsset("fox.ply", true), getAsset("cat_head.obj", true)}) {
    a.printThyName();
    ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
    VertexPositionGeometry& geom = *a.geometry;

    size_t nV = mesh.nVertices();
    std::vector<Vertex> waypoints = {mesh.vertex(0), mesh.vertex(nV / 3), mesh.vertex(2 * nV / 3)};

    std::unique_ptr<FlipEdgeNetwork> network =
        FlipEdgeNetwork::constructFromPiecewiseDijkstraPath(mesh, geom, waypoints, false);

    ASSERT_NE(network, nullptr);
    ASSERT_EQ(network->paths.size(), 1);
    EXPECT_GT(network->paths[0]->size(), 0);
  }
}

TEST_F(FlipGeodesicSuite, ConstructClosedLoop) {
  for (const MeshAsset& a : {getAsset("fox.ply", true), getAsset("cat_head.obj", true)}) {
    a.printThyName();
    ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
    VertexPositionGeometry& geom = *a.geometry;

    size_t nV = mesh.nVertices();
    std::vector<Vertex> waypoints = {mesh.vertex(0), mesh.vertex(nV / 3), mesh.vertex(2 * nV / 3)};

    std::unique_ptr<FlipEdgeNetwork> network =
        FlipEdgeNetwork::constructFromPiecewiseDijkstraPath(mesh, geom, waypoints, /*closed=*/true);

    ASSERT_NE(network, nullptr);
    ASSERT_EQ(network->paths.size(), 1);
    EXPECT_TRUE(network->paths[0]->isClosed);
  }
}


// ============================================================
// =============== Shortening
// ============================================================

TEST_F(FlipGeodesicSuite, ShortenedPathIsStraight) {
  for (const MeshAsset& a : {getAsset("fox.ply", true), getAsset("cat_head.obj", true)}) {
    a.printThyName();
    ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
    VertexPositionGeometry& geom = *a.geometry;

    std::unique_ptr<FlipEdgeNetwork> network =
        FlipEdgeNetwork::constructFromDijkstraPath(mesh, geom, mesh.vertex(0), mesh.vertex(mesh.nVertices() / 2));
    ASSERT_NE(network, nullptr);

    network->iterativeShorten();

    EXPECT_GE(network->minAngle(), M_PI - network->EPS_ANGLE);
  }
}

TEST_F(FlipGeodesicSuite, ShorteningDoesNotIncreaseLength) {
  for (const MeshAsset& a : {getAsset("fox.ply", true), getAsset("cat_head.obj", true)}) {
    a.printThyName();
    ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
    VertexPositionGeometry& geom = *a.geometry;

    std::unique_ptr<FlipEdgeNetwork> network =
        FlipEdgeNetwork::constructFromDijkstraPath(mesh, geom, mesh.vertex(0), mesh.vertex(mesh.nVertices() / 2));
    ASSERT_NE(network, nullptr);

    double lengthBefore = network->length();
    network->iterativeShorten();
    double lengthAfter = network->length();

    EXPECT_LE(lengthAfter, lengthBefore + 1e-6);
  }
}

// ============================================================
// =============== Piecewise Dijkstra overlap removal
// ============================================================

// See: https://github.com/nmwsharp/geometry-central/pull/237

// Small 5-triangle mesh for overlap-removal tests
namespace {
std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>> buildOverlapTestMesh() {
  std::vector<Vector3> verts = {{0, 0, 0}, {0, 3, 0}, {0, 3, 1}, {0, 5, 0}, {0, 5, 5}, {0, 3, 2}, {0, 5, 0}};
  std::vector<std::vector<size_t>> faces = {{0, 1, 6}, {1, 3, 4}, {2, 1, 4}, {2, 4, 5}, {1, 2, 6}};
  return makeManifoldSurfaceMeshAndGeometry(faces, verts);
}
} // namespace

TEST_F(FlipGeodesicSuite, PiecewiseDijkstraSingleDetour) {
  std::unique_ptr<ManifoldSurfaceMesh> meshPtr;
  std::unique_ptr<VertexPositionGeometry> geomPtr;
  std::tie(meshPtr, geomPtr) = buildOverlapTestMesh();
  ManifoldSurfaceMesh& mesh = *meshPtr;
  VertexPositionGeometry& geom = *geomPtr;

  std::vector<Vertex> path = {mesh.vertex(0), mesh.vertex(2), mesh.vertex(3)};
  std::unique_ptr<FlipEdgeNetwork> network =
      FlipEdgeNetwork::constructFromPiecewiseDijkstraPath(mesh, geom, path, false);

  ASSERT_EQ(network->paths.size(), 1);
  EXPECT_EQ(network->paths[0]->getHalfedgeList().size(), 2);
}

TEST_F(FlipGeodesicSuite, PiecewiseDijkstraMultipleOverlapAtStart) {
  std::unique_ptr<ManifoldSurfaceMesh> meshPtr;
  std::unique_ptr<VertexPositionGeometry> geomPtr;
  std::tie(meshPtr, geomPtr) = buildOverlapTestMesh();
  ManifoldSurfaceMesh& mesh = *meshPtr;
  VertexPositionGeometry& geom = *geomPtr;

  std::vector<Vertex> path = {mesh.vertex(1), mesh.vertex(5), mesh.vertex(3)};
  std::unique_ptr<FlipEdgeNetwork> network =
      FlipEdgeNetwork::constructFromPiecewiseDijkstraPath(mesh, geom, path, false);

  ASSERT_EQ(network->paths.size(), 1);
  EXPECT_EQ(network->paths[0]->getHalfedgeList().size(), 1);
}

TEST_F(FlipGeodesicSuite, PiecewiseDijkstraMultipleOverlapAtEnd) {
  std::unique_ptr<ManifoldSurfaceMesh> meshPtr;
  std::unique_ptr<VertexPositionGeometry> geomPtr;
  std::tie(meshPtr, geomPtr) = buildOverlapTestMesh();
  ManifoldSurfaceMesh& mesh = *meshPtr;
  VertexPositionGeometry& geom = *geomPtr;

  std::vector<Vertex> path = {mesh.vertex(0), mesh.vertex(2), mesh.vertex(5), mesh.vertex(1)};
  std::unique_ptr<FlipEdgeNetwork> network =
      FlipEdgeNetwork::constructFromPiecewiseDijkstraPath(mesh, geom, path, false);

  ASSERT_EQ(network->paths.size(), 1);
  EXPECT_EQ(network->paths[0]->getHalfedgeList().size(), 1);
}
