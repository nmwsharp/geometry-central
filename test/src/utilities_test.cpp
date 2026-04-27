#include "geometrycentral/utilities/knn.h"
#include "geometrycentral/utilities/vector3.h"

#include "geometrycentral/surface/mesh_ray_tracer.h"
#include "geometrycentral/surface/surface_mesh_factories.h"

#include "gtest/gtest.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

// ============================================================
// =============== NearestNeighborFinder tests
// ============================================================

TEST(UtilitiesSuite, KNNReturnsClosestPoint) {
  std::vector<Vector3> points = {{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {10, 0, 0}};
  NearestNeighborFinder finder(points);

  std::vector<size_t> result = finder.kNearest({1.1, 0, 0}, 1);
  ASSERT_EQ(result.size(), 1);
  EXPECT_EQ(result[0], 1); // closest to {1,0,0}
}

TEST(UtilitiesSuite, KNNNeighborsExcludesSource) {
  std::vector<Vector3> points = {{0, 0, 0}, {1, 0, 0}, {2, 0, 0}};
  NearestNeighborFinder finder(points);

  std::vector<size_t> result = finder.kNearestNeighbors(1, 1);
  ASSERT_EQ(result.size(), 1);
  EXPECT_NE(result[0], 1); // source index should not appear
}

TEST(UtilitiesSuite, KNNRadiusSearch) {
  std::vector<Vector3> points = {{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {10, 0, 0}};
  NearestNeighborFinder finder(points);

  std::vector<size_t> result = finder.radiusSearch({0, 0, 0}, 1.5);
  EXPECT_EQ(result.size(), 2); // {0,0,0} and {1,0,0} are within radius 1.5
}

// ============================================================
// =============== MeshRayTracer tests
// ============================================================

// Single triangle in the z=0 plane: vertices (0,0,0), (2,0,0), (0,2,0)
namespace {
std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>> buildRayTracerTestMesh() {
  std::vector<Vector3> verts = {{0, 0, 0}, {2, 0, 0}, {0, 2, 0}};
  std::vector<std::vector<size_t>> faces = {{0, 1, 2}};
  return makeManifoldSurfaceMeshAndGeometry(faces, verts);
}
} // namespace

TEST(UtilitiesSuite, RayTracerHit) {
  std::unique_ptr<ManifoldSurfaceMesh> meshPtr;
  std::unique_ptr<VertexPositionGeometry> geomPtr;
  std::tie(meshPtr, geomPtr) = buildRayTracerTestMesh();

  MeshRayTracer tracer(*geomPtr);
  RayHitResult result = tracer.trace({0.5, 0.5, -1}, {0, 0, 1});

  ASSERT_TRUE(result.hit);
  EXPECT_NEAR(result.t, 1.0, 1e-6);
  EXPECT_EQ(result.faceIndex, 0);
}

TEST(UtilitiesSuite, RayTracerMiss) {
  std::unique_ptr<ManifoldSurfaceMesh> meshPtr;
  std::unique_ptr<VertexPositionGeometry> geomPtr;
  std::tie(meshPtr, geomPtr) = buildRayTracerTestMesh();

  MeshRayTracer tracer(*geomPtr);
  RayHitResult result = tracer.trace({5, 5, -1}, {0, 0, 1});

  EXPECT_FALSE(result.hit);
}

TEST(UtilitiesSuite, RayTracerBarycentrics) {
  std::unique_ptr<ManifoldSurfaceMesh> meshPtr;
  std::unique_ptr<VertexPositionGeometry> geomPtr;
  std::tie(meshPtr, geomPtr) = buildRayTracerTestMesh();

  MeshRayTracer tracer(*geomPtr);
  // Ray aimed at vertex 0 (origin) should return faceCoords near (1,0,0)
  RayHitResult result = tracer.trace({0.01, 0.01, -1}, {0, 0, 1});

  ASSERT_TRUE(result.hit);
  double sumBary = result.faceCoords.x + result.faceCoords.y + result.faceCoords.z;
  EXPECT_NEAR(sumBary, 1.0, 1e-6);
}
