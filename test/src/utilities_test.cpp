#include "geometrycentral/utilities/knn.h"
#include "geometrycentral/utilities/vector3.h"

#include "gtest/gtest.h"

using namespace geometrycentral;

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
