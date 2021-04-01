#include "geometrycentral/pointcloud/point_cloud.h"
#include "geometrycentral/pointcloud/point_cloud_heat_solver.h"
#include "geometrycentral/pointcloud/point_cloud_io.h"
#include "geometrycentral/pointcloud/point_position_frame_geometry.h"
#include "geometrycentral/pointcloud/point_position_geometry.h"
#include "geometrycentral/pointcloud/point_position_normal_geometry.h"
#include "geometrycentral/pointcloud/sample_cloud.h"

#include "gtest/gtest.h"

#include <iostream>
#include <string>
#include <unordered_set>


using namespace geometrycentral;
using namespace geometrycentral::surface;
using namespace geometrycentral::pointcloud;
using std::cout;
using std::endl;


// helpers
namespace {
std::mt19937 mt(42);

double unitRandSeeded() {
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  return dist(mt);
};

// Generate a random point cloud
// NOTE: technically, generating random clouds gives us a tiny probability of a degenerate condition causing some test
// to spontaneously fail. However, this is very unlikely, and the RNG is seeeded, so I think it's alright.
std::tuple<std::unique_ptr<PointCloud>, PointData<Vector3>> generateRandomCloud(size_t nPts) {

  std::unique_ptr<PointCloud> cloud(new PointCloud(nPts));
  PointData<Vector3> pos(*cloud);
  for (Point p : cloud->points()) {
    pos[p] = Vector3{unitRandSeeded(), unitRandSeeded(), unitRandSeeded()};
  }
  return std::make_tuple(std::move(cloud), pos);
}

} // namespace

class PointCloudSuite : public ::testing::Test {};

// ============================================================
// =============== Basic data structure tests
// ============================================================

TEST_F(PointCloudSuite, PointCloudConstructor) {
  PointCloud cloud(256);
  EXPECT_EQ(cloud.nPoints(), 256);
}

TEST_F(PointCloudSuite, IteratePoints) {
  size_t N = 256;
  PointCloud cloud(N);
  size_t sum = 0;
  for (Point p : cloud.points()) {
    EXPECT_GE(p.getIndex(), 0);
    sum += p.getIndex();
  }

  // make sure we touched each once
  EXPECT_EQ(sum, N * (N - 1) / 2);
}

TEST_F(PointCloudSuite, ContainerBasics) {
  size_t N = 256;
  PointCloud cloud(N);

  PointData<int> vals(cloud);
  EXPECT_EQ(vals.size(), N);

  // Store and retrieve some data
  for (Point p : cloud.points()) {
    vals[p] = p.getIndex();
  }
  for (Point p : cloud.points()) {
    EXPECT_EQ(vals[p], p.getIndex());
  }

  // Fill
  vals = PointData<int>(cloud, 77);
  for (Point p : cloud.points()) {
    EXPECT_EQ(vals[p], 77);
  }
}


// ============================================================
// =============== Geometry tests
// ============================================================

TEST_F(PointCloudSuite, PositionGeometryBasic) {
  size_t N = 256;
  std::unique_ptr<PointCloud> cloud;
  PointData<Vector3> pos;
  std::tie(cloud, pos) = generateRandomCloud(N);

  PointPositionGeometry geom(*cloud, pos);

  for (Point p : cloud->points()) {
    EXPECT_GE(geom.positions[p].x, 0.);
  }
}

TEST_F(PointCloudSuite, PositionNormalGeometryBasic) {
  size_t N = 256;
  std::unique_ptr<PointCloud> cloud;
  PointData<Vector3> pos;
  std::tie(cloud, pos) = generateRandomCloud(N);

  // Random normals
  PointData<Vector3> knownNormals(*cloud);
  for (Point p : cloud->points()) {
    knownNormals[p] = normalize(Vector3{unitRandSeeded(), unitRandSeeded(), unitRandSeeded()});
  }

  PointPositionNormalGeometry geom(*cloud, pos, knownNormals);

  geom.requireNormals();
  for (Point p : cloud->points()) {
    EXPECT_EQ(geom.normals[p], knownNormals[p]);
  }
}

TEST_F(PointCloudSuite, PositionFrameGeometryBasic) {
  size_t N = 256;
  std::unique_ptr<PointCloud> cloud;
  PointData<Vector3> pos;
  std::tie(cloud, pos) = generateRandomCloud(N);

  // Random frames
  PointData<Vector3> knownNormals(*cloud);
  PointData<std::array<Vector3, 2>> knownBases(*cloud);
  PointData<std::array<Vector3, 3>> knownFrames(*cloud);
  for (Point p : cloud->points()) {
    knownNormals[p] = normalize(Vector3{unitRandSeeded(), unitRandSeeded(), unitRandSeeded()});
    std::array<Vector3, 2> b = knownNormals[p].buildTangentBasis();
    knownBases[p] = {b[1], -b[0]}; // mix it up so it's not the default basis

    knownFrames[p][0] = knownBases[p][0];
    knownFrames[p][1] = knownBases[p][1];
    knownFrames[p][2] = knownNormals[p];
  }

  PointPositionFrameGeometry geom(*cloud, pos, knownFrames);

  geom.requireNormals();
  geom.requireTangentBasis();
  for (Point p : cloud->points()) {
    EXPECT_EQ(geom.normals[p], knownNormals[p]);
    EXPECT_EQ(geom.tangentBasis[p][0], knownBases[p][0]);
    EXPECT_EQ(geom.tangentBasis[p][1], knownBases[p][1]);
  }
}


TEST_F(PointCloudSuite, GeometryQuantity_PointIndices) {
  // Make the geometry
  size_t N = 256;
  std::unique_ptr<PointCloud> cloud;
  PointData<Vector3> pos;
  std::tie(cloud, pos) = generateRandomCloud(N);
  PointPositionGeometry geom(*cloud, pos);

  geom.requirePointIndices();
  for (Point p : cloud->points()) {
    EXPECT_GE(geom.pointIndices[p], p.getIndex());
  }
}

TEST_F(PointCloudSuite, GeometryQuantity_Neighbors) {
  // Make the geometry
  size_t N = 256;
  std::unique_ptr<PointCloud> cloud;
  PointData<Vector3> pos;
  std::tie(cloud, pos) = generateRandomCloud(N);
  PointPositionGeometry geom(*cloud, pos);

  geom.requireNeighbors();
  for (Point p : cloud->points()) {
    std::vector<Point>& neigh = geom.neighbors->neighbors[p];
    size_t M = neigh.size();
    EXPECT_EQ(M, geom.kNeighborSize);
    std::unordered_set<Point> seenNeigh;
    seenNeigh.insert(p);
    for (size_t iN = 0; iN < M; iN++) {
      Point pN = neigh[iN];
      EXPECT_EQ(seenNeigh.find(pN), seenNeigh.end()); // shouldn't have seen
      seenNeigh.insert(pN);
    }
  }
}

TEST_F(PointCloudSuite, GeometryQuantity_Normals) {
  // Make the geometry
  size_t N = 256;
  std::unique_ptr<PointCloud> cloud;
  PointData<Vector3> pos;
  std::tie(cloud, pos) = generateRandomCloud(N);
  PointPositionGeometry geom(*cloud, pos);

  geom.requireNormals();
  for (Point p : cloud->points()) {
    EXPECT_NEAR(geom.normals[p].norm(), 1.0, 1e-6);
  }
}


TEST_F(PointCloudSuite, GeometryQuantity_Basis) {
  // Make the geometry
  size_t N = 256;
  std::unique_ptr<PointCloud> cloud;
  PointData<Vector3> pos;
  std::tie(cloud, pos) = generateRandomCloud(N);
  PointPositionGeometry geom(*cloud, pos);

  geom.requireTangentBasis();
  geom.requireNormals();
  for (Point p : cloud->points()) {
    Vector3 cval = cross(geom.tangentBasis[p][0], geom.tangentBasis[p][1]);
    EXPECT_NEAR(norm(cval - geom.normals[p]), 0., 1e-6);
  }
}

TEST_F(PointCloudSuite, GeometryQuantity_TangentCoords) {
  // Make the geometry
  size_t N = 256;
  std::unique_ptr<PointCloud> cloud;
  PointData<Vector3> pos;
  std::tie(cloud, pos) = generateRandomCloud(N);
  PointPositionGeometry geom(*cloud, pos);

  geom.requireTangentCoordinates();
  geom.requireNeighbors();
  for (Point p : cloud->points()) {
    std::vector<Vector2>& coords = geom.tangentCoordinates[p];
    std::vector<Point>& neigh = geom.neighbors->neighbors[p];
    EXPECT_EQ(coords.size(), neigh.size());
    for (Vector2 c : coords) {
      EXPECT_TRUE(isfinite(c));
    }
  }
}

TEST_F(PointCloudSuite, GeometryQuantity_TangentTransport) {
  // Make the geometry
  size_t N = 256;
  std::unique_ptr<PointCloud> cloud;
  PointData<Vector3> pos;
  std::tie(cloud, pos) = generateRandomCloud(N);
  PointPositionGeometry geom(*cloud, pos);

  geom.requireTangentTransport();
  geom.requireNeighbors();
  for (Point p : cloud->points()) {
    std::vector<Vector2>& trans = geom.tangentTransport[p];
    std::vector<Point>& neigh = geom.neighbors->neighbors[p];
    EXPECT_EQ(trans.size(), neigh.size());
    for (Vector2 t : trans) {
      EXPECT_NEAR(t.norm(), 1.0, 1e-6);
    }
  }
}

TEST_F(PointCloudSuite, GeometryQuantity_TuftedTriangulation) {
  // Make the geometry
  size_t N = 256;
  std::unique_ptr<PointCloud> cloud;
  PointData<Vector3> pos;
  std::tie(cloud, pos) = generateRandomCloud(N);
  PointPositionGeometry geom(*cloud, pos);

  geom.requireTuftedTriangulation();

  // Check mesh connectivity
  EXPECT_EQ(geom.tuftedMesh->nVertices(), N);
  geom.tuftedMesh->validateConnectivity();
  EXPECT_TRUE(geom.tuftedMesh->isEdgeManifold());

  // Check geometry
  geom.tuftedGeom->requireEdgeCotanWeights();
  for (Edge e : geom.tuftedMesh->edges()) {
    EXPECT_GT(geom.tuftedGeom->edgeCotanWeights[e], -1e-4);
  }
}

TEST_F(PointCloudSuite, GeometryQuantity_Laplacian) {
  // Make the geometry
  size_t N = 256;
  std::unique_ptr<PointCloud> cloud;
  PointData<Vector3> pos;
  std::tie(cloud, pos) = generateRandomCloud(N);
  PointPositionGeometry geom(*cloud, pos);

  geom.requireLaplacian();

  checkFinite(geom.laplacian);
  checkSymmetric(geom.laplacian);

  EXPECT_NEAR(geom.laplacian.sum(), 0., 1e-5);
}

TEST_F(PointCloudSuite, GeometryQuantity_ConnectionLaplacian) {
  // Make the geometry
  size_t N = 256;
  std::unique_ptr<PointCloud> cloud;
  PointData<Vector3> pos;
  std::tie(cloud, pos) = generateRandomCloud(N);
  PointPositionGeometry geom(*cloud, pos);

  geom.requireConnectionLaplacian();

  checkFinite(geom.connectionLaplacian);
  checkHermitian(geom.connectionLaplacian);
}

TEST_F(PointCloudSuite, GeometryQuantity_Gradient) {
  // Make the geometry
  size_t N = 256;
  std::unique_ptr<PointCloud> cloud;
  PointData<Vector3> pos;
  std::tie(cloud, pos) = generateRandomCloud(N);
  PointPositionGeometry geom(*cloud, pos);

  geom.requireGradient();

  checkFinite(geom.gradient);

  // Make sure constant vector is in the kernel
  Vector<std::complex<double>> ones = Vector<std::complex<double>>::Ones(N);
  Vector<std::complex<double>> gradOne = geom.gradient * ones;
  double maxMag = gradOne.array().abs().maxCoeff();

  EXPECT_NEAR(maxMag, 0., 1e-5);
}

// ============================================================
// =============== Algorithm tests
// ============================================================

TEST_F(PointCloudSuite, HeatSolverDistance) {
  size_t N = 256;
  std::unique_ptr<PointCloud> cloud;
  PointData<Vector3> pos;
  std::tie(cloud, pos) = generateRandomCloud(N);
  PointPositionGeometry geom(*cloud, pos);

  Point pSource = cloud->point(7);

  PointCloudHeatSolver solver(*cloud, geom);
  PointData<double> dist = solver.computeDistance(pSource);
  for (Point p : cloud->points()) {
    EXPECT_TRUE(std::isfinite(dist[p]));
    EXPECT_GE(dist[p], 0.);
  }
  EXPECT_EQ(dist[pSource], 0.);
}

TEST_F(PointCloudSuite, HeatSolverScalarExtend) {
  size_t N = 256;
  std::unique_ptr<PointCloud> cloud;
  PointData<Vector3> pos;
  std::tie(cloud, pos) = generateRandomCloud(N);
  PointPositionGeometry geom(*cloud, pos);

  PointCloudHeatSolver solver(*cloud, geom);

  Point pSource1 = cloud->point(7);
  Point pSource2 = cloud->point(12);
  Point pSource3 = cloud->point(13);
  double X1 = 1.;
  double X2 = -12;
  double X3 = 8;
  PointData<double> extend = solver.extendScalars(
      {std::make_tuple(pSource1, X1), std::make_tuple(pSource2, X2), std::make_tuple(pSource3, X3)});
  for (Point p : cloud->points()) {
    EXPECT_TRUE(std::isfinite(extend[p]));
    EXPECT_LT(std::fabs(extend[p]), 12. + 1e-3);
  }

  EXPECT_NEAR(extend[pSource1], X1, 1e-1);
  EXPECT_NEAR(extend[pSource2], X2, 1e-1);
  EXPECT_NEAR(extend[pSource3], X3, 1e-1);
}

TEST_F(PointCloudSuite, HeatSolverTransport) {
  size_t N = 256;
  std::unique_ptr<PointCloud> cloud;
  PointData<Vector3> pos;
  std::tie(cloud, pos) = generateRandomCloud(N);
  PointPositionGeometry geom(*cloud, pos);

  Point pSource1 = cloud->point(7);
  Vector2 X1{1., 0};

  PointCloudHeatSolver solver(*cloud, geom);


  { // Single source
    PointData<Vector2> transport = solver.transportTangentVector(pSource1, X1);
    for (Point p : cloud->points()) {
      EXPECT_TRUE(isfinite(transport[p]));
      EXPECT_NEAR(norm(transport[p]), 1., 1e-5);
    }
    EXPECT_NEAR(norm(transport[pSource1] - X1), 0., 1e-2);
  }

  { // Multiple sources
    Point pSource2 = cloud->point(12);
    Point pSource3 = cloud->point(13);
    Vector2 X2{2., 2.};
    Vector2 X3{0., 4.};
    PointData<Vector2> transport = solver.transportTangentVectors(
        {std::make_tuple(pSource1, X1), std::make_tuple(pSource2, X2), std::make_tuple(pSource3, X3)});
    for (Point p : cloud->points()) {
      EXPECT_TRUE(isfinite(transport[p]));
      EXPECT_GT(norm(transport[p]), 0.);
      EXPECT_LT(norm(transport[p]), 4. + 1e-3);
    }

    EXPECT_NEAR(norm(transport[pSource1] - X1), 0., 1e-1);
    EXPECT_NEAR(norm(transport[pSource2] - X2), 0., 1e-1);
    EXPECT_NEAR(norm(transport[pSource3] - X3), 0., 1e-1);
  }
}

TEST_F(PointCloudSuite, HeatSolverLogmap) {
  size_t N = 256;
  std::unique_ptr<PointCloud> cloud;
  PointData<Vector3> pos;
  std::tie(cloud, pos) = generateRandomCloud(N);
  PointPositionGeometry geom(*cloud, pos);

  Point pSource = cloud->point(7);

  PointCloudHeatSolver solver(*cloud, geom);
  PointData<Vector2> logmap = solver.computeLogMap(pSource);
  for (Point p : cloud->points()) {
    EXPECT_TRUE(isfinite(logmap[p]));
  }
  EXPECT_EQ(norm(logmap[pSource]), 0.);
}

// ============================================================
// =============== Utility tests
// ============================================================

TEST_F(PointCloudSuite, ReadWrite_ply) {
  // Make the geometry
  size_t N = 256;
  std::unique_ptr<PointCloud> cloud;
  PointData<Vector3> pos;
  std::tie(cloud, pos) = generateRandomCloud(N);
  PointPositionGeometry geom(*cloud, pos);

  // Write it out
  std::string filename = "test_cloud.ply";
  writePointCloud(*cloud, geom, filename);

  // Read it back in
  std::unique_ptr<PointCloud> newCloud;
  std::unique_ptr<PointPositionGeometry> newGeom;
  std::tie(newCloud, newGeom) = readPointCloud(filename);

  // Make sure we got what we expected
  EXPECT_EQ(cloud->nPoints(), newCloud->nPoints());

  for (size_t iP = 0; iP < N; iP++) {
    EXPECT_NEAR(geom.positions[iP].x, newGeom->positions[iP].x, 1e-8);
    EXPECT_NEAR(geom.positions[iP].y, newGeom->positions[iP].y, 1e-8);
    EXPECT_NEAR(geom.positions[iP].z, newGeom->positions[iP].z, 1e-8);
  }
}

TEST_F(PointCloudSuite, ReadWrite_obj) {
  // Make the geometry
  size_t N = 256;
  std::unique_ptr<PointCloud> cloud;
  PointData<Vector3> pos;
  std::tie(cloud, pos) = generateRandomCloud(N);
  PointPositionGeometry geom(*cloud, pos);

  // Write it out
  std::string filename = "test_cloud.obj";
  writePointCloud(*cloud, geom, filename);

  // Read it back in
  std::unique_ptr<PointCloud> newCloud;
  std::unique_ptr<PointPositionGeometry> newGeom;
  std::tie(newCloud, newGeom) = readPointCloud(filename);

  // Make sure we got what we expected
  EXPECT_EQ(cloud->nPoints(), newCloud->nPoints());

  for (size_t iP = 0; iP < N; iP++) {
    EXPECT_NEAR(geom.positions[iP].x, newGeom->positions[iP].x, 1e-8);
    EXPECT_NEAR(geom.positions[iP].y, newGeom->positions[iP].y, 1e-8);
    EXPECT_NEAR(geom.positions[iP].z, newGeom->positions[iP].z, 1e-8);
  }
}
