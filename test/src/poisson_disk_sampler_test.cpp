
#include <gtest/gtest.h>

#include <geometrycentral/surface/poisson_disk_sampler.h>
#include <geometrycentral/surface/simple_polygon_mesh.h>
#include <geometrycentral/surface/surface_mesh_factories.h>


using namespace geometrycentral;
using namespace geometrycentral::surface;


size_t sampleSquareDisk(double width, double sampling_distance) {
  double const PTS[4][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}};
  unsigned const TRIS[2][3] = {{0, 1, 2}, {0, 2, 3}};

  SimplePolygonMesh simpleMesh;
  for (const auto& t : TRIS) simpleMesh.polygons.push_back({t[0], t[1], t[2]});

  for (const auto& p : PTS) simpleMesh.vertexCoordinates.push_back({width * p[0], width * p[1], width * p[2]});

  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geometry;
  std::tie(mesh, geometry) = makeManifoldSurfaceMeshAndGeometry(simpleMesh.polygons, simpleMesh.vertexCoordinates);

  PoissonDiskOptions options;
  options.minDist = sampling_distance;

  // make tests reproducible
  geometrycentral::util_mersenne_twister.seed(101);

  PoissonDiskSampler sampler(*mesh, *geometry);
  auto samples = sampler.sample(options);
  return samples.size();
}


class PoissonDiskSamplerSuite : public ::testing::Test {};

TEST_F(PoissonDiskSamplerSuite, PoissonDiskSamplerConstructor) {
  // PoissonDiskSampler doesn't allow to set a random seed.
  // To prevent flaky test failures we 'average' over 10 iterations.
  size_t n1 = 0, n2 = 0;
  for (size_t iter = 0; iter < 10; iter++) {
    n1 += sampleSquareDisk(1.0, 0.1);
    n2 += sampleSquareDisk(100.0, 100.0 * 0.1);
  }

  EXPECT_GT(n1, 600);
  EXPECT_LT(n1, 720);
  EXPECT_NEAR(n1, n2, 100);
}
