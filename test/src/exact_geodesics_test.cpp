#include "geometrycentral/surface/exact_geodesics.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "load_test_meshes.h"

#include "gtest/gtest.h"

#include <algorithm>
#include <vector>

using namespace geometrycentral;
using namespace geometrycentral::surface;

class ExactGeodesicSuite : public MeshAssetSuite {};

namespace {

// A distance that about half of the vertices are closer than
double medianDistance(SurfaceMesh& mesh, const VertexData<double>& dist) {
  std::vector<double> d;
  for (Vertex v : mesh.vertices()) d.push_back(dist[v]);
  std::nth_element(d.begin(), d.begin() + d.size() / 2, d.end());
  return d[d.size() / 2];
}

} // namespace

// The distance field does not depend on whether the caller holds the curvatures used by the saddle test
TEST_F(ExactGeodesicSuite, ResultIndependentOfRequiredQuantities) {
  for (const MeshAsset& a : {getAsset("spot.ply", true), getAsset("bob_small.ply", true)}) {
    a.printThyName();
    SurfaceMesh& mesh = *a.mesh;
    VertexPositionGeometry& geometry = *a.geometry;
    Vertex source = mesh.vertex(0);

    GeodesicAlgorithmExact mmp(mesh, geometry);
    mmp.propagate(source);
    VertexData<double> dist = mmp.getDistanceFunction();

    geometry.requireVertexGaussianCurvatures();
    GeodesicAlgorithmExact mmpHeld(mesh, geometry);
    mmpHeld.propagate(source);
    VertexData<double> distHeld = mmpHeld.getDistanceFunction();
    geometry.unrequireVertexGaussianCurvatures();

    for (Vertex v : mesh.vertices()) {
      EXPECT_EQ(dist[v], distHeld[v]);
    }
  }
}

// With a distance bound and no stop points, propagation stops at the bound: vertices closer than the bound
// get their exact distance, and some vertices beyond it are left unreached
TEST_F(ExactGeodesicSuite, MaxPropagationDistanceWithoutStopPoints) {
  for (const MeshAsset& a : {getAsset("spot.ply", true), getAsset("bob_small.ply", true)}) {
    a.printThyName();
    SurfaceMesh& mesh = *a.mesh;
    VertexPositionGeometry& geometry = *a.geometry;
    Vertex source = mesh.vertex(0);

    GeodesicAlgorithmExact full(mesh, geometry);
    full.propagate(source);
    VertexData<double> ref = full.getDistanceFunction();
    double bound = medianDistance(mesh, ref);

    GeodesicAlgorithmExact mmp(mesh, geometry);
    mmp.propagate(source, bound);
    VertexData<double> dist = mmp.getDistanceFunction();

    size_t unreached = 0;
    for (Vertex v : mesh.vertices()) {
      if (ref[v] < bound) {
        EXPECT_NEAR(dist[v], ref[v], 1e-12 * bound);
      }
      if (dist[v] >= GEODESIC_INF) unreached++;
    }
    EXPECT_GT(unreached, 0u);
  }
}

// With stop points, propagation stops once they are covered; with a distance bound too, once both hold
TEST_F(ExactGeodesicSuite, MaxPropagationDistanceWithStopPoints) {
  for (const MeshAsset& a : {getAsset("spot.ply", true), getAsset("bob_small.ply", true)}) {
    a.printThyName();
    SurfaceMesh& mesh = *a.mesh;
    VertexPositionGeometry& geometry = *a.geometry;
    Vertex source = mesh.vertex(0);

    GeodesicAlgorithmExact full(mesh, geometry);
    full.propagate(source);
    VertexData<double> ref = full.getDistanceFunction();
    double bound = medianDistance(mesh, ref);

    // A stop vertex well inside the bound
    Vertex stop;
    for (Vertex v : mesh.vertices()) {
      if (v != source && ref[v] < 0.5 * bound && (stop == Vertex() || ref[v] > ref[stop])) stop = v;
    }
    ASSERT_TRUE(stop != Vertex());

    GeodesicAlgorithmExact mmpStop(mesh, geometry);
    mmpStop.propagate(source, GEODESIC_INF, {stop});
    EXPECT_NEAR(mmpStop.getDistance(stop), ref[stop], 1e-12 * bound);

    GeodesicAlgorithmExact mmpBoth(mesh, geometry);
    mmpBoth.propagate(source, bound, {stop});
    VertexData<double> dist = mmpBoth.getDistanceFunction();
    for (Vertex v : mesh.vertices()) {
      if (ref[v] < bound) {
        EXPECT_NEAR(dist[v], ref[v], 1e-12 * bound);
      }
    }
  }
}
