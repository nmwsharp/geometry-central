// FlipOutContext (tuned FlipOut) must return the shipped FlipOut's geodesics.

#include "geometrycentral/surface/flip_geodesics.h"
#include "geometrycentral/surface/flip_geodesics_context.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "load_test_meshes.h"

#include "gtest/gtest.h"

#include <random>

using namespace geometrycentral;
using namespace geometrycentral::surface;

class FlipOutContextSuite : public MeshAssetSuite {};

namespace {

struct Shipped {
  double length;
  std::vector<uint32_t> halfedges;
  std::vector<SurfacePoint> polyline;
};

Shipped shippedFlipOut(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, Vertex a, Vertex b, bool polyline) {
  Shipped r;
  auto net = FlipEdgeNetwork::constructFromDijkstraPath(mesh, geom, a, b);
  net->iterativeShorten();
  r.length = net->length();
  for (Halfedge he : net->paths[0]->getHalfedgeList()) r.halfedges.push_back((uint32_t)he.getIndex());
  if (polyline) r.polyline = net->getPathPolyline().front();
  return r;
}

std::vector<std::pair<size_t, size_t>> randomPairs(size_t nV, size_t n, unsigned seed) {
  std::mt19937 rng(seed);
  std::uniform_int_distribution<size_t> dist(0, nV - 1);
  std::vector<std::pair<size_t, size_t>> pairs;
  while (pairs.size() < n) {
    size_t v0 = dist(rng), v1 = dist(rng);
    if (v0 != v1) pairs.emplace_back(v0, v1);
  }
  return pairs;
}

} // namespace

// Same lengths (to rounding of the sum) and the same intrinsic halfedge sequence as the shipped
// path, on closed meshes and meshes with boundary
TEST_F(FlipOutContextSuite, MatchesShippedFlipOut) {
  bool sawBoundary = false;
  for (const MeshAsset& a : {getAsset("spot.ply", true), getAsset("bob_small.ply", true), getAsset("fox.ply", true),
                             getAsset("cat_head.obj", true)}) {
    a.printThyName();
    ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
    VertexPositionGeometry& geom = *a.geometry;
    sawBoundary = sawBoundary || mesh.hasBoundary();

    FlipOutContext ctx(mesh, geom);
    size_t nFlips = 0;
    for (auto& pr : randomPairs(mesh.nVertices(), 40, 7)) {
      Vertex v0 = mesh.vertex(pr.first), v1 = mesh.vertex(pr.second);
      Shipped s = shippedFlipOut(mesh, geom, v0, v1, false);
      double L = ctx.geodesicLength(v0, v1);
      EXPECT_NEAR(L, s.length, 1e-12 * s.length);
      EXPECT_EQ(ctx.lastPathHalfedges(), s.halfedges);
      nFlips += s.halfedges.size();
    }
    EXPECT_GT(nFlips, 0u);
  }
  EXPECT_TRUE(sawBoundary) << "the test set should include a mesh with boundary";
}

// Queries do not leave traces: a repeated query, and a fresh context, give the same answer
TEST_F(FlipOutContextSuite, RewindRestoresTriangulation) {
  auto asset = getAsset("spot.ply", true);
  ManifoldSurfaceMesh& mesh = *asset.manifoldMesh;
  VertexPositionGeometry& geom = *asset.geometry;

  auto pairs = randomPairs(mesh.nVertices(), 30, 11);
  FlipOutContext ctx(mesh, geom);
  std::vector<double> first;
  std::vector<std::vector<uint32_t>> firstPaths;
  for (auto& pr : pairs) {
    first.push_back(ctx.geodesicLength(mesh.vertex(pr.first), mesh.vertex(pr.second)));
    firstPaths.push_back(ctx.lastPathHalfedges());
  }
  EXPECT_GT(ctx.nFlips, 0);
  for (size_t i = 0; i < pairs.size(); i++) {
    double again = ctx.geodesicLength(mesh.vertex(pairs[i].first), mesh.vertex(pairs[i].second));
    EXPECT_EQ(again, first[i]);
    EXPECT_EQ(ctx.lastPathHalfedges(), firstPaths[i]);
    FlipOutContext fresh(mesh, geom);
    EXPECT_EQ(fresh.geodesicLength(mesh.vertex(pairs[i].first), mesh.vertex(pairs[i].second)), first[i]);
  }
}

// Degenerate inputs
TEST_F(FlipOutContextSuite, SameVertexAndNeighbors) {
  auto asset = getAsset("spot.ply", true);
  ManifoldSurfaceMesh& mesh = *asset.manifoldMesh;
  VertexPositionGeometry& geom = *asset.geometry;
  FlipOutContext ctx(mesh, geom);

  Vertex v = mesh.vertex(42);
  EXPECT_EQ(ctx.geodesicLength(v, v), 0.);
  EXPECT_TRUE(ctx.lastPathHalfedges().empty());

  geom.requireEdgeLengths();
  for (Halfedge he : v.outgoingHalfedges()) {
    double L = ctx.geodesicLength(v, he.twin().vertex());
    EXPECT_EQ(L, geom.edgeLengths[he.edge()]);
    ASSERT_EQ(ctx.lastPathHalfedges().size(), 1u);
    EXPECT_EQ(ctx.lastPathHalfedges()[0], (uint32_t)he.getIndex());
  }
  geom.unrequireEdgeLengths();
}

// The traced path equals the shipped polyline
TEST_F(FlipOutContextSuite, PolylineMatchesShipped) {
  for (const MeshAsset& a : {getAsset("spot.ply", true), getAsset("cat_head.obj", true)}) {
    a.printThyName();
    ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
    VertexPositionGeometry& geom = *a.geometry;
    FlipOutContext ctx(mesh, geom);
    size_t nTraced = 0;
    for (auto& pr : randomPairs(mesh.nVertices(), 12, 3)) {
      Vertex v0 = mesh.vertex(pr.first), v1 = mesh.vertex(pr.second);
      Shipped s = shippedFlipOut(mesh, geom, v0, v1, true);
      std::vector<SurfacePoint> pts;
      double L = ctx.geodesic(v0, v1, pts);
      EXPECT_NEAR(L, s.length, 1e-12 * s.length);
      ASSERT_EQ(pts.size(), s.polyline.size());
      EXPECT_EQ(pts.front(), SurfacePoint(v0));
      EXPECT_EQ(pts.back(), SurfacePoint(v1));
      for (size_t i = 0; i < pts.size(); i++) {
        EXPECT_EQ(pts[i].type, s.polyline[i].type);
        Vector3 p = pts[i].interpolate(geom.inputVertexPositions);
        Vector3 q = s.polyline[i].interpolate(geom.inputVertexPositions);
        EXPECT_LT((p - q).norm(), 1e-9);
      }
      nTraced += pts.size();
    }
    EXPECT_GT(nTraced, 0u);
  }
}
