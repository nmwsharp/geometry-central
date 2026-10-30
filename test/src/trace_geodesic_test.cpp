// traceGeodesicBetweenVertices: a segment between two vertices that passes exactly through a vertex of a
// zero-area face, or starts in such a face's corner of angle pi, ends at its end vertex with its length;
// on a clean mesh the result is the trimmed face-by-face trace, unchanged; FlipEdgeNetwork::getPathPolyline through
// a zero-area face runs from the path's first vertex to its last.

#include "geometrycentral/surface/flip_geodesics.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/trace_geodesic.h"

#include "gtest/gtest.h"

#include <cmath>
#include <random>

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace {

// A planar patch with a zero-area face: M = (1, 0) lies on the edge A-B, the cap (A, B, M) has a corner of
// angle pi at M; C above and D below. Vertices A 0, M 1, B 2, C 3, D 4.
std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>> capPatch() {
  std::vector<Vector3> verts = {{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {1, 1, 0}, {1, -1, 0}};
  std::vector<std::vector<size_t>> faces = {{0, 1, 3}, {1, 2, 3}, {0, 2, 1}, {2, 0, 4}};
  return makeManifoldSurfaceMeshAndGeometry(faces, verts);
}

Halfedge outgoing(Vertex from, Vertex to) {
  for (Halfedge he : from.outgoingHalfedges())
    if (he.tipVertex() == to) return he;
  return Halfedge();
}

void expectSegment(VertexPositionGeometry& geom, const std::vector<SurfacePoint>& path, Vertex a, Vertex b,
                   double length) {
  ASSERT_GE(path.size(), 2u);
  EXPECT_TRUE(path.front().type == SurfacePointType::Vertex && path.front().vertex == a);
  EXPECT_TRUE(path.back().type == SurfacePointType::Vertex && path.back().vertex == b);
  const double L = surfacePathLength(geom, path); // NaN if two consecutive points share no face
  EXPECT_NEAR(L, length, 1e-12 * length);
  // the polyline is straight: every point on the segment a -> b
  const Vector3 pa = geom.inputVertexPositions[a], pb = geom.inputVertexPositions[b];
  for (const SurfacePoint& p : path) {
    const Vector3 x = p.interpolate(geom.inputVertexPositions);
    EXPECT_NEAR(norm(cross(x - pa, pb - pa)) / norm(pb - pa), 0., 1e-12);
  }
}

} // namespace

TEST(TraceGeodesicBetweenVertices, ThroughVertexOfZeroAreaFace) {
  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geom;
  std::tie(mesh, geom) = capPatch();
  geom->requireHalfedgeVectorsInVertex();
  Vertex C = mesh->vertex(3), D = mesh->vertex(4);
  // C -> D runs straight down through M
  Vector2 dir = geom->halfedgeVectorsInVertex[outgoing(C, mesh->vertex(1))].normalize();
  std::vector<SurfacePoint> path = traceGeodesicBetweenVertices(*geom, C, 2. * dir, D);
  expectSegment(*geom, path, C, D, 2.);
}

TEST(TraceGeodesicBetweenVertices, FromCornerOfAnglePi) {
  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geom;
  std::tie(mesh, geom) = capPatch();
  geom->requireHalfedgeVectorsInVertex();
  geom->requireCornerAngles();
  Vertex M = mesh->vertex(1), D = mesh->vertex(4);
  Halfedge ma = outgoing(M, mesh->vertex(0)); // the cap's corner at M spans from M -> A to M -> B
  ASSERT_NEAR(geom->cornerAngles[ma.corner()], M_PI, 1e-12);
  // M -> D leaves through the middle of that corner; M's angle sum is 2 pi, so the tangent space is unscaled
  Vector2 dir = geom->halfedgeVectorsInVertex[ma].normalize() * Vector2::fromAngle(M_PI / 2.);
  std::vector<SurfacePoint> path = traceGeodesicBetweenVertices(*geom, M, dir, D);
  expectSegment(*geom, path, M, D, 1.);
}

TEST(TraceGeodesicBetweenVertices, CleanMeshUnchanged) {
  // FlipOut paths on the bunny: each intrinsic edge traced as before (traceGeodesic, trimTraceResult, end
  // vertex appended), bit for bit
  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geom;
  std::tie(mesh, geom) = readManifoldSurfaceMesh(std::string(GC_TEST_ASSETS_ABS_PATH) + "/bob_small.ply");
  std::mt19937 rng(7);
  std::uniform_int_distribution<size_t> pick(0, mesh->nVertices() - 1);
  size_t traced = 0;
  for (int q = 0; q < 10; q++) {
    Vertex a = mesh->vertex(pick(rng)), b = mesh->vertex(pick(rng));
    if (a == b) continue;
    auto net = FlipEdgeNetwork::constructFromDijkstraPath(*mesh, *geom, a, b);
    net->iterativeShorten();
    for (Halfedge he : net->paths[0]->getHalfedgeList()) {
      if (net->tri->edgeIsOriginal[he.edge()]) continue;
      Vertex va = net->tri->vertexLocations[he.vertex()].vertex, vb = net->tri->vertexLocations[he.twin().vertex()].vertex;
      // SignpostIntrinsicTriangulation::halfedgeVector
      net->tri->requireVertexAngleSums();
      const double scale = net->tri->vertexAngleSums[he.vertex()] / (he.vertex().isBoundary() ? M_PI : 2. * M_PI);
      Vector2 vec = Vector2::fromAngle(net->tri->signpostAngle[he] / scale) * net->tri->edgeLengths[he.edge()];
      TraceOptions options;
      options.includePath = true;
      options.maxIters = mesh->nFaces() * 10;
      TraceGeodesicResult old = traceGeodesic(*geom, SurfacePoint(va), vec, options);
      ASSERT_TRUE(trimTraceResult(old, vb));
      old.pathPoints.push_back(SurfacePoint(vb));
      std::vector<SurfacePoint> now = traceGeodesicBetweenVertices(*geom, va, vec, vb, options);
      ASSERT_EQ(now.size(), old.pathPoints.size());
      for (size_t k = 0; k < now.size(); k++) {
        EXPECT_TRUE(now[k] == old.pathPoints[k]);
        if (now[k].type == SurfacePointType::Edge) EXPECT_EQ(now[k].tEdge, old.pathPoints[k].tEdge);
      }
      traced++;
    }
  }
  EXPECT_GT(traced, 10u);
}

TEST(TraceGeodesicBetweenVertices, FlipOutPolylineThroughZeroAreaFace) {
  // FlipOut's shortest path C -> D is the straight segment through M; its polyline must run from C to D on the
  // surface, every step inside one face, with the path's length
  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geom;
  std::tie(mesh, geom) = capPatch();
  Vertex C = mesh->vertex(3), D = mesh->vertex(4);
  auto net = FlipEdgeNetwork::constructFromDijkstraPath(*mesh, *geom, C, D);
  net->iterativeShorten();
  std::vector<std::vector<SurfacePoint>> polyline = net->getPathPolyline();
  ASSERT_EQ(polyline.size(), 1u);
  EXPECT_NEAR(net->length(), 2., 1e-6); // FlipOut stops within its angle tolerance of straight
  expectSegment(*geom, polyline[0], C, D, 2.);
}
