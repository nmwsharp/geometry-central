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

// ============================================================
// =============== Closed-loop shortening (termination)
// ============================================================
//
// A contractible closed loop has no nontrivial geodesic representative, so flip-shortening can
// collapse it to a single self-edge and then spin forever unfolding/reshortening it. The fix in
// processSingleEdgeLoop must (a) make such loops TERMINATE, while (b) NOT regressing
// non-contractible loops, which can legitimately pass through a single self-edge on the way to a
// real geodesic (and need the unfold to straighten that last corner).

namespace {

// Deterministic [-1,1] hash (no <random> needed), shared by the mesh builders below.
double flipLoopPseudo(int i) {
  double s = std::sin((double)i * 12.9898) * 43758.5453;
  return 2.0 * (s - std::floor(s)) - 1.0;
}

// Flat triangulated disk with a sharp cone apex at the center. Intrinsically Euclidean away from
// the apex, so EVERY closed loop is contractible. Returns the mesh plus the ring-1 loop edges.
std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
buildConeDisk(int R, int nSeg, double coneH) {
  auto ringVert = [&](int r, int k) -> size_t {
    return 1 + (size_t)(r - 1) * nSeg + (size_t)(((k % nSeg) + nSeg) % nSeg);
  };
  std::vector<Vector3> pos;
  pos.push_back(Vector3{0.0, 0.0, coneH}); // apex
  for (int r = 1; r <= R; r++)
    for (int k = 0; k < nSeg; k++) {
      double ang = 2.0 * M_PI * k / nSeg;
      pos.push_back(Vector3{(double)r * std::cos(ang), (double)r * std::sin(ang), 0.0});
    }
  std::vector<std::vector<size_t>> faces;
  for (int k = 0; k < nSeg; k++) faces.push_back({0, ringVert(1, k), ringVert(1, k + 1)});
  for (int r = 1; r < R; r++)
    for (int k = 0; k < nSeg; k++) {
      size_t a = ringVert(r, k), b = ringVert(r, k + 1), c = ringVert(r + 1, k), d = ringVert(r + 1, k + 1);
      faces.push_back({a, c, d});
      faces.push_back({a, d, b});
    }
  return makeManifoldSurfaceMeshAndGeometry(faces, pos);
}

// Torus, optionally with one vertex spiked radially (a cone point) and the tube pinched. On
// cone-concentrated/pinched tori, a non-contractible ring collapses to a sub-pi self-edge whose
// unfold is required to reach the geodesic -- the case that distinguishes a correct fix from a
// blanket no-op.
std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
buildTorus(int nU, int nV, double Rmaj, double rmin, double jitter, int coneVert, double coneAmt, double pinch) {
  auto vid = [&](int i, int j) { return (size_t)((i % nU) * nV + (j % nV)); };
  std::vector<Vector3> pos(nU * nV);
  for (int i = 0; i < nU; i++)
    for (int j = 0; j < nV; j++) {
      double u = 2 * M_PI * i / nU, v = 2 * M_PI * j / nV;
      double rm = rmin * (1.0 + pinch * std::cos(u));
      double jit = 1.0 + jitter * flipLoopPseudo(i * 131 + j * 17);
      double Rr = (Rmaj + rm * std::cos(v)) * jit;
      Vector3 P{Rr * std::cos(u), Rr * std::sin(u), rm * std::sin(v) * jit};
      size_t id = vid(i, j);
      if ((int)id == coneVert) P *= (1.0 + coneAmt);
      pos[id] = P;
    }
  std::vector<std::vector<size_t>> faces;
  for (int i = 0; i < nU; i++)
    for (int j = 0; j < nV; j++)
      faces.push_back({vid(i, j), vid(i + 1, j), vid(i + 1, j + 1)}),
          faces.push_back({vid(i, j), vid(i + 1, j + 1), vid(i, j + 1)});
  return makeManifoldSurfaceMeshAndGeometry(faces, pos);
}

// Mark the edge va--vb in the path set; returns false if not found.
bool markLoopEdge(ManifoldSurfaceMesh& mesh, EdgeData<bool>& inPath, size_t va, size_t vb) {
  for (Halfedge he : mesh.vertex(va).outgoingHalfedges())
    if (he.tipVertex() == mesh.vertex(vb)) {
      inPath[he.edge()] = true;
      return true;
    }
  return false;
}

} // namespace

// A contractible closed loop must TERMINATE (not spin). We cap iterations generously so that a
// regression cannot hang the test suite; a correct run drains the queue in a handful of steps.
TEST_F(FlipGeodesicSuite, ContractibleClosedLoopTerminates) {
  const int nSeg = 6;
  std::unique_ptr<ManifoldSurfaceMesh> meshPtr;
  std::unique_ptr<VertexPositionGeometry> geomPtr;
  std::tie(meshPtr, geomPtr) = buildConeDisk(/*R=*/2, nSeg, /*coneH=*/4.0);
  ManifoldSurfaceMesh& mesh = *meshPtr;
  VertexPositionGeometry& geom = *geomPtr;

  EdgeData<bool> inPath(mesh, false);
  for (int k = 0; k < nSeg; k++) ASSERT_TRUE(markLoopEdge(mesh, inPath, 1 + k, 1 + (k + 1) % nSeg));
  VertexData<bool> extraMarked(mesh, false);

  std::unique_ptr<FlipEdgeNetwork> network =
      FlipEdgeNetwork::constructFromEdgeSet(mesh, geom, inPath, extraMarked);
  ASSERT_NE(network, nullptr);
  ASSERT_EQ(network->paths.size(), 1);
  ASSERT_TRUE(network->paths[0]->isClosed);
  network->addAllWedgesToAngleQueue(); // queue the closed-loop seam wedge (ctor skips the first segment)

  const size_t cap = 100000;
  network->iterativeShorten(cap);

  // Termination signature: the queue drained and we did NOT exhaust the cap (a spin would do both
  // the opposite). This assertion cannot hang even if the guard is removed -- the cap bounds it.
  EXPECT_TRUE(network->wedgeAngleQueue.empty());
  EXPECT_LT(network->nShortenIters, (long long)cap);
}

// A non-contractible loop on a plain torus must converge to a closed geodesic.
TEST_F(FlipGeodesicSuite, NonContractibleClosedLoopReachesGeodesic) {
  const int nU = 12, nV = 8;
  std::unique_ptr<ManifoldSurfaceMesh> meshPtr;
  std::unique_ptr<VertexPositionGeometry> geomPtr;
  std::tie(meshPtr, geomPtr) = buildTorus(nU, nV, 3.0, 1.0, 0.0, -1, 0.0, 0.0);
  ManifoldSurfaceMesh& mesh = *meshPtr;
  VertexPositionGeometry& geom = *geomPtr;

  EdgeData<bool> inPath(mesh, false);
  auto vid = [&](int i, int j) { return (size_t)((i % nU) * nV + (j % nV)); };
  for (int j = 0; j < nV; j++) ASSERT_TRUE(markLoopEdge(mesh, inPath, vid(0, j), vid(0, j + 1)));
  VertexData<bool> extraMarked(mesh, false);

  std::unique_ptr<FlipEdgeNetwork> network =
      FlipEdgeNetwork::constructFromEdgeSet(mesh, geom, inPath, extraMarked);
  ASSERT_NE(network, nullptr);
  ASSERT_TRUE(network->paths[0]->isClosed);
  network->addAllWedgesToAngleQueue(); // queue the closed-loop seam wedge (ctor skips the first segment)

  network->iterativeShorten(100000);
  EXPECT_TRUE(network->wedgeAngleQueue.empty());
  EXPECT_GE(network->minAngle(), M_PI - network->EPS_ANGLE); // straight => geodesic
}

// The discriminating case: a non-contractible ring on a cone-concentrated, pinched torus whose
// final corner collapses to a sub-pi single self-edge. Reaching the geodesic REQUIRES unfolding
// that self-edge -- so a blanket "never unfold" guard leaves it stuck at a kink (minAngle < pi),
// while the round-trip-progress guard straightens it. This test fails on the naive fix.
TEST_F(FlipGeodesicSuite, NonContractibleLoopThroughSelfEdgeReachesGeodesic) {
  const int nU = 4, nV = 4;
  std::unique_ptr<ManifoldSurfaceMesh> meshPtr;
  std::unique_ptr<VertexPositionGeometry> geomPtr;
  std::tie(meshPtr, geomPtr) = buildTorus(nU, nV, 2.0, 1.2, 0.8, /*coneVert=*/6, /*coneAmt=*/3.5, /*pinch=*/0.7);
  ManifoldSurfaceMesh& mesh = *meshPtr;
  VertexPositionGeometry& geom = *geomPtr;

  EdgeData<bool> inPath(mesh, false);
  auto vid = [&](int i, int j) { return (size_t)((i % nU) * nV + (j % nV)); };
  for (int j = 0; j < nV; j++) ASSERT_TRUE(markLoopEdge(mesh, inPath, vid(0, j), vid(0, j + 1)));
  VertexData<bool> extraMarked(mesh, false);

  std::unique_ptr<FlipEdgeNetwork> network =
      FlipEdgeNetwork::constructFromEdgeSet(mesh, geom, inPath, extraMarked);
  ASSERT_NE(network, nullptr);
  ASSERT_TRUE(network->paths[0]->isClosed);
  network->addAllWedgesToAngleQueue(); // queue the closed-loop seam wedge (ctor skips the first segment)

  network->iterativeShorten(100000);
  EXPECT_TRUE(network->wedgeAngleQueue.empty());
  EXPECT_GE(network->minAngle(), M_PI - network->EPS_ANGLE);
}
