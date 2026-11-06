#include "geometrycentral/surface/boundary_first_flattening.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/numerical/linear_algebra_utilities.h"
#include "geometrycentral/numerical/linear_solvers.h"

#include "gtest/gtest.h"

#include <cmath>
#include <memory>
#include <vector>

using namespace geometrycentral;
using namespace geometrycentral::surface;

// ============================================================
// =============== BFF tests
// ============================================================
//
// A flat unit disk as a polar grid: center vertex, nR rings of nS vertices. The "slit" variant
// is the same disk cut open along the ray theta = 0 (ring vertices on the ray are duplicated),
// so the center is a boundary vertex and the mesh is still a topological disk. Scale factors u
// for one interior cone at the center come from the Yamabe equation on the UNCUT disk with
// u = 0 on the outer circle. Lengths scale by e^{+u}, matching computeBoundaryPositions().

namespace {

const int nR = 20, nS = 64;
const double coneAngle = 0.5; // target curvature at the center

struct Disk {
  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geom;
  std::vector<size_t> source; // vertex -> uncut vertex id
};

size_t uncutId(int k, int j) { return k == 0 ? 0 : 1 + (k - 1) * nS + (j % nS); }

Disk makeDisk(bool slit) {
  const int perRing = slit ? nS + 1 : nS;
  auto id = [&](int k, int j) -> size_t { return k == 0 ? 0 : 1 + (k - 1) * perRing + (slit ? j : j % nS); };
  std::vector<Vector3> pos(1 + nR * perRing);
  Disk d;
  d.source.resize(pos.size());
  pos[0] = Vector3{0, 0, 0};
  for (int k = 1; k <= nR; k++) {
    for (int j = 0; j < perRing; j++) {
      double r = double(k) / nR, t = 2 * M_PI * j / nS;
      pos[id(k, j)] = Vector3{r * std::cos(t), r * std::sin(t), 0};
      d.source[id(k, j)] = uncutId(k, j);
    }
  }
  std::vector<std::vector<size_t>> faces;
  for (int j = 0; j < nS; j++) {
    faces.push_back({id(0, 0), id(1, j), id(1, j + 1)});
    for (int k = 1; k < nR; k++) {
      faces.push_back({id(k, j), id(k + 1, j), id(k + 1, j + 1)});
      faces.push_back({id(k, j), id(k + 1, j + 1), id(k, j + 1)});
    }
  }
  std::tie(d.mesh, d.geom) = makeManifoldSurfaceMeshAndGeometry(faces, pos);
  return d;
}

// Yamabe on the uncut disk: target curvature K~ = Omega + L u (L = PSD cotan Laplacian, lengths ~ e^{+u}),
// K~ = coneAngle at the center and 0 at other interior vertices, u = 0 on the boundary.
Vector<double> coneScaleFactors(Disk& d) {
  d.geom->requireCotanLaplacian();
  d.geom->requireVertexAngleSums();
  size_t n = d.mesh->nVertices();
  Vector<bool> interior(n);
  Vector<double> rhs(n);
  for (Vertex v : d.mesh->vertices()) {
    size_t i = v.getIndex();
    interior(i) = !v.isBoundary();
    double target = (i == 0) ? coneAngle : 0.;
    rhs(i) = v.isBoundary() ? 0. : target - (2 * M_PI - d.geom->vertexAngleSums[v]);
  }
  SparseMatrix<double> L = d.geom->cotanLaplacian;
  BlockDecompositionResult<double> b = blockDecomposeSquare(L, interior);
  Vector<double> rI, rB;
  decomposeVector(b, rhs, rI, rB);
  PositiveDefiniteSolver<double> solver(b.AA);
  Vector<double> uI = solver.solve(rI);
  return reassembleVector(b, uI, Vector<double>(Vector<double>::Zero(rB.size())));
}

// Boundary values in BFF's internal order (boundary vertices in vertex order).
Vector<double> boundaryVector(ManifoldSurfaceMesh& mesh, const VertexData<double>& f) {
  std::vector<double> out;
  for (Vertex v : mesh.vertices())
    if (v.isBoundary()) out.push_back(f[v]);
  Vector<double> r(out.size());
  for (size_t i = 0; i < out.size(); i++) r(i) = out[i];
  return r;
}

double polylineLength(const VertexData<Vector2>& uv, ManifoldSurfaceMesh& mesh, const std::vector<size_t>& ids) {
  double len = 0;
  for (size_t i = 1; i < ids.size(); i++) len += norm(uv[mesh.vertex(ids[i])] - uv[mesh.vertex(ids[i - 1])]);
  return len;
}

class BFFSuite : public ::testing::Test {
protected:
  void SetUp() override {
    uncut = makeDisk(false);
    slit = makeDisk(true);
    Vector<double> u = coneScaleFactors(uncut);
    uSlit = VertexData<double>(*slit.mesh);
    for (Vertex v : slit.mesh->vertices()) uSlit[v] = u(slit.source[v.getIndex()]);
  }
  Disk uncut, slit;
  VertexData<double> uSlit;
};

} // namespace

// Dirichlet-to-Neumann and Neumann-to-Dirichlet must be inverse maps (up to the constant NtD removes).
TEST_F(BFFSuite, DirichletNeumannRoundTrip) {
  BFF bff(*slit.mesh, *slit.geom);
  Vector<double> uB = boundaryVector(*slit.mesh, uSlit);
  Vector<double> back = bff.neumannToDirichlet(bff.dirichletToNeumann(uB));
  Vector<double> expected = uB.array() - uB.mean();
  EXPECT_LT((back - expected).cwiseAbs().maxCoeff(), 1e-8);
}

// Prescribing u, or prescribing the k that u implies, must give the same flattening
// (flattenFromExteriorAngles uses mean-zero u, i.e. a uniform scale of e^{-mean u}).
TEST_F(BFFSuite, ScaleFactorsMatchExteriorAngles) {
  BFF bff(*slit.mesh, *slit.geom);
  Vector<double> uB = boundaryVector(*slit.mesh, uSlit);
  Vector<double> kB = bff.dirichletToNeumann(uB);
  VertexData<double> kData(*slit.mesh, 0.);
  size_t iB = 0;
  for (Vertex v : slit.mesh->vertices())
    if (v.isBoundary()) kData[v] = kB(iB++);

  VertexData<Vector2> fromU = bff.flattenFromScaleFactors(uSlit);
  VertexData<Vector2> fromK = bff.flattenFromExteriorAngles(kData);
  double s = std::exp(uB.mean());
  double maxDiff = 0;
  for (Vertex v : slit.mesh->vertices()) maxDiff = std::max(maxDiff, norm(fromU[v] - s * fromK[v]));
  EXPECT_LT(maxDiff, 1e-8);
}

// The flattened slit disk is a flat cone sector: notch angle = coneAngle, the two sides of the
// slit have equal length, and the outer circle (u = 0) keeps its length.
TEST_F(BFFSuite, ConeFromScaleFactorsIsIsometricOnTheCut) {
  BFF bff(*slit.mesh, *slit.geom);
  VertexData<Vector2> uv = bff.flattenFromScaleFactors(uSlit);
  ManifoldSurfaceMesh& m = *slit.mesh;

  double angleSum = 0;
  for (Corner c : m.vertex(0).adjacentCorners()) {
    Vector2 p = uv[c.vertex()], a = uv[c.halfedge().next().vertex()], b = uv[c.halfedge().next().next().vertex()];
    angleSum += std::abs(angle(a - p, b - p));
  }
  EXPECT_NEAR(2 * M_PI - angleSum, coneAngle, 1e-8);

  const int perRing = nS + 1;
  std::vector<size_t> sideA{0}, sideB{0}, outer;
  for (int k = 1; k <= nR; k++) {
    sideA.push_back(1 + (k - 1) * perRing + 0);
    sideB.push_back(1 + (k - 1) * perRing + nS);
  }
  for (int j = 0; j <= nS; j++) outer.push_back(1 + (nR - 1) * perRing + j);
  double a = polylineLength(uv, m, sideA), b = polylineLength(uv, m, sideB);
  EXPECT_NEAR(a, b, 1e-8 * a);

  double outer3D = 0;
  for (size_t i = 1; i < outer.size(); i++)
    outer3D += norm(slit.geom->vertexPositions[m.vertex(outer[i])] - slit.geom->vertexPositions[m.vertex(outer[i - 1])]);
  EXPECT_NEAR(polylineLength(uv, m, outer), outer3D, 1e-3 * outer3D);
}
