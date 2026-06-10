#include "geometrycentral/surface/marching_triangles.h"
#include "geometrycentral/surface/surface_mesh_factories.h"

#include "gtest/gtest.h"

#include <memory>
#include <tuple>
#include <vector>

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace {

std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>> buildSingleTriangleMesh() {
  std::vector<Vector3> verts = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
  std::vector<std::vector<size_t>> faces = {{0, 1, 2}};
  return makeManifoldSurfaceMeshAndGeometry(faces, verts);
}

std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>> buildTwoTriangleMesh() {
  std::vector<Vector3> verts = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};
  std::vector<std::vector<size_t>> faces = {{0, 1, 2}, {0, 2, 3}};
  return makeManifoldSurfaceMeshAndGeometry(faces, verts);
}

} // namespace

TEST(MarchingTrianglesSuite, ExactIsoVertexCreatesSegment) {
  std::unique_ptr<ManifoldSurfaceMesh> meshPtr;
  std::unique_ptr<VertexPositionGeometry> geomPtr;
  std::tie(meshPtr, geomPtr) = buildSingleTriangleMesh();

  ManifoldSurfaceMesh& mesh = *meshPtr;
  VertexData<double> values(mesh);
  values[mesh.vertex(0)] = 0.;
  values[mesh.vertex(1)] = 1.;
  values[mesh.vertex(2)] = -1.;

  std::vector<std::vector<SurfacePoint>> curves = marchingTriangles(*geomPtr, values, 0.);

  ASSERT_EQ(curves.size(), 1);
  ASSERT_EQ(curves[0].size(), 2);

  bool hasIsoVertex = false;
  bool hasCrossing = false;
  for (const SurfacePoint& p : curves[0]) {
    if (p.type == SurfacePointType::Vertex && p.vertex == mesh.vertex(0)) {
      hasIsoVertex = true;
    }
    if (p.type == SurfacePointType::Edge) {
      Vertex vA = p.edge.firstVertex();
      Vertex vB = p.edge.secondVertex();
      bool onOppositeEdge = (vA == mesh.vertex(1) && vB == mesh.vertex(2)) ||
                            (vA == mesh.vertex(2) && vB == mesh.vertex(1));
      if (onOppositeEdge) {
        hasCrossing = true;
        EXPECT_NEAR(p.tEdge, 0.5, 1e-12);
      }
    }
  }

  EXPECT_TRUE(hasIsoVertex);
  EXPECT_TRUE(hasCrossing);
}

TEST(MarchingTrianglesSuite, SharedEdgeCrossingsStitchIntoOneCurve) {
  std::unique_ptr<ManifoldSurfaceMesh> meshPtr;
  std::unique_ptr<VertexPositionGeometry> geomPtr;
  std::tie(meshPtr, geomPtr) = buildTwoTriangleMesh();

  ManifoldSurfaceMesh& mesh = *meshPtr;
  VertexData<double> values(mesh);
  values[mesh.vertex(0)] = -0.8929741240814866;
  values[mesh.vertex(1)] = -0.8929741240814866;
  values[mesh.vertex(2)] = 0.1403065708384703;
  values[mesh.vertex(3)] = 0.1403065708384703;

  std::vector<std::vector<SurfacePoint>> curves = marchingTriangles(*geomPtr, values, -0.891336467690912);

  ASSERT_EQ(curves.size(), 1);
  EXPECT_EQ(curves[0].size(), 3);
}

TEST(MarchingTrianglesSuite, FullyIsoValuedFaceIsSkipped) {
  std::unique_ptr<ManifoldSurfaceMesh> meshPtr;
  std::unique_ptr<VertexPositionGeometry> geomPtr;
  std::tie(meshPtr, geomPtr) = buildSingleTriangleMesh();

  ManifoldSurfaceMesh& mesh = *meshPtr;
  VertexData<double> values(mesh, 0.);

  std::vector<std::vector<SurfacePoint>> curves = marchingTriangles(*geomPtr, values, 0.);

  EXPECT_TRUE(curves.empty());
}
