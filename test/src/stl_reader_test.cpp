// A few simple tests for reading STL ascii and binary
// files
#include <memory>

#include "geometrycentral/surface/meshio.h"
#include "gtest/gtest.h"


using namespace geometrycentral;
using namespace geometrycentral::surface;


TEST(StlReaderTests, LoadAsciiStlFile) {
  std::string fullPath = std::string(GC_TEST_ASSETS_ABS_PATH) + "/" + "stl_box_ascii.stl";

  std::unique_ptr<geometrycentral::surface::SurfaceMesh> mesh;
  geometrycentral::surface::ManifoldSurfaceMesh* manifoldMesh = nullptr;
  std::unique_ptr<geometrycentral::surface::VertexPositionGeometry> geometry;
  std::tie(mesh, geometry) = readSurfaceMesh(fullPath);
  ASSERT_EQ(mesh->nVertices(), 8);
  ASSERT_EQ(mesh->nFaces(), 12);
}

TEST(StlReaderTests, LoadBinaryStlFile) {
  std::string fullPath = std::string(GC_TEST_ASSETS_ABS_PATH) + "/" + "stl_box_binary.stl";

  std::unique_ptr<geometrycentral::surface::SurfaceMesh> mesh;
  geometrycentral::surface::ManifoldSurfaceMesh* manifoldMesh = nullptr;
  std::unique_ptr<geometrycentral::surface::VertexPositionGeometry> geometry;
  std::tie(mesh, geometry) = readSurfaceMesh(fullPath);
  ASSERT_EQ(mesh->nVertices(), 8);
  ASSERT_EQ(mesh->nFaces(), 12);
}

