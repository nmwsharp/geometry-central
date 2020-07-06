#include "geometrycentral/surface/simple_polygon_mesh.h"

#include "load_test_meshes.h"

#include "happly.h"

#include "gtest/gtest.h"

#include <iostream>
#include <string>
#include <unordered_set>


using namespace geometrycentral;
using namespace geometrycentral::surface;
using std::cout;
using std::endl;

class SimplePolygonSuite : public MeshAssetSuite {};

// ============================================================
// =============== SimplePolygonMesh tests
// ============================================================

TEST_F(SimplePolygonSuite, BasicTest) {
  for (MeshAsset& a : allMeshes()) {
    a.printThyName();

    SimplePolygonMesh simpleMesh(a.sourcePath);

    ASSERT_GT(simpleMesh.nFaces(), 0);
    ASSERT_GT(simpleMesh.nVertices(), 0);
  }
}


TEST_F(SimplePolygonSuite, PLYLoad) {

  std::string fullPath = std::string(GC_TEST_ASSETS_ABS_PATH) + "/" + "spot.ply";
  happly::PLYData plyIn(fullPath, true);
  plyIn.validate();

  std::vector<std::array<double, 3>> vPos = plyIn.getVertexPositions();
  std::vector<std::vector<size_t>> fInd = plyIn.getFaceIndices();

  for (std::vector<size_t>& face : fInd) {
    EXPECT_EQ(face.size(), 3);
  }
}

