#include "geometrycentral/surface/simple_polygon_mesh.h"

#include "load_test_meshes.h"

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

