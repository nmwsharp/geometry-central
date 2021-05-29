

#include <gtest/gtest.h>
#include <iostream>

#include "geometrycentral/surface/rich_surface_mesh_data.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/utilities/eigen_interop_helpers.h"

#include "load_test_meshes.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

class MappingUtilityTest : public MeshAssetSuite {};

TEST_F(MappingUtilityTest, MappingTest) {
  for (MeshAsset& a : allMeshes(true)) {
    SurfaceMesh& mesh = *a.mesh;
    VertexPositionGeometry& geom = *a.geometry;

    // Test Map N x 1 -> N x 3
    auto pos1 = EigenMap<double, 3>(geom.vertexPositions);
    ASSERT_EQ(geom.vertexPositions.size(), pos1.rows());
    ASSERT_EQ(3, pos1.cols());

    // Test Flattened Map N x 1 -> 3N x 1
    auto pos2 = FlattenedEigenMap<double, 3>(geom.vertexPositions);
    ASSERT_EQ(3 * geom.vertexPositions.size(), pos2.rows());
    ASSERT_EQ(1, pos2.cols());

    for (std::size_t i = 0; i < mesh.nVertices(); ++i) {
      std::size_t xidx = 3 * i;
      std::size_t yidx = 3 * i + 1;
      std::size_t zidx = 3 * i + 2;
      // check x values
      ASSERT_EQ(geom.vertexPositions[i][0], pos1(i, 0));
      ASSERT_EQ(geom.vertexPositions[i][0], pos2(xidx));

      // check y values
      ASSERT_EQ(geom.vertexPositions[i][1], pos1(i, 1));
      ASSERT_EQ(geom.vertexPositions[i][1], pos2(yidx));

      // check z values
      ASSERT_EQ(geom.vertexPositions[i][2], pos1(i, 2));
      ASSERT_EQ(geom.vertexPositions[i][2], pos2(zidx));

      geom.vertexPositions[i][0] += 1;
      ASSERT_EQ(geom.vertexPositions[i][0], pos1(i, 0));
      ASSERT_EQ(geom.vertexPositions[i][0], pos2(xidx));
      ASSERT_EQ(geom.vertexPositions[i][1], pos1(i, 1));
      ASSERT_EQ(geom.vertexPositions[i][1], pos2(yidx));
      ASSERT_EQ(geom.vertexPositions[i][2], pos1(i, 2));
      ASSERT_EQ(geom.vertexPositions[i][2], pos2(zidx));

      pos1(i, 1) -= 3.14;
      ASSERT_EQ(geom.vertexPositions[i][0], pos1(i, 0));
      ASSERT_EQ(geom.vertexPositions[i][0], pos2(xidx));
      ASSERT_EQ(geom.vertexPositions[i][1], pos1(i, 1));
      ASSERT_EQ(geom.vertexPositions[i][1], pos2(yidx));
      ASSERT_EQ(geom.vertexPositions[i][2], pos1(i, 2));
      ASSERT_EQ(geom.vertexPositions[i][2], pos2(zidx));

      pos2(zidx) /= 4;
      ASSERT_EQ(geom.vertexPositions[i][0], pos1(i, 0));
      ASSERT_EQ(geom.vertexPositions[i][0], pos2(xidx));
      ASSERT_EQ(geom.vertexPositions[i][1], pos1(i, 1));
      ASSERT_EQ(geom.vertexPositions[i][1], pos2(yidx));
      ASSERT_EQ(geom.vertexPositions[i][2], pos1(i, 2));
      ASSERT_EQ(geom.vertexPositions[i][2], pos2(zidx));
    }
  }
}


TEST_F(MappingUtilityTest, MappingConstCorrectnessTest) {
  for (MeshAsset& a : allMeshes(true)) {
    SurfaceMesh& mesh = *a.mesh;
    const VertexData<Vector3>& data((a.geometry)->vertexPositions);

    const auto pos1 = EigenMap<double, 3>(data);
    ASSERT_EQ(data.size(), pos1.rows());
    ASSERT_EQ(3, pos1.cols());

    const auto pos2 = FlattenedEigenMap<double, 3>(data);
    ASSERT_EQ(3 * data.size(), pos2.rows());
    ASSERT_EQ(1, pos2.cols());

    for (std::size_t i = 0; i < mesh.nVertices(); ++i) {
      std::size_t xidx = 3 * i;
      std::size_t yidx = 3 * i + 1;
      std::size_t zidx = 3 * i + 2;
      // check x values
      ASSERT_EQ(data[i][0], pos1(i, 0));
      ASSERT_EQ(data[i][0], pos2(xidx));

      // check y values
      ASSERT_EQ(data[i][1], pos1(i, 1));
      ASSERT_EQ(data[i][1], pos2(yidx));

      // check z values
      ASSERT_EQ(data[i][2], pos1(i, 2));
      ASSERT_EQ(data[i][2], pos2(zidx));
    }
  }
}
