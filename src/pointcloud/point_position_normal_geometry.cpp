#include "geometrycentral/pointcloud/point_position_normal_geometry.h"

namespace geometrycentral {
namespace pointcloud {


PointPositionNormalGeometry::PointPositionNormalGeometry(PointCloud& cloud, const PointData<Vector3>& positions_,
                                                         const PointData<Vector3>& normals_)
    : PointPositionGeometry(cloud, positions_) {
  normals = normals_;
  normalsQ.clearable = false;
}

void PointPositionNormalGeometry::computeNormals() {
  // do nothing; already populated
}

} // namespace pointcloud
} // namespace geometrycentral
