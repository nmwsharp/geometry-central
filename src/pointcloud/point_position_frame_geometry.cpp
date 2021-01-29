#include "geometrycentral/pointcloud/point_position_frame_geometry.h"

namespace geometrycentral {
namespace pointcloud {


PointPositionFrameGeometry::PointPositionFrameGeometry(PointCloud& cloud, const PointData<Vector3>& positions_,
                                                       const PointData<std::array<Vector3, 3>>& frames)
    : PointPositionGeometry(cloud, positions_) {

  normals = PointData<Vector3>(cloud);
  tangentBasis = PointData<std::array<Vector3, 2>>(cloud);

  for (Point p : cloud.points()) {
    tangentBasis[p][0] = frames[p][0];
    tangentBasis[p][1] = frames[p][1];
    normals[p] = frames[p][2];
  }

  normalsQ.clearable = false;
  tangentBasisQ.clearable = false;
}

void PointPositionFrameGeometry::computeNormals() {
  // do nothing; already populated
}
void PointPositionFrameGeometry::computeTangentBasis() {
  // do nothing; already populated
}

} // namespace pointcloud
} // namespace geometrycentral
