#pragma once

#include "geometrycentral/pointcloud/point_position_geometry.h"

namespace geometrycentral {
namespace pointcloud {

// Extends the point position geometry by also specifying known normals

class PointPositionNormalGeometry : public PointPositionGeometry {

public:
  PointPositionNormalGeometry(PointCloud& cloud, const PointData<Vector3>& positions, const PointData<Vector3>& normals);

  // Normals are stored in the `normals` field already present in the parent class


protected:
  // Override the normal computation to do nothing; we just use the existing normals.
  virtual void computeNormals() override;
};

} // namespace pointcloud
} // namespace geometrycentral
