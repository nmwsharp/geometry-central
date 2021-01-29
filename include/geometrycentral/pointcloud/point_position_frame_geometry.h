#pragma once

#include "geometrycentral/pointcloud/point_position_geometry.h"

namespace geometrycentral {
namespace pointcloud {

// Extends the point position geometry by also specifying known normals & tangent frames

class PointPositionFrameGeometry : public PointPositionGeometry {

public:
  // Frame input should be {basisX, basisY, normal}, forming an orthonormal right-handed coordinate system.
  PointPositionFrameGeometry(PointCloud& cloud, const PointData<Vector3>& positions,
                             const PointData<std::array<Vector3, 3>>& frames);

  // Normals & tangents are stored in the `normals` and `tangentBasis` fields already present in the parent class

protected:
  // Override the normal/tangent computation to do nothing; we just use the existing normals.
  void computeNormals() override;
  void computeTangentBasis() override;
};

} // namespace pointcloud
} // namespace geometrycentral
