#pragma once

#include "geometrycentral/pointcloud/point_cloud.h"
#include "geometrycentral/utilities/vector3.h"

namespace geometrycentral {
namespace pointcloud {

// Represents a set of neighborhoods for a point cloud. Immutable after construction.
//
// TODO upgrade this to be much cleverer and use a flat representation w/ custom iterators and containers.
//
class Neighborhoods {
public:
  // == Constructors
  Neighborhoods(PointCloud& cloud, const PointData<Vector3>& positions, unsigned int nNeighbors);
  // TODO constructor for fixed-ball neighborhood

  // == Members
  PointCloud& cloud;

  PointData<std::vector<Point>> neighbors;

  // == Functions

protected:
  // TODO
  // PointData<size_t> listIndEnd;
  // std::vector<size_t> flatNeighborList;
};


} // namespace pointcloud
} // namespace geometrycentral
