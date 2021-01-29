#pragma once

#include "geometrycentral/pointcloud/point_cloud.h"
#include "geometrycentral/pointcloud/point_position_geometry.h"

namespace geometrycentral {
namespace pointcloud {

// === Readers ===

// Read from a file by name. Type can be optionally inferred from filename.
std::tuple<std::unique_ptr<PointCloud>, std::unique_ptr<PointPositionGeometry>> readPointCloud(std::string filename,
                                                                                               std::string type = "");

// Same as above, but from an istream. Must specify type.
std::tuple<std::unique_ptr<PointCloud>, std::unique_ptr<PointPositionGeometry>> readPointCloud(std::istream& in,
                                                                                               std::string type);


// === Writers ===

// Write to file by name. Type can be optionally inferred from filename.
void writePointCloud(PointCloud& cloud, PointPositionGeometry& geometry, std::string filename, std::string type = "");

// Same as above, to to an ostream. Must specify type.
void writePointCloud(PointCloud& cloud, PointPositionGeometry& geometry, std::ostream& out, std::string type);

} // namespace pointcloud
} // namespace geometrycentral
