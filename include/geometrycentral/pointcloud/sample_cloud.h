#pragma once

#include "geometrycentral/pointcloud/point_cloud.h"
#include "geometrycentral/pointcloud/point_position_geometry.h"
#include "geometrycentral/surface/embedded_geometry_interface.h"
#include "geometrycentral/surface/surface_point.h"

#include <tuple>

namespace geometrycentral {
namespace pointcloud {

std::tuple<std::unique_ptr<PointCloud>, PointData<Vector3>, PointData<surface::SurfacePoint>>
uniformlySamplePointsOnSurface(surface::SurfaceMesh& mesh, surface::EmbeddedGeometryInterface& geom, size_t nPts);


}
} // namespace geometrycentral
