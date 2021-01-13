#pragma once

#include "geometrycentral/surface/embedded_geometry_interface.h"
#include "geometrycentral/pointcloud/geometry3D.h"
#include "geometrycentral/pointcloud/point_cloud.h"
#include "geometrycentral/surface/surface_point.h"

#include <tuple>

namespace geometrycentral {
namespace pointcloud {

std::tuple<std::unique_ptr<PointCloud>, std::unique_ptr<Geometry3D>, PointData<surface::SurfacePoint>>
uniformlySamplePointsOnSurface(surface::SurfaceMesh& mesh, surface::EmbeddedGeometryInterface& geom, size_t nPts);


}
} // namespace geometrycentral
