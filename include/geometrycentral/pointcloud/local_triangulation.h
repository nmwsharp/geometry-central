#pragma once

#include "geometrycentral/pointcloud/geometry3D.h"
#include "geometrycentral/pointcloud/point_cloud.h"

namespace geometrycentral {
namespace pointcloud {

// In the local neighborhood of each triangle, construct the triangles in the local planar Delaunay triangulation which
// touch the center point.
// The optional perturbation heuristic can help to ensure that some triangles are generated in utterly degenerate cases,
// such as many exactly coincident points.
PointData<std::vector<std::array<Point, 3>>> buildLocalTriangulations(PointCloud& cloud, Geometry3D& geom,
                                                                      bool withDegeneracyHeuristic = true);

// Convert a local neighbor indexed list to global indices
PointData<std::vector<std::array<size_t, 3>>>
handleToInds(PointCloud& cloud, const PointData<std::vector<std::array<Point, 3>>>& handleResult);

std::vector<std::vector<size_t>> handleToFlatInds(PointCloud& cloud,
                                                  const PointData<std::vector<std::array<Point, 3>>>& handleResult);


} // namespace pointcloud
} // namespace geometrycentral

