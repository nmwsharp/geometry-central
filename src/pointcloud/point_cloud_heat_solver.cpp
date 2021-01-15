#include "geometrycentral/pointcloud/point_cloud_heat_solver.h"


namespace geometrycentral {
namespace pointcloud {

PointCloudHeatSolver::PointCloudHeatSolver(PointCloud& cloud_, Geometry3D& geom_, double tCoef_)
    : tCoef(tCoef_), cloud(cloud_), geom(geom_) {

  geom.requireNeighbors();
  geom.requireTuftedTriangulation();
  // geom.tuftedGeom->requireEdgeLengths();
  // geom.tuftedGeom->requireCotanLaplacian();
  // geom.requireLaplacian();
  // geom.requireConnectionLaplacian();

  heatDistanceWorker.reset(new surface::HeatMethodDistanceSolver(*geom.tuftedGeom, tCoef));
  shortTime = heatDistanceWorker->shortTime;
}

// === Heat method for distance
// For distance, we basically just run the mesh version on the tufted triangulation of the point cloud.
PointData<double> PointCloudHeatSolver::computeDistance(const Point& sourcePoint) {
  std::vector<Point> v{sourcePoint};
  return computeDistance(v);
}
PointData<double> PointCloudHeatSolver::computeDistance(const std::vector<Point>& sourcePoints) {
  std::vector<surface::Vertex> sourceVerts;
  for(Point p  : sourcePoints) {
    sourceVerts.push_back(geom.tuftedMesh->vertex(p.getIndex()));
  }
  return PointData<double>(cloud, heatDistanceWorker->computeDistance(sourceVerts).raw());
}

} // namespace pointcloud
} // namespace geometrycentral
