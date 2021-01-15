#pragma once

#include "geometrycentral/numerical/linear_solvers.h"
#include "geometrycentral/pointcloud/geometry3D.h"
#include "geometrycentral/pointcloud/point_cloud.h"
#include "geometrycentral/surface/edge_length_geometry.h"
#include "geometrycentral/surface/heat_method_distance.h"
#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/utilities/vector2.h"
#include "geometrycentral/utilities/vector3.h"

namespace geometrycentral {
namespace pointcloud {

class PointCloudHeatSolver {

public:
  // === Constructor
  PointCloudHeatSolver(PointCloud& cloud, Geometry3D& geom, double tCoef = 1.0);

  // === Methods

  // Solve for distance from a single vertex (or collection of vertices)
  PointData<double> computeDistance(const Point& sourcePoint);
  PointData<double> computeDistance(const std::vector<Point>& sourcePoints);

  // Compute parallel transport along shortest geodesics from sources at vertices
  PointData<Vector2> computeParallelTransport(const Point& sourcePoint, const Vector2& sourceVector);
  PointData<Vector2> computeParallelTransport(const std::vector<Point>& sourcePoints,
                                              const std::vector<Vector2>& sourceVectors);

  // Compute the logarithmic map from a source point
  PointData<Vector2> computeLogMap(const Point& sourcePoint);

  // === Options and parameters

  const double tCoef; // the time parameter used for heat flow, measured as time = tCoef * mean_edge_length^2
                      // default: 1.0

private:
  // === Members

  // Input cloud and geometry
  PointCloud& cloud;
  Geometry3D& geom;

  // Parameters
  double shortTime; // the actual time used for heat flow computed from tCoef

  // Lazy strategy: bootstrap off the mesh version of the solver for distance solves on the tufted mesh
  std::unique_ptr<surface::HeatMethodDistanceSolver> heatDistanceWorker;

  // Solvers
  std::unique_ptr<PositiveDefiniteSolver<std::complex<double>>> vectorHeatSolver;
};

} // namespace pointcloud
} // namespace geometrycentral
