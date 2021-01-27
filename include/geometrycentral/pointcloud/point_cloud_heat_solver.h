#pragma once

#include "geometrycentral/numerical/linear_solvers.h"
#include "geometrycentral/pointcloud/point_cloud.h"
#include "geometrycentral/pointcloud/point_position_geometry.h"
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
  PointCloudHeatSolver(PointCloud& cloud, PointPositionGeometry& geom, double tCoef = 1.0);

  // === Methods

  // Solve for distance from a single vertex (or collection of points)
  PointData<double> computeDistance(const Point& sourcePoint);
  PointData<double> computeDistance(const std::vector<Point>& sourcePoints);

  // === Scalar Extension
  PointData<double> extendScalars(const std::vector<std::tuple<Point, double>>& sources);

  // Compute parallel transport along shortest geodesics from sources at points
  PointData<Vector2> transportTangentVector(const Point& sourcePoint, const Vector2& sourceVector);
  PointData<Vector2> transportTangentVectors(const std::vector<std::tuple<Point, Vector2>>& sources);

  // Compute the logarithmic map from a source point
  PointData<Vector2> computeLogMap(const Point& sourcePoint);

  // === Options and parameters

  const double tCoef; // the time parameter used for heat flow, measured as time = tCoef * mean_edge_length^2
                      // default: 1.0

private:
  // === Members

  // Input cloud and geometry
  PointCloud& cloud;
  PointPositionGeometry& geom;

  // Parameters
  double shortTime; // the actual time used for heat flow computed from tCoef

  // Populate solvers lazily as needed
  void ensureHaveHeatDistanceWorker();
  void ensureHaveVectorHeatSolver();

  // Lazy strategy: bootstrap off the mesh version of the solver for distance solves on the tufted mesh
  std::unique_ptr<surface::HeatMethodDistanceSolver> heatDistanceWorker;

  // Solvers
  std::unique_ptr<PositiveDefiniteSolver<double>> vectorHeatSolver;
};

} // namespace pointcloud
} // namespace geometrycentral
