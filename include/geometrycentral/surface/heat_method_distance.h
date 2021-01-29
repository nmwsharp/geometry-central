#pragma once

#include "geometrycentral/numerical/linear_solvers.h"
#include "geometrycentral/surface/edge_length_geometry.h"
#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/utilities/vector2.h"
#include "geometrycentral/utilities/vector3.h"

namespace geometrycentral {


namespace pointcloud {
class PointCloudHeatSolver; // forward declare to friend below
}

namespace surface {

// One-off function to compute distance from a vertex
VertexData<double> heatMethodDistance(IntrinsicGeometryInterface& geom, Vertex v);


// Stateful class. Allows efficient repeated solves
class HeatMethodDistanceSolver {

public:
  // === Constructor
  HeatMethodDistanceSolver(IntrinsicGeometryInterface& geom, double tCoef = 1.0, bool useRobustLaplacian = false);

  // === Methods

  // Solve for distance from a single vertex
  VertexData<double> computeDistance(const Vertex& sourceVert);

  // Solve for distance from a collection of vertices
  VertexData<double> computeDistance(const std::vector<Vertex>& sourceVerts);

  // Solve for distance from a single surface point
  VertexData<double> computeDistance(const SurfacePoint& sourcePoint);

  // Solve for distance from a collection of surface points
  VertexData<double> computeDistance(const std::vector<SurfacePoint>& sourcePoints);

  // Solve for distance from a custom right hand side
  // (returns WITHOUT performing constant shift to 0)
  Vector<double> computeDistanceRHS(const Vector<double>& rhs);

  // === Options and parameters

  const double tCoef; // the time parameter used for heat flow, measured as time = tCoef * mean_edge_length^2
                      // default: 1.0
  const bool useRobustLaplacian;

private:
  // === Members

  // Input mesh and geometry
  SurfaceMesh& mesh;
  IntrinsicGeometryInterface& geom;

  // Tufted cover mesh & geometry (these will only be popualted if useRobustLaplacian = true)
  std::unique_ptr<SurfaceMesh> tuftedMesh;
  std::unique_ptr<EdgeLengthGeometry> tuftedIntrinsicGeom;

  // Parameters
  double shortTime; // the actual time used for heat flow computed from tCoef

  // Solvers
  std::unique_ptr<PositiveDefiniteSolver<double>> heatSolver;
  std::unique_ptr<PositiveDefiniteSolver<double>> poissonSolver;

  // Helpers

  // Return either the input mesh/geometry, or the tufted mesh/geometry, based on whether useRobustLaplacian=true
  SurfaceMesh& getMesh();
  IntrinsicGeometryInterface& getGeom();

  // The point cloud version bootstraps off of this code
  friend class pointcloud::PointCloudHeatSolver;
};


} // namespace surface
} // namespace geometrycentral
