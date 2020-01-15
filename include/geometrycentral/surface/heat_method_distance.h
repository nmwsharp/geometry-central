#pragma once

#include "geometrycentral/utilities/vector3.h"

#include "geometrycentral/utilities/vector2.h"

#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/surface/intrinsic_geometry_interface.h"

#include "geometrycentral/numerical/linear_solvers.h"

namespace geometrycentral {
namespace surface {

// One-off function to compute distance from a vertex
VertexData<double> heatMethodDistance(IntrinsicGeometryInterface& geom, Vertex v);


// Stateful class. Allows efficient repeated solves

enum class ComputeTriangulation { Original = 0, IntrinsicDelaunay, IntrinsicDelaunayRefine };

class HeatMethodDistanceSolver {

public:
  // === Constructor
  HeatMethodDistanceSolver(IntrinsicGeometryInterface& geom, double tCoef = 1.0);


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


  // what triangulation to perform the computation on
  // TODO not supported yet
  const ComputeTriangulation computeTri = ComputeTriangulation::Original;


private:
  
  // === Members

  // Basics
  HalfedgeMesh& mesh;
  IntrinsicGeometryInterface& geom;

  // Parameters
  double shortTime;   // the actual time used for heat flow computed from tCoef

  // Solvers
  std::unique_ptr<PositiveDefiniteSolver<double>> heatSolver;
  std::unique_ptr<PositiveDefiniteSolver<double>> poissonSolver;
  
};


} // namespace surface
} // namespace geometrycentral
