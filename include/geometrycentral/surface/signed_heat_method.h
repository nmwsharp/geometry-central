#pragma once

#include "geometrycentral/numerical/linear_solvers.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/surface_point.h"

namespace geometrycentral {
namespace surface {

#ifndef SHM_H
#define SHM_H
enum class LevelSetConstraint { None = 0, ZeroSet, Multiple };
#endif

struct Curve {
  std::vector<SurfacePoint> nodes;
  bool isSigned = true;
};

class SignedHeatMethodSolver {

public:
  // === Constructor
  SignedHeatMethodSolver(IntrinsicGeometryInterface& geom, double tCoef = 1.0);

  VertexData<double> computeDistance(const std::vector<Curve>& curves,
                                     const std::vector<SurfacePoint>& points = std::vector<SurfacePoint>(),
                                     int levelSetConstraint = LevelSetConstraint::ZeroSet,
                                     bool useSoftLevelSetConstraint = false, double softLevelSetWeight = 0.,
                                     bool solvePiecewiseDistance = false);

  VertexData<double> computeDistance(const std::vector<Curve>& curves,
                                     int levelSetConstraint = LevelSetConstraint::ZeroSet,
                                     bool useSoftLevelSetConstraint = false, double softLevelSetWeight = 0.,
                                     bool solvePiecewiseDistance = false);

  VertexData<double> computeDistance(const std::vector<SurfacePoint>& points);

  // === Options and parameters
  const double tCoef; // the time parameter used for heat flow, measured as time = tCoef * mean_edge_length^2
                      // default: 1.0

private:
  // === Members
  SurfaceMesh& mesh;
  IntrinsicGeometryInterface& geom;

  // Parameters
  double shortTime;

  // Solvers
  std::unique_ptr<LinearSolver<std::complex<double>>> vectorHeatSolver;
  std::unique_ptr<PositiveDefiniteSolver<double>> poissonSolver;
  SparseMatrix<double> massMat;

  // Helpers
  void ensureHaveVectorHeatSolver();
  void ensureHavePoissonSolver();
}