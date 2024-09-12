#pragma once

#include "geometrycentral/numerical/linear_solvers.h"
#include "geometrycentral/surface/embedded_geometry_interface.h"
#include "geometrycentral/surface/signed_heat_method.h"
#include "geometrycentral/surface/surface_mesh.h"

namespace geometrycentral {
namespace surface {

class PolygonMeshHeatSolver {

public:
  // === Constructor
  PolygonMeshHeatSolver(EmbeddedGeometryInterface& geom, double tCoef = 1.0);

  // === Methods

  // Solve for distance from a single vertex (or collection of points)
  VertexData<double> computeDistance(const Vertex& sourceVert);
  VertexData<double> computeDistance(const std::vector<Vertex>& sourceVerts);

  // Scalar Extension
  VertexData<double> extendScalars(const std::vector<std::tuple<Vertex, double>>& sources);

  // Compute parallel transport along shortest geodesics from sources at points
  VertexData<Vector2> transportTangentVector(const Vertex& sourceVert, const Vector2& sourceVector);
  VertexData<Vector2> transportTangentVectors(const std::vector<std::tuple<Vertex, Vector2>>& sources);

  // Solve for signed distance from a curve comprising a sequence of vertices.
  VertexData<double> computeSignedDistance(const std::vector<std::vector<Vertex>>& curves,
                                           const LevelSetConstraint& levelSetConstraint = LevelSetConstraint::ZeroSet);

  // === Options and parameters

  const double tCoef; // the time parameter used for heat flow, measured as time = tCoef * maxDiagonalLength
                      // default: 1.0

private:
  // === Members

  // Input mesh and geometry
  SurfaceMesh& mesh;
  EmbeddedGeometryInterface& geom;

  // Parameters
  double shortTime;

  // Solvers
  void ensureHaveScalarHeatSolver();
  void ensureHaveVectorHeatSolver();
  void ensureHavePoissonSolver();
  std::unique_ptr<PositiveDefiniteSolver<std::complex<double>>> vectorHeatSolver;
  std::unique_ptr<PositiveDefiniteSolver<double>> scalarHeatSolver, poissonSolver;
  SparseMatrix<double> massMat, laplaceMat;

  // Helpers
  void buildSignedCurveSource(const std::vector<Vertex>& curve, Vector<std::complex<double>>& X0) const;
  double computeAverageValue(const std::vector<std::vector<Vertex>>& curves, const Vector<double>& u);
};

} // namespace surface
} // namespace geometrycentral