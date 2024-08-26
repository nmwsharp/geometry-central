#pragma once

#include "geometrycentral/numerical/linear_solvers.h"
#include "geometrycentral/surface/barycentric_vector.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/surface_point.h"

namespace geometrycentral {

#ifndef SHM_H
#define SHM_H
enum class LevelSetConstraint { None = 0, ZeroSet, Multiple };
#endif

namespace surface {

struct Curve {
  std::vector<SurfacePoint> nodes;
  bool isSigned = true;
};

struct SignedHeatOptions {
  bool preserveSourceNormals = false;
  LevelSetConstraint levelSetConstraint = LevelSetConstraint::ZeroSet;
  double softLevelSetWeight = -1.;
};

class SignedHeatMethodSolver {

public:
  // === Constructor
  SignedHeatMethodSolver(IntrinsicGeometryInterface& geom, double tCoef = 1.0);

  VertexData<double> computeDistance(const std::vector<Curve>& curves,
                                     const std::vector<SurfacePoint>& points = std::vector<SurfacePoint>(),
                                     const SignedHeatOptions& options = SignedHeatOptions());

  VertexData<double> computeDistance(const std::vector<Curve>& curves,
                                     const SignedHeatOptions& options = SignedHeatOptions());

  VertexData<double> computeDistance(const std::vector<SurfacePoint>& points);

  // === Options and parameters
  void setDiffusionTimeCoefficient(double tCoef = 1.0);

private:
  // === Members
  SurfaceMesh& mesh;
  IntrinsicGeometryInterface& geom;

  // Parameters
  double shortTime;
  double meanNodeDistance;
  bool timeUpdated = false;

  // Solvers
  std::unique_ptr<LinearSolver<std::complex<double>>> vectorHeatSolver;
  std::unique_ptr<PositiveDefiniteSolver<double>> poissonSolver;
  SparseMatrix<double> massMat, doubleVectorOp;

  // Helpers
  Vector<std::complex<double>> integrateVectorHeatFlow(const std::vector<Curve>& curves,
                                                       const std::vector<SurfacePoint>& points,
                                                       const SignedHeatOptions& options) const;
  VertexData<double> integrateVectorField(const Vector<std::complex<double>>& Xt, const std::vector<Curve>& curves,
                                          const std::vector<SurfacePoint>& points, const SignedHeatOptions& options);

  void ensureHaveVectorHeatSolver();
  void ensureHavePoissonSolver();

  SparseMatrix<double> crouzeixRaviartDoubleConnectionLaplacian() const;
  SparseMatrix<double> crouzeixRaviartDoubleMassMatrix() const;

  void buildSignedCurveSource(const Curve& curve, Vector<std::complex<double>>& X0) const;
  void buildUnsignedVertexSource(const Vertex& v, Vector<std::complex<double>>& X0, double weight = 1.) const;
  void buildUnsignedPointSource(const SurfacePoint& point, Vector<std::complex<double>>& X0) const;
  double lengthOfSegment(const SurfacePoint& pA, const SurfacePoint& pB) const;
  SurfacePoint midSegmentSurfacePoint(const SurfacePoint& pA, const SurfacePoint& pB) const;
  std::complex<double> projectedNormal(const SurfacePoint& pA, const SurfacePoint& pB, const Edge& e) const;
  BarycentricVector barycentricVectorInFace(const Halfedge& he, const Face& f) const;
  FaceData<BarycentricVector> sampleAtFaceBarycenters(const Vector<std::complex<double>>& Xt);
  Vector<double> integrateWithZeroSetConstraint(const Vector<double>& rhs, const std::vector<Curve>& curves,
                                                const SignedHeatOptions& options);
  Vector<double> integrateWithLevelSetConstraints(const Vector<double>& rhs, const std::vector<Curve>& curves,
                                                  const SignedHeatOptions& options);
  double computeAverageValueOnSource(const Vector<double>& phi, const std::vector<Curve>& curves) const;
};

} // namespace surface
} // namespace geometrycentral