#pragma once

#include "geometrycentral/numerical/linear_algebra_utilities.h"
#include "geometrycentral/numerical/linear_solvers.h"
#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"

namespace geometrycentral {
namespace surface {

VertexData<Vector2> parameterizeBFF(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geo);
VertexData<Vector2> parameterizeBFFfromScaleFactors(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geo,
                                                    const VertexData<double>& boundaryScaleFactors);
VertexData<Vector2> parameterizeBFFfromExteriorAngles(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geo,
                                                      const VertexData<double>& exteriorAngles);

class BFF {
public:
  BFF(ManifoldSurfaceMesh& mesh_, IntrinsicGeometryInterface& geo_);

  VertexData<Vector2> flatten();
  VertexData<Vector2> flattenFromScaleFactors(const VertexData<double>& uBdy);
  VertexData<Vector2> flattenFromExteriorAngles(const VertexData<double>& kBdy);
  VertexData<Vector2> flattenFromBoth(const Vector<double>& uBdy, const Vector<double>& kBdy);

  Vector<double> dirichletToNeumann(const Vector<double>& uBdy);
  Vector<double> neumannToDirichlet(const Vector<double>& kBdy);

  std::array<Vector<double>, 2> computeBoundaryPositions(const Vector<double>& uBdy, const Vector<double>& kBdy);

protected:
  ManifoldSurfaceMesh& mesh;
  IntrinsicGeometryInterface& geo;
  VertexData<size_t> vIdx;
  VertexData<int> iIdx, bIdx;


  SparseMatrix<double> L, Lii, Lib, Lbb;
  std::unique_ptr<PositiveDefiniteSolver<double>> Liisolver;
  std::unique_ptr<PositiveDefiniteSolver<double>> Lsolver;
  Vector<double> Omegai, Omegab;

  Vector<bool> isInterior;
  size_t nInterior, nBoundary, nVertices;

  BlockDecompositionResult<double> Ldecomp;

  void ensureHaveLSolver();
};

} // namespace surface
} // namespace geometrycentral
