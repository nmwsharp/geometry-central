#pragma once

#include "geometrycentral/numerical/linear_algebra_types.h"
#include "geometrycentral/numerical/linear_solvers.h"
#include "geometrycentral/surface/intrinsic_triangulation.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/simple_polygon_mesh.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/utilities/elementary_geometry.h"

#include <deque>

namespace geometrycentral {
namespace surface {

enum class TransferMethod { Pointwise = 0, L2 };

class AttributeTransfer {

public:
  // Constructor
  AttributeTransfer(CommonSubdivision& cs, VertexPositionGeometry& geomA);
  AttributeTransfer(IntrinsicTriangulation& intTri);

  // Members
  CommonSubdivision& cs;

  SparseMatrix<double> P_A;           // maps scalars at vertices to CS
  SparseMatrix<double> P_B;           // maps scalars at vertices to CS
  SparseMatrix<double> M_CS_Galerkin; // galerkin mass matrix on common subdivision

  std::unique_ptr<SquareSolver<double>> AtoB_L2_Solver;
  std::unique_ptr<SquareSolver<double>> BtoA_L2_Solver;


  // Methods

  // High-level
  VertexData<double> transferAtoB(const VertexData<double>& valuesOnA, TransferMethod method);
  VertexData<double> transferBtoA(const VertexData<double>& valuesOnB, TransferMethod method);

  // Low-level
  VertexData<double> transferAtoB_Pointwise(const VertexData<double>& valuesOnA);
  VertexData<double> transferAtoB_L2(const VertexData<double>& valuesOnA);

  VertexData<double> transferBtoA_Pointwise(const VertexData<double>& valuesOnB);
  VertexData<double> transferBtoA_L2(const VertexData<double>& valuesOnB);

  // Prepare data
  std::pair<SparseMatrix<double>, SparseMatrix<double>> constructAtoBMatrices() const;
  std::pair<SparseMatrix<double>, SparseMatrix<double>> constructBtoAMatrices() const;
};


// One-off functions
VertexData<double> transferAtoB(CommonSubdivision& cs, VertexPositionGeometry& geomA,
                                const VertexData<double>& valuesOnA, TransferMethod method);
VertexData<double> transferAtoB(IntrinsicTriangulation& intTri, const VertexData<double>& valuesOnA,
                                TransferMethod method);

VertexData<double> transferBtoA(CommonSubdivision& cs, VertexPositionGeometry& geomA,
                                const VertexData<double>& valuesOnB, TransferMethod method);
VertexData<double> transferBtoA(IntrinsicTriangulation& intTri, const VertexData<double>& valuesOnB,
                                TransferMethod method);


} // namespace surface
} // namespace geometrycentral
