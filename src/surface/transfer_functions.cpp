#include "geometrycentral/surface/transfer_functions.h"


AttributeTransfer::AttributeTransfer(CommonSubdivision& cs_, VertexPositionGeometry& geomA) : cs(cs_) {
  // Make sure the common subdivision has full mesh connectivity
  // (this is lazy and expensive; we could actually get away without it)
  if (!cs.mesh) cs.constructMesh();

  // === Prepare data
  P_A = cs.interpolationMatrixA();
  P_B = cs.interpolationMatrixB();
  M_CS_Galerkin = cs.vertexGalerkinMassMatrixFromPositionsA(geomA.vertexPositions);
}

AttributeTransfer::AttributeTransfer(CommonSubdivision& cs_, IntrinsicTriangulation& intTri) : cs(cs_) {
  // Make sure the common subdivision has full mesh connectivity
  // (this is lazy and expensive; we could actually get away without it)
  if (!cs.mesh) cs.constructMesh();

  // === Prepare data
  P_A = cs.interpolationMatrixA();
  P_B = cs.interpolationMatrixB();
  M_CS_Galerkin = cs.vertexGalerkinMassMatrixFromLengthsB(intTri.edgeLengths);
}


VertexData<double> AttributeTransfer::transferAtoB(const VertexData<double>& valuesOnA, TransferMethod method) {

  switch (method) {
  case TransferMethod::Pointwise: {
    return transferAtoB_Pointwise(valuesOnA);
  }
  case TransferMethod::L2: {
    return transferAtoB_L2(valuesOnA);
  }
  }

  return VertexData<double>(); // unreachable
}
VertexData<double> AttributeTransfer::transferBtoA(const VertexData<double>& valuesOnB, TransferMethod method) {

  switch (method) {
  case TransferMethod::Pointwise: {
    return transferBtoA_Pointwise(valuesOnB);
  }
  case TransferMethod::L2: {
    return transferBtoA_L2(valuesOnB);
  }
  }

  return VertexData<double>(); // unreachable
}


VertexData<double> AttributeTransfer::transferAtoB_Pointwise(const VertexData<double>& valuesOnA) {

  // TODO could implement this more simply using the intrinsic triangulation
  // object, but this way uses only data in common subdivision.

  VertexData<double> result(cs.meshB);
  for (Vertex v : cs.mesh->vertices()) {
    CommonSubdivisionPoint& p = *cs.sourcePoints[v];
    if (p.posB.type == SurfacePointType::Vertex) {
      result[p.posB.vertex] = p.posA.interpolate(valuesOnA);
    }
  }
  return result;
}

VertexData<double> AttributeTransfer::transferAtoB_L2(const VertexData<double>& valuesOnA) {
  if (!AtoB_L2_Solver) {
    SparseMatrix<double> mat = P_B.transpose() * M_CS_Galerkin * P_B;
    AtoB_L2_Solver.reset(new Solver<double>(mat));
  }
  Vector<double> vec = P_B.transpose() * M_CS_Galerkin * P_A * valuesOnA.toVector();
  Vector<double> result = AtoB_L2_Solver->solve(vec);
  return VertexData<double>(cs.meshB, result);
}

VertexData<double> AttributeTransfer::transferBtoA_Pointwise(const VertexData<double>& valuesOnB) {

  // TODO could implement this more simply using the intrinsic triangulation
  // object, but this way uses only data in common subdivision.

  VertexData<double> result(cs.meshA);
  for (Vertex v : cs.mesh->vertices()) {
    CommonSubdivisionPoint& p = *cs.sourcePoints[v];
    if (p.posA.type == SurfacePointType::Vertex) {
      result[p.posA.vertex] = p.posB.interpolate(valuesOnB);
    }
  }
  return result;
}

VertexData<double> AttributeTransfer::transferBtoA_L2(const VertexData<double>& valuesOnB) {
  if (!BtoA_L2_Solver) {
    SparseMatrix<double> mat = P_A.transpose() * M_CS_Galerkin * P_A;
    BtoA_L2_Solver.reset(new Solver<double>(mat));
  }
  Vector<double> vec = P_A.transpose() * M_CS_Galerkin * P_B * valuesOnB.toVector();

  Vector<double> result = BtoA_L2_Solver->solve(vec);
  return VertexData<double>(cs.meshA, result);
}

VertexData transferAtoB(CommonSubdivision& cs, VertexPositionGeometry& geomA, const VertexData<double>& valuesOnA,
                        TransferMethod method) {
  AttributeTransfer transfer(cs, geomA);
  return transfer.transferAtoB(valuesOnA, method);
}

VertexData transferAtoB(IntrinsicTriangulation& intTri, const VertexData<double>& valuesOnA, TransferMethod method) {
  std::unique_ptr<CommonSubdivision> cs = intTri.extractCommonSubdivision();
  AttributeTransfer transfer(*cs, intTri);
  return transfer.transferAtoB(valuesOnA, method);
}

VertexData transferBtoA(CommonSubdivision& cs, VertexPositionGeometry& geomA, const VertexData<double>& valuesOnB,
                        TransferMethod method) {
  AttributeTransfer transfer(cs, geomA);
  return transfer.transferBtoA(valuesOnB, method);
}

VertexData transferBtoA(IntrinsicTriangulation& geomA, const VertexData<double>& valuesOnB, TransferMethod method) {
  AttributeTransfer transfer(*cs, intTri);
  return transfer.transferBtoA(valuesOnB, method);
}
