#include "geometrycentral/surface/boundary_first_flattening.h"

namespace geometrycentral {
namespace surface {

VertexData<Vector2> parameterizeBFF(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geo) {
  BFF bff(mesh, geo);
  return bff.flatten();
}

VertexData<Vector2> parameterizeBFFfromScaleFactors(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geo,
                                                    const VertexData<double>& boundaryScaleFactors) {
  BFF bff(mesh, geo);
  return bff.flattenFromScaleFactors(boundaryScaleFactors);
}

VertexData<Vector2> parameterizeBFFfromExteriorAngles(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geo,
                                                      const VertexData<double>& exteriorAngles) {
  BFF bff(mesh, geo);
  return bff.flattenFromExteriorAngles(exteriorAngles);
}

BFF::BFF(ManifoldSurfaceMesh& mesh_, IntrinsicGeometryInterface& geo_) : mesh(mesh_), geo(geo_) {

  GC_SAFETY_ASSERT(mesh.eulerCharacteristic() == 2 && mesh.nBoundaryLoops() == 1,
                   "Input to BFF must be a topological disk");

  vIdx = mesh.getVertexIndices();
  iIdx = VertexData<int>(mesh, -1);
  bIdx = VertexData<int>(mesh, -1);

  geo.requireCotanLaplacian();

  L = geo.cotanLaplacian;

  isInterior = Vector<bool>(mesh.nVertices());
  size_t iB = 0, iI = 0;
  for (Vertex v : mesh.vertices()) {
    if (v.isBoundary()) {
      bIdx[v] = iB++;
      isInterior(vIdx[v]) = false;
    } else {
      iIdx[v] = iI++;
      isInterior(vIdx[v]) = true;
    }
  }
  nVertices = mesh.nVertices();
  nInterior = mesh.nInteriorVertices();
  nBoundary = nVertices - nInterior;

  geo.requireCotanLaplacian();
  SparseMatrix<double> L = geo.cotanLaplacian;
  shiftDiagonal(L, 1e-12);

  Ldecomp = blockDecomposeSquare(L, isInterior);

  Lii = Ldecomp.AA;
  Lib = Ldecomp.AB;
  Lbb = Ldecomp.BB;

  // TODO: extract this factorization from a full factorization of L
  Liisolver.reset(new PositiveDefiniteSolver<double>(Lii));

  geo.requireVertexAngleSums();
  Omegai = Vector<double>(nInterior);
  Omegab = Vector<double>(nBoundary);
  for (Vertex v : mesh.vertices()) {
    if (v.isBoundary()) {
      Omegab(bIdx[v]) = M_PI - geo.vertexAngleSums[v];
    } else {
      Omegai(iIdx[v]) = 2 * M_PI - geo.vertexAngleSums[v];
    }
  }
}

VertexData<Vector2> BFF::flatten() {
  VertexData<double> u(mesh, 0);
  return flattenFromScaleFactors(u);
}

VertexData<Vector2> BFF::flattenFromScaleFactors(const VertexData<double>& uData) {
  // Extract boundary values
  Vector<double> uBdy, ignore;
  decomposeVector(Ldecomp, uData.toVector(), ignore, uBdy);

  // Compute complementary data
  Vector<double> kBdy = dirichletToNeumann(uBdy);

  // Flatten
  return flattenFromBoth(uBdy, kBdy);
}

VertexData<Vector2> BFF::flattenFromExteriorAngles(const VertexData<double>& kData) {
  // Extract boundary values
  Vector<double> kBdy, ignore;
  decomposeVector(Ldecomp, kData.toVector(), ignore, kBdy);

  GC_SAFETY_ASSERT(abs(kBdy.sum() - 2 * M_PI) < 1e-3,
                   "BFF error: target exterior angles must sum to 2 pi, but the input sums to " +
                       std::to_string(kBdy.sum()));

  // Compute complementary data
  Vector<double> uBdy = neumannToDirichlet(kBdy);

  // Flatten
  return flattenFromBoth(uBdy, kBdy);
}

VertexData<Vector2> BFF::flattenFromBoth(const Vector<double>& uBdy, const Vector<double>& kBdy) {

  Vector<double> boundaryX, boundaryY;
  std::tie(boundaryX, boundaryY) = tuple_cat(computeBoundaryPositions(uBdy, kBdy));

  Vector<double> interiorX = Liisolver->solve(-Lib * boundaryX);
  Vector<double> interiorY = Liisolver->solve(-Lib * boundaryY);

  VertexData<Vector2> parm(mesh);
  for (Vertex v : mesh.vertices()) {
    if (v.isBoundary()) {
      size_t iV = bIdx[v];
      parm[v] = Vector2{boundaryX(iV), boundaryY(iV)};
    } else {
      size_t iV = iIdx[v];
      parm[v] = Vector2{interiorX(iV), interiorY(iV)};
    }
  }

  return parm;
}

Vector<double> BFF::dirichletToNeumann(const Vector<double>& uBdy) {
  return Omegab - (Lib.transpose() * Liisolver->solve(Omegai - Lib * uBdy)) - Lbb * uBdy;
}

Vector<double> BFF::neumannToDirichlet(const Vector<double>& kBdy) {
  // Convert Neumann data to Dirichlet data by solving the Poisson equation and reading off values
  ensureHaveLSolver();
  Vector<double> rhs = reassembleVector(Ldecomp, Omegai, Vector<double>(Omegab - kBdy));
  Vector<double> fullSolution = -Lsolver->solve(rhs);
  Vector<double> uBdy, ignore;
  decomposeVector(Ldecomp, fullSolution, ignore, uBdy);
  double uMean = uBdy.mean(); // Ensure that u has mean 0
  for (int i = 0; i < uBdy.size(); i++) uBdy(i) -= uMean;
  return uBdy;
}

std::array<Vector<double>, 2> BFF::computeBoundaryPositions(const Vector<double>& uBdy, const Vector<double>& kBdy) {
  
  geo.requireEdgeLengths();

  double phi = 0;

  std::vector<Eigen::Triplet<double>> Ntriplets;
  Vector<double> targetLength(nBoundary);
  DenseMatrix<double> T(2, nBoundary);

  for (BoundaryLoop b : mesh.boundaryLoops()) {
    for (Halfedge he : b.adjacentHalfedges()) {
      size_t iV = bIdx[he.vertex()];
      Edge e = he.edge();
      GC_SAFETY_ASSERT(iV < nBoundary, "invalid boundary vertex index");

      targetLength(iV) = geo.edgeLengths[e] * exp(0.5 * (uBdy(bIdx[he.tailVertex()]) + uBdy(bIdx[he.tipVertex()])));

      T(0, iV) = cos(phi);
      T(1, iV) = sin(phi);

      Ntriplets.emplace_back(iV, iV, geo.edgeLengths[e]);

      phi += kBdy(bIdx[he.tipVertex()]);
    }
  }

  SparseMatrix<double> Ninv(nBoundary, nBoundary);
  Ninv.setFromTriplets(std::begin(Ntriplets), std::end(Ntriplets));

  Vector<double> roundedLength =
      targetLength - Ninv * T.transpose() * (T * Ninv * T.transpose()).inverse() * (T * targetLength);

  std::array<Vector<double>, 2> bdyPositions{Vector<double>(nBoundary), Vector<double>(nBoundary)};
  double x = 0, y = 0;
  for (BoundaryLoop b : mesh.boundaryLoops()) {
    for (Halfedge he : b.adjacentHalfedges()) {
      size_t iV = bIdx[he.vertex()];

      // flip the param, since we orbit boundary loops in
      // the opposite direction it ends up backward
      bdyPositions[0](iV) = -x;
      bdyPositions[1](iV) = y;

      x += roundedLength(iV) * T(0, iV);
      y += roundedLength(iV) * T(1, iV);
    }
  }

  return bdyPositions;
}

void BFF::ensureHaveLSolver() {
  if (!Lsolver) {
    Lsolver.reset(new PositiveDefiniteSolver<double>(L));
  }
}
} // namespace surface
} // namespace geometrycentral
