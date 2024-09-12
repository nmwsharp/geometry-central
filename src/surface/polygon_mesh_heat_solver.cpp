#include "geometrycentral/surface/polygon_mesh_heat_solver.h"

namespace geometrycentral {
namespace surface {

PolygonMeshHeatSolver::PolygonMeshHeatSolver(EmbeddedGeometryInterface& geom_, double tCoef_)
    : tCoef(tCoef_), mesh(geom_.mesh), geom(geom_)

{
  // Though de Goes et al. suggest using the maximum length over all polygon diagonals, using the mean length seems more
  // accurate.
  geom.requireEdgeLengths();
  double meanEdgeLength = 0.;
  for (Edge e : mesh.edges()) meanEdgeLength += geom.edgeLengths[e];
  meanEdgeLength /= mesh.nEdges();
  shortTime = tCoef * meanEdgeLength * meanEdgeLength;
  geom.unrequireEdgeLengths();

  geom.requirePolygonVertexLumpedMassMatrix();
  geom.requirePolygonLaplacian();
  massMat = geom.polygonVertexLumpedMassMatrix;
  laplaceMat = geom.polygonLaplacian;
  geom.unrequirePolygonVertexLumpedMassMatrix();
  geom.unrequirePolygonLaplacian();
}

void PolygonMeshHeatSolver::ensureHaveScalarHeatSolver() {
  if (scalarHeatSolver != nullptr) return;

  SparseMatrix<double> heatOp = massMat + shortTime * laplaceMat;
  scalarHeatSolver.reset(new PositiveDefiniteSolver<double>(heatOp));
}

void PolygonMeshHeatSolver::ensureHaveVectorHeatSolver() {
  if (vectorHeatSolver != nullptr) return;

  geom.requirePolygonVertexConnectionLaplacian();

  SparseMatrix<std::complex<double>>& L = geom.polygonVertexConnectionLaplacian;
  SparseMatrix<std::complex<double>> vectorOp = massMat.cast<std::complex<double>>() + shortTime * L;

  vectorHeatSolver.reset(new PositiveDefiniteSolver<std::complex<double>>(vectorOp));

  geom.unrequirePolygonVertexConnectionLaplacian();
}


void PolygonMeshHeatSolver::ensureHavePoissonSolver() {
  if (poissonSolver != nullptr) return;

  poissonSolver.reset(new PositiveDefiniteSolver<double>(laplaceMat));
}

VertexData<double> PolygonMeshHeatSolver::computeDistance(const Vertex& sourceVert) {
  std::vector<Vertex> v{sourceVert};
  return computeDistance(v);
}

VertexData<double> PolygonMeshHeatSolver::computeDistance(const std::vector<Vertex>& sourceVerts) {
  GC_SAFETY_ASSERT(sourceVerts.size() != 0, "must have at least one source");

  geom.requirePolygonGradientMatrix();
  geom.requirePolygonDivergenceMatrix();
  geom.requireVertexIndices();

  // Flow heat.
  size_t V = mesh.nVertices();
  size_t F = mesh.nFaces();
  Vector<double> rho = Vector<double>::Zero(V);
  for (Vertex v : sourceVerts) {
    size_t vIdx = geom.vertexIndices[v];
    rho[vIdx] += 1;
  }

  ensureHaveScalarHeatSolver();
  Vector<double> X = scalarHeatSolver->solve(rho);

  // Normalize gradient.
  Vector<double> Y = geom.polygonGradientMatrix * X; // 3|F|
  for (size_t i = 0; i < F; i++) {
    Vector3 g = {Y[3 * i], Y[3 * i + 1], Y[3 * i + 2]};
    g /= g.norm();
    for (int j = 0; j < 3; j++) Y[3 * i + j] = g[j];
  }

  // Integrate.
  ensureHavePoissonSolver();
  Vector<double> div = -geom.polygonDivergenceMatrix * Y;
  Vector<double> distances = poissonSolver->solve(div);

  // Shift solution.
  distances -= distances.minCoeff() * Vector<double>::Ones(V);

  geom.unrequirePolygonGradientMatrix();
  geom.unrequirePolygonDivergenceMatrix();
  geom.unrequireVertexIndices();

  return VertexData<double>(mesh, distances);
}

VertexData<double> PolygonMeshHeatSolver::extendScalars(const std::vector<std::tuple<Vertex, double>>& sources) {
  GC_SAFETY_ASSERT(sources.size() != 0, "must have at least one source");

  ensureHaveScalarHeatSolver();

  size_t V = mesh.nVertices();
  Vector<double> rhsVals = Vector<double>::Zero(V);
  Vector<double> rhsOnes = Vector<double>::Zero(V);
  for (size_t i = 0; i < sources.size(); i++) {
    size_t ind = std::get<0>(sources[i]).getIndex();
    double val = std::get<1>(sources[i]);
    rhsVals[ind] = val;
    rhsOnes[ind] = 1.;
  }

  Vector<double> interpVals = scalarHeatSolver->solve(rhsVals);
  Vector<double> interpOnes = scalarHeatSolver->solve(rhsOnes);
  Vector<double> resultArr = (interpVals.array() / interpOnes.array());

  VertexData<double> result(mesh, resultArr);
  return result;
}

VertexData<Vector2> PolygonMeshHeatSolver::transportTangentVector(const Vertex& sourceVert,
                                                                  const Vector2& sourceVector) {

  return transportTangentVectors({std::make_tuple(sourceVert, sourceVector)});
}

VertexData<Vector2>
PolygonMeshHeatSolver::transportTangentVectors(const std::vector<std::tuple<Vertex, Vector2>>& sources) {
  GC_SAFETY_ASSERT(sources.size() != 0, "must have at least one source");

  ensureHaveVectorHeatSolver();
  geom.requireVertexIndices();

  // Construct rhs
  size_t V = mesh.nVertices();
  Vector<std::complex<double>> dirRHS = Vector<std::complex<double>>::Zero(V);
  std::vector<std::tuple<Vertex, double>> magnitudeSources;
  bool normsAllSame = true;
  double firstNorm = std::get<1>(sources[0]).norm();
  for (size_t i = 0; i < sources.size(); i++) {
    Vertex v = std::get<0>(sources[i]);
    size_t ind = geom.vertexIndices[v];
    Vector2 vec = std::get<1>(sources[i]);
    dirRHS(ind) += std::complex<double>(vec);
    magnitudeSources.emplace_back(v, vec.norm());

    // Check if all norms same
    double thisNorm = vec.norm();
    if (std::abs(firstNorm - thisNorm) > std::fmax(firstNorm, thisNorm) * 1e-10) {
      normsAllSame = false;
    }
  }
  geom.unrequireVertexIndices();

  // Transport
  Vector<std::complex<double>> vecSolution = vectorHeatSolver->solve(dirRHS);

  // Set scale
  VertexData<Vector2> solution(mesh);
  if (normsAllSame) {
    vecSolution = (vecSolution.array() / vecSolution.array().abs()) * firstNorm;
    for (size_t i = 0; i < V; i++) solution[i] = Vector2::fromComplex(vecSolution[i]);
  } else {
    VertexData<double> interpMags = extendScalars(magnitudeSources);
    for (size_t i = 0; i < V; i++) {
      Vector2 dir = Vector2::fromComplex(vecSolution[i]).normalize();
      solution[i] = dir * interpMags[i];
    }
  }

  return solution;
}

VertexData<double> PolygonMeshHeatSolver::computeSignedDistance(const std::vector<std::vector<Vertex>>& curves,
                                                                const LevelSetConstraint& levelSetConstraint) {

  size_t V = mesh.nVertices();
  Vector<std::complex<double>> X0 = Vector<std::complex<double>>::Zero(V);
  geom.requireVertexIndices();
  geom.requireVertexTangentBasis();
  geom.requireVertexNormals();
  for (const auto& curve : curves) buildSignedCurveSource(curve, X0);
  geom.unrequireVertexNormals();
  GC_SAFETY_ASSERT(X0.norm() != 0, "must have at least one source");

  ensureHaveVectorHeatSolver();
  Vector<std::complex<double>> Xt = vectorHeatSolver->solve(X0);

  // Average onto faces, and normalize.
  size_t F = mesh.nFaces();
  Vector<double> Y(3 * F);
  geom.requireFaceIndices();
  for (Face f : mesh.faces()) {
    Vector3 Yf = {0, 0, 0};
    for (Vertex v : f.adjacentVertices()) {
      size_t vIdx = geom.vertexIndices[v];
      Yf += std::real(Xt[vIdx]) * geom.vertexTangentBasis[v][0];
      Yf += std::imag(Xt[vIdx]) * geom.vertexTangentBasis[v][1];
    }
    Yf /= Yf.norm();
    size_t fIdx = geom.faceIndices[f];
    for (int j = 0; j < 3; j++) {
      Y[3 * fIdx + j] = Yf[j];
    }
  }
  geom.unrequireFaceIndices();
  geom.unrequireVertexTangentBasis();

  geom.requirePolygonDivergenceMatrix();
  Vector<double> divYt = geom.polygonDivergenceMatrix * Y;
  geom.unrequirePolygonDivergenceMatrix();

  Vector<double> phi;
  if (levelSetConstraint == LevelSetConstraint::None) {
    ensureHavePoissonSolver();
    phi = poissonSolver->solve(divYt);
    // Shift so that average value is zero along input curves.
    double shift = computeAverageValue(curves, phi);
    phi -= shift * Vector<double>::Ones(V);
  } else if (levelSetConstraint == LevelSetConstraint::ZeroSet) {
    Vector<bool> setAMembership = Vector<bool>::Ones(V);
    for (const auto& curve : curves) {
      for (const Vertex& v : curve) {
        setAMembership[geom.vertexIndices[v]] = false;
      }
    }
    int nB = V - setAMembership.cast<int>().sum();
    Vector<double> bcVals = Vector<double>::Zero(nB);
    BlockDecompositionResult<double> decomp = blockDecomposeSquare(laplaceMat, setAMembership, true);
    Vector<double> rhsValsA, rhsValsB;
    decomposeVector(decomp, divYt, rhsValsA, rhsValsB);
    Vector<double> combinedRHS = rhsValsA;
    Vector<double> Aresult = solvePositiveDefinite(decomp.AA, combinedRHS);
    phi = reassembleVector(decomp, Aresult, bcVals);
  } else if (levelSetConstraint == LevelSetConstraint::Multiple) {
    std::vector<Eigen::Triplet<double>> triplets;
    SparseMatrix<double> A;
    size_t m = 0;
    for (const auto& curve : curves) {
      size_t nNodes = curve.size();
      if (nNodes < 2) continue;
      size_t v0 = geom.vertexIndices[curve[0]];
      for (size_t i = 1; i < nNodes; i++) {
        size_t vIdx = geom.vertexIndices[curve[i]];
        triplets.emplace_back(m, v0, 1);
        triplets.emplace_back(m, vIdx, -1);
        m++;
      }
    }
    A.resize(m, V);
    A.setFromTriplets(triplets.begin(), triplets.end());
    SparseMatrix<double> Z(m, m);
    SparseMatrix<double> LHS1 = horizontalStack<double>({laplaceMat, A.transpose()});
    SparseMatrix<double> LHS2 = horizontalStack<double>({A, Z});
    SparseMatrix<double> LHS = verticalStack<double>({LHS1, LHS2});
    Vector<double> RHS = Vector<double>::Zero(V + m);
    RHS.head(V) = divYt;
    Vector<double> soln = solveSquare(LHS, RHS);
    phi = soln.head(V);
    double shift = computeAverageValue(curves, phi);
    phi -= shift * Vector<double>::Ones(V);
  }

  geom.unrequireVertexIndices();

  return VertexData<double>(mesh, phi);
}

void PolygonMeshHeatSolver::buildSignedCurveSource(const std::vector<Vertex>& curve,
                                                   Vector<std::complex<double>>& X0) const {

  // Encode curve input by expressing curve normals in vertex tangent spaces.
  size_t nNodes = curve.size();
  for (size_t i = 0; i < nNodes - 1; i++) {
    Vertex vA = curve[i];
    Vertex vB = curve[i + 1];
    size_t vIdxA = geom.vertexIndices[vA];
    size_t vIdxB = geom.vertexIndices[vB];
    Vector3 pA = geom.vertexPositions[vA];
    Vector3 pB = geom.vertexPositions[vB];
    Vector3 tangent = pB - pA;
    Vector3 vnA = geom.vertexNormals[vA];
    Vector3 vnB = geom.vertexNormals[vB];
    Vector3 xAxisA = geom.vertexTangentBasis[vA][0];
    Vector3 yAxisA = geom.vertexTangentBasis[vA][1];
    Vector3 xAxisB = geom.vertexTangentBasis[vB][0];
    Vector3 yAxisB = geom.vertexTangentBasis[vB][1];
    Vector3 tangentA = dot(xAxisA, tangent) * xAxisA + dot(yAxisA, tangent) * yAxisA;
    Vector3 tangentB = dot(xAxisB, tangent) * xAxisB + dot(yAxisB, tangent) * yAxisB;
    Vector3 normalA = cross(vnA, tangentA);
    Vector3 normalB = cross(vnB, tangentB);
    X0[vIdxA] += std::complex<double>(dot(geom.vertexTangentBasis[vA][0], normalA),
                                      dot(geom.vertexTangentBasis[vA][1], normalA));
    X0[vIdxB] += std::complex<double>(dot(geom.vertexTangentBasis[vB][0], normalB),
                                      dot(geom.vertexTangentBasis[vB][1], normalB));
  }
}

double PolygonMeshHeatSolver::computeAverageValue(const std::vector<std::vector<Vertex>>& curves,
                                                  const Vector<double>& u) {

  geom.requireVertexIndices();

  double shift = 0.;
  double totalLength = 0.;
  for (const auto& curve : curves) {
    size_t nNodes = curve.size();
    for (size_t i = 0; i < nNodes - 1; i++) {
      Vertex vA = curve[i];
      Vertex vB = curve[i + 1];
      size_t vIdxA = geom.vertexIndices[vA];
      size_t vIdxB = geom.vertexIndices[vB];
      Vector3 pA = geom.vertexPositions[vA];
      Vector3 pB = geom.vertexPositions[vB];
      double length = (pB - pA).norm();
      shift += 0.5 * length * (u[vIdxA] + u[vIdxB]);
      totalLength += length;
    }
  }
  shift /= totalLength;

  geom.unrequireVertexIndices();

  return shift;
}

} // namespace surface
} // namespace geometrycentral