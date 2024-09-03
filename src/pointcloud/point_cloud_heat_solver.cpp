#include "geometrycentral/pointcloud/point_cloud_heat_solver.h"


namespace geometrycentral {
namespace pointcloud {


PointCloudHeatSolver::PointCloudHeatSolver(PointCloud& cloud_, PointPositionGeometry& geom_, double tCoef_)
    : tCoef(tCoef_), cloud(cloud_), geom(geom_) {

  // WARNING: these routines mix indexing throughout; ensure the cloud is compressed
  GC_SAFETY_ASSERT(cloud.isCompressed(), "cloud must be compressed");

  geom.requireNeighbors();
  geom.requireTuftedTriangulation();
  geom.tuftedGeom->requireEdgeLengths();
  geom.requireTangentCoordinates();
  geom.requireNeighbors();
  // geom.tuftedGeom->requireCotanLaplacian();
  // geom.requireLaplacian();


  // Compute a timescale
  double meanEdgeLength = 0.;
  for (surface::Edge e : geom.tuftedMesh->edges()) {
    meanEdgeLength += geom.tuftedGeom->edgeLengths[e];
  }
  meanEdgeLength /= geom.tuftedMesh->nEdges();
  shortTime = tCoef * meanEdgeLength * meanEdgeLength;
}

void PointCloudHeatSolver::ensureHaveHeatDistanceWorker() {
  if (heatDistanceWorker != nullptr) return;

  heatDistanceWorker.reset(new surface::HeatMethodDistanceSolver(*geom.tuftedGeom, tCoef));
}

void PointCloudHeatSolver::ensureHaveVectorHeatSolver() {
  if (vectorHeatSolver != nullptr) return;

  geom.requireConnectionLaplacian();
  geom.tuftedGeom->requireVertexLumpedMassMatrix();

  heatDistanceWorker.reset(
      new surface::HeatMethodDistanceSolver(*geom.tuftedGeom, tCoef, false)); // we already have the tufted IDT

  SparseMatrix<double>& Lconn = geom.connectionLaplacian;
  SparseMatrix<double>& massMat = geom.tuftedGeom->vertexLumpedMassMatrix;

  // Build the operator
  SparseMatrix<double> vectorOp = complexToReal(massMat.cast<std::complex<double>>().eval()) + shortTime * Lconn;

  // Note: since tufted Laplacian is always Delaunay, the connection Laplacian is SPD, and we can use Cholesky
  vectorHeatSolver.reset(new PositiveDefiniteSolver<double>(vectorOp));

  geom.unrequireConnectionLaplacian();
  geom.tuftedGeom->unrequireVertexLumpedMassMatrix();
}

// === Heat method for distance
// For distance, we basically just run the mesh version on the tufted triangulation of the point cloud.
PointData<double> PointCloudHeatSolver::computeDistance(const Point& sourcePoint) {
  std::vector<Point> v{sourcePoint};
  return computeDistance(v);
}
PointData<double> PointCloudHeatSolver::computeDistance(const std::vector<Point>& sourcePoints) {
  ensureHaveHeatDistanceWorker();

  std::vector<surface::Vertex> sourceVerts;
  for (Point p : sourcePoints) {
    sourceVerts.push_back(geom.tuftedMesh->vertex(p.getIndex()));
  }
  return PointData<double>(cloud, heatDistanceWorker->computeDistance(sourceVerts).raw());
}

PointData<double> PointCloudHeatSolver::extendScalars(const std::vector<std::tuple<Point, double>>& sources) {
  ensureHaveHeatDistanceWorker();

  GC_SAFETY_ASSERT(sources.size() != 0, "must have at least one source");

  ensureHaveVectorHeatSolver();

  size_t N = cloud.nPoints();
  Vector<double> rhsVals = Vector<double>::Zero(N);
  Vector<double> rhsOnes = Vector<double>::Zero(N);
  for (size_t i = 0; i < sources.size(); i++) {
    size_t ind = std::get<0>(sources[i]).getIndex();
    double val = std::get<1>(sources[i]);
    rhsOnes(ind) = 1.;
    rhsVals(ind) = val;
  }

  Vector<double> interpVals = heatDistanceWorker->heatSolver->solve(rhsVals);
  Vector<double> interpOnes = heatDistanceWorker->heatSolver->solve(rhsOnes);
  Vector<double> resultArr = (interpVals.array() / interpOnes.array());

  PointData<double> result(cloud, resultArr);
  return result;
}

// == Parallel transport

PointData<Vector2> PointCloudHeatSolver::transportTangentVector(const Point& sourcePoint, const Vector2& sourceVector) {

  // std::vector<Point> s{sourcePoint};
  // std::vector<Vector2> v{sourceVector};

  return transportTangentVectors({std::make_tuple(sourcePoint, sourceVector)});
}

PointData<Vector2>
PointCloudHeatSolver::transportTangentVectors(const std::vector<std::tuple<Point, Vector2>>& sources) {

  GC_SAFETY_ASSERT(sources.size() != 0, "must have at least one source");

  ensureHaveVectorHeatSolver();

  // Consturct rhs
  size_t N = cloud.nPoints();
  Vector<std::complex<double>> dirRHS = Vector<std::complex<double>>::Zero(N);
  bool normsAllSame = true;
  double firstNorm = std::get<1>(sources[0]).norm();
  for (size_t i = 0; i < sources.size(); i++) {
    size_t ind = std::get<0>(sources[i]).getIndex();
    Vector2 vec = std::get<1>(sources[i]);
    dirRHS(ind) += std::complex<double>(vec);

    // Check if all norms same
    double thisNorm = vec.norm();
    if (std::abs(firstNorm - thisNorm) > std::fmax(firstNorm, thisNorm) * 1e-10) {
      normsAllSame = false;
    }
  }


  // Transport
  // Result is 2N packed complex
  Vector<double> dirInterpPacked = vectorHeatSolver->solve(complexToReal(dirRHS));

  // Normalize
  Vector<std::complex<double>> dirInterp = Vector<std::complex<double>>::Zero(N);
  for (size_t i = 0; i < N; i++) {
    Vector2 val{dirInterpPacked(2 * i), dirInterpPacked(2 * i + 1)};
    dirInterp(i) = normalizeCutoff(val);
  }

  // Set scale
  if (normsAllSame) {
    // Set the scale by projection
    dirInterp *= firstNorm;
  } else {
    // Interpolate magnitudes
    ensureHaveHeatDistanceWorker();

    Vector<double> rhsNorm = Vector<double>::Zero(N);
    Vector<double> rhsOnes = Vector<double>::Zero(N);
    for (size_t i = 0; i < sources.size(); i++) {
      size_t ind = std::get<0>(sources[i]).getIndex();
      Vector2 vec = std::get<1>(sources[i]);
      rhsOnes(ind) = 1.;
      rhsNorm(ind) = norm(vec);
    }

    Vector<double> interpNorm = heatDistanceWorker->heatSolver->solve(rhsNorm);
    Vector<double> interpOnes = heatDistanceWorker->heatSolver->solve(rhsOnes);

    dirInterp = dirInterp.array() * (interpNorm.array() / interpOnes.array());
  }


  // Copy to a container
  PointData<Vector2> result(cloud);
  for (size_t i = 0; i < N; i++) {
    result[i] = Vector2::fromComplex(dirInterp(i));
  }

  return result;
}

// Compute the logarithmic map from a source point
PointData<Vector2> PointCloudHeatSolver::computeLogMap(const Point& sourcePoint) {
  ensureHaveHeatDistanceWorker();
  ensureHaveVectorHeatSolver();
  geom.requireTangentTransport();

  size_t N = cloud.nPoints();
  PointData<Vector2> logmapResult(cloud);

  { // == Get the direction component from transport of outward vectors

    /* Transport angle strategy
    Vector<std::complex<double>> rhsHorizontal = Vector<std::complex<double>>::Zero(N);
    Vector<std::complex<double>> rhsOutward = Vector<std::complex<double>>::Zero(N);

    rhsHorizontal[sourcePoint.getIndex()] = 1.;

    std::vector<Point>& neighbors = geom.neighbors->neighbors[sourcePoint];
    for (size_t iN = 0; iN < neighbors.size(); iN++) {
      Point neigh = neighbors[iN];
      Vector2 neighOutward = geom.tangentCoordinates[sourcePoint][iN];
      Vector2 outward = geom.tangentTransport[sourcePoint][iN] * neighOutward;
      if(geom.tangentTransportOrientFlip[sourcePoint][iN]) {
        outward = outward.conj();
      }
      rhsOutward[neigh.getIndex()] = outward;
    }

    // Transport
    Vector<std::complex<double>> dirHorizontal = realToComplex(vectorHeatSolver->solve(complexToReal(rhsHorizontal)));
    Vector<std::complex<double>> dirOutward = realToComplex(vectorHeatSolver->solve(complexToReal(rhsOutward)));

    // Store directional component of logmap
    for (size_t i = 0; i < N; i++) {
      logmapResult[i] = normalizeCutoff(Vector2::fromComplex(dirOutward(i) / dirHorizontal(i)));
    }
    */

    // Unit circle strategy
    Vector<double> rhsX = Vector<double>::Zero(N);
    Vector<double> rhsY = Vector<double>::Zero(N);

    std::vector<Point>& neighbors = geom.neighbors->neighbors[sourcePoint];
    for (size_t iN = 0; iN < neighbors.size(); iN++) {
      Point neigh = neighbors[iN];
      Vector2 neighOutward = geom.tangentCoordinates[sourcePoint][iN];
      rhsX[neigh.getIndex()] = neighOutward.x;
      rhsY[neigh.getIndex()] = neighOutward.y;
    }

    // Transport
    Vector<double> dirX = heatDistanceWorker->heatSolver->solve(rhsX);
    Vector<double> dirY = heatDistanceWorker->heatSolver->solve(rhsY);

    // Store directional component of logmap
    for (size_t i = 0; i < N; i++) {
      logmapResult[i] = normalize(Vector2{dirX(i), dirY(i)});
    }
  }

  { // == Get the magnitude component as the distance

    // This is different from what is presented in the Vector Heat Method paper; computing distance entirely on the
    // tufted cover seem to be more robust.

    PointData<double> dist = computeDistance(sourcePoint);

    logmapResult *= dist;
  }

  geom.unrequireTangentTransport();
  return logmapResult;
}

PointData<double> PointCloudHeatSolver::computeSignedDistance(const std::vector<std::vector<Point>>& curves,
                                                              const PointData<Vector3>& cloudNormals,
                                                              const SignedHeatOptions& options) {

  GC_SAFETY_ASSERT(curves.size() != 0, "must have at least one source");

  ensureHaveVectorHeatSolver();

  size_t N = cloud.nPoints();

  // Build curve sources.
  geom.requirePointIndices();
  geom.requireTangentBasis();
  Vector<double> X0 = Vector<double>::Zero(2 * N);
  std::vector<Eigen::Triplet<double>> triplets;
  SparseMatrix<double> C;
  size_t nConstraints = 0;
  PointData<int> curveIndices(cloud, -1);
  for (const auto& curve : curves) {
    size_t nNodes = curve.size();
    // Iterate over curve segments, and accumulate each segment's contribution onto each of its endpoints.
    for (size_t i = 0; i < nNodes - 1; i++) {
      const Point& vA = curve[i];
      const Point& vB = curve[(i + 1) % nNodes];
      Vector3 pA = geom.positions[vA];
      Vector3 pB = geom.positions[vB];
      Vector3 segment = 0.5 * (pB - pA);
      // Project curve segment onto each tangent plane of its endpoints.
      // Use each tangent plane's normal to determine the curve's in-plane normal, expressed in each tangent basis.
      Vector3 surfaceNormalA = cloudNormals[vA].normalize();
      Vector3 surfaceNormalB = cloudNormals[vB].normalize();
      Vector3 xAxisA = geom.tangentBasis[vA][0];
      Vector3 yAxisA = geom.tangentBasis[vA][1];
      Vector3 xAxisB = geom.tangentBasis[vB][0];
      Vector3 yAxisB = geom.tangentBasis[vB][1];
      Vector3 curveTangentA = dot(xAxisA, segment) * xAxisA + dot(yAxisA, segment) * yAxisA;
      Vector3 curveTangentB = dot(xAxisB, segment) * xAxisB + dot(yAxisB, segment) * yAxisB;
      Vector3 inPlaneNormalA = cross(surfaceNormalA, curveTangentA);
      Vector3 inPlaneNormalB = cross(surfaceNormalB, curveTangentB);
      Vector2 curveNormalA = {dot(geom.tangentBasis[vA][0], inPlaneNormalA),
                              dot(geom.tangentBasis[vA][1], inPlaneNormalA)};
      Vector2 curveNormalB = {dot(geom.tangentBasis[vB][0], inPlaneNormalB),
                              dot(geom.tangentBasis[vB][1], inPlaneNormalB)};
      // Project onto each point's local coordinate system.
      size_t vIdxA = geom.pointIndices[vA];
      size_t vIdxB = geom.pointIndices[vB];
      for (int j = 0; j < 2; j++) {
        X0[2 * vIdxA + j] = curveNormalA[j];
        X0[2 * vIdxB + j] = curveNormalB[j];
      }
      if (options.preserveSourceNormals) {
        if (curveIndices[vA] < 0) {
          curveIndices[vA] = nConstraints;
          nConstraints++;
        }
        if (curveIndices[vB] < 0) {
          curveIndices[vB] = nConstraints;
          nConstraints++;
        }
        triplets.emplace_back(curveIndices[vA], 2 * vIdxA, dot(xAxisA, curveTangentA));
        triplets.emplace_back(curveIndices[vA], 2 * vIdxA + 1, dot(yAxisA, curveTangentA));
        triplets.emplace_back(curveIndices[vB], 2 * vIdxB, dot(xAxisB, curveTangentB));
        triplets.emplace_back(curveIndices[vB], 2 * vIdxB + 1, dot(yAxisB, curveTangentB));
      }
    }
  }

  // Diffuse vectors.
  Vector<double> Xt;
  if (!options.preserveSourceNormals) {
    Xt = vectorHeatSolver->solve(X0);
  } else {
    C.resize(nConstraints, 2 * N);
    C.setFromTriplets(triplets.begin(), triplets.end());

    geom.requireConnectionLaplacian();
    geom.tuftedGeom->requireVertexLumpedMassMatrix();
    SparseMatrix<double>& Lconn = geom.connectionLaplacian;
    SparseMatrix<double>& massMat = geom.tuftedGeom->vertexLumpedMassMatrix;
    SparseMatrix<double> vectorOp = complexToReal(massMat.cast<std::complex<double>>().eval()) + shortTime * Lconn;
    geom.unrequireConnectionLaplacian();
    geom.tuftedGeom->unrequireVertexLumpedMassMatrix();

    // Solve system.
    SparseMatrix<double> Z(nConstraints, nConstraints);
    SparseMatrix<double> LHS1 = horizontalStack<double>({vectorOp, C.transpose()});
    SparseMatrix<double> LHS2 = horizontalStack<double>({C, Z});
    SparseMatrix<double> LHS = verticalStack<double>({LHS1, LHS2});
    Vector<double> RHS = Vector<double>::Zero(2 * N + nConstraints);
    RHS.head(2 * N) = X0;
    Vector<double> soln = solveSquare(LHS, RHS);
    Xt = soln.head(2 * N);
  }

  // Normalize vectors & compute divergence on the tufted cover.
  geom.tuftedGeom->requireHalfedgeCotanWeights();
  geom.tuftedGeom->requireVertexIndices();
  Vector<double> divYt = Vector<double>::Zero(N);
  for (surface::Face f : geom.tuftedMesh->faces()) {
    // Compute averaged vector on each face.
    Vector3 grad = {0, 0, 0};
    for (surface::Vertex v : f.adjacentVertices()) {
      size_t idx = geom.tuftedGeom->vertexIndices[v];
      Vector3 X = Xt[2 * idx] * geom.tangentBasis[idx][0] + Xt[2 * idx + 1] * geom.tangentBasis[idx][1];
      grad += X;
    }
    grad = grad.normalize();
    // Accumulate divergence.
    for (surface::Halfedge he : f.adjacentHalfedges()) {
      double cotTheta = geom.tuftedGeom->halfedgeCotanWeights[he];
      size_t vA = geom.tuftedGeom->vertexIndices[he.tailVertex()];
      size_t vB = geom.tuftedGeom->vertexIndices[he.tipVertex()];
      Vector3 heVec = geom.positions[vB] - geom.positions[vA];
      double val = cotTheta * dot(heVec, grad);
      divYt[vA] += val;
      divYt[vB] -= val;
    }
  }
  geom.tuftedGeom->unrequireHalfedgeCotanWeights();
  geom.tuftedGeom->unrequireVertexIndices();
  geom.unrequireTangentBasis();

  // Integrate.
  Vector<double> phi = Vector<double>::Zero(N);
  switch (options.levelSetConstraint) {
  case (LevelSetConstraint::None): {
    ensureHaveHeatDistanceWorker();
    phi = heatDistanceWorker->poissonSolver->solve(divYt);
    // Shift to average zero along (all of) the source geometry.
    double shift = 0.;
    double totalLength = 0.;
    for (const auto& curve : curves) {
      size_t nNodes = curve.size();
      for (size_t i = 0; i < nNodes - 1; i++) {
        const Point& vA = curve[i];
        const Point& vB = curve[(i + 1) % nNodes];
        double length = (geom.positions[vB] - geom.positions[vA]).norm();
        size_t vIdxA = geom.pointIndices[vA];
        size_t vIdxB = geom.pointIndices[vB];
        shift += length * 0.5 * (phi[vIdxA] + phi[vIdxB]);
        totalLength += length;
      }
    }
    shift /= totalLength;
    phi -= shift * Vector<double>::Ones(N);
    break;
  }
  case (LevelSetConstraint::ZeroSet): {
    geom.tuftedGeom->requireCotanLaplacian();
    SparseMatrix<double>& L = geom.tuftedGeom->cotanLaplacian;
    geom.tuftedGeom->unrequireCotanLaplacian();

    if (options.softLevelSetWeight < 0) {
      Vector<bool> setAMembership = Vector<bool>::Ones(N);
      for (const auto& curve : curves) {
        size_t nNodes = curve.size();
        for (size_t i = 0; i < nNodes; i++) {
          setAMembership[geom.pointIndices[curve[i]]] = false;
        }
      }
      size_t nB = N - setAMembership.cast<size_t>().sum();
      Vector<double> bcVals = Vector<double>::Zero(nB);
      BlockDecompositionResult<double> decomp = blockDecomposeSquare(L, setAMembership, true);
      Vector<double> rhsValsA, rhsValsB;
      decomposeVector(decomp, divYt, rhsValsA, rhsValsB);
      Vector<double> combinedRHS = rhsValsA;
      shiftDiagonal(decomp.AA, 1e-8);
      PositiveDefiniteSolver<double> constrainedSolver(decomp.AA);
      Vector<double> Aresult = constrainedSolver.solve(combinedRHS);
      phi = reassembleVector(decomp, Aresult, bcVals);
    } else {
      std::vector<Eigen::Triplet<double>> triplets;
      SparseMatrix<double> A; // constraint matrix
      size_t m = 0;
      for (const auto& curve : curves) {
        size_t nNodes = curve.size();
        bool isClosed = (curve[0] == curve[nNodes - 1]);
        size_t uB = isClosed ? nNodes - 1 : nNodes;
        for (size_t i = 1; i < uB; i++) {
          triplets.emplace_back(m, geom.pointIndices[curve[i]], 1);
          m++;
        }
      }
      A.resize(m, N);
      A.setFromTriplets(triplets.begin(), triplets.end());
      SparseMatrix<double> LHS = L + options.softLevelSetWeight * A.transpose() * A;
      phi = solvePositiveDefinite(LHS, divYt);
    }
    break;
  }
  case (LevelSetConstraint::Multiple): {
    geom.tuftedGeom->requireCotanLaplacian();
    SparseMatrix<double>& L = geom.tuftedGeom->cotanLaplacian;
    geom.tuftedGeom->unrequireCotanLaplacian();

    std::vector<Eigen::Triplet<double>> triplets;
    SparseMatrix<double> A; // constraint matrix
    size_t m = 0;
    for (const auto& curve : curves) {
      size_t nNodes = curve.size();
      bool isClosed = (curve[0] == curve[nNodes - 1]);
      size_t uB = isClosed ? nNodes - 1 : nNodes;
      const Point& p0 = curve[0];
      for (size_t i = 1; i < uB; i++) {
        triplets.emplace_back(m, geom.pointIndices[curve[i]], -1);
        triplets.emplace_back(m, geom.pointIndices[p0], 1);
        m++;
      }
    }
    A.resize(m, N);
    A.setFromTriplets(triplets.begin(), triplets.end());
    if (options.softLevelSetWeight < 0) {
      SparseMatrix<double> Z(m, m);
      SparseMatrix<double> LHS1 = horizontalStack<double>({L, A.transpose()});
      SparseMatrix<double> LHS2 = horizontalStack<double>({A, Z});
      SparseMatrix<double> LHS = verticalStack<double>({LHS1, LHS2});
      Vector<double> RHS = Vector<double>::Zero(N + m);
      RHS.head(N) = divYt;
      Vector<double> soln = solveSquare(LHS, RHS);
      phi = soln.head(N);
    } else {
      SparseMatrix<double> LHS = L + options.softLevelSetWeight * A.transpose() * A;
      phi = solvePositiveDefinite(LHS, divYt);
    }
    break;
  }
  }
  geom.unrequirePointIndices();
  PointData<double> GSDF(cloud, phi);

  return GSDF;
}

} // namespace pointcloud
} // namespace geometrycentral
