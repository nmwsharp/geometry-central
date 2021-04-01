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

} // namespace pointcloud
} // namespace geometrycentral
