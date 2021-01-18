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

  heatDistanceWorker.reset(new surface::HeatMethodDistanceSolver(*geom.tuftedGeom, tCoef));

  SparseMatrix<std::complex<double>>& Lconn = geom.connectionLaplacian;
  SparseMatrix<double>& massMat = geom.tuftedGeom->vertexLumpedMassMatrix;

  // Build the operator
  SparseMatrix<std::complex<double>> vectorOp = massMat.cast<std::complex<double>>() + shortTime * Lconn;

  // Note: since tufted Laplacian is always Delaunay, the connection Laplacian is SPD, and we can use Cholesky
  vectorHeatSolver.reset(new PositiveDefiniteSolver<std::complex<double>>(vectorOp));

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

// == Parallel transport

PointData<Vector2> PointCloudHeatSolver::computeParallelTransport(const Point& sourcePoint,
                                                                  const Vector2& sourceVector) {

  std::vector<Point> s{sourcePoint};
  std::vector<Vector2> v{sourceVector};

  return computeParallelTransport(s, v);
}

PointData<Vector2> PointCloudHeatSolver::computeParallelTransport(const std::vector<Point>& sourcePoints,
                                                                  const std::vector<Vector2>& sourceVectors) {

  if (sourcePoints.size() != sourceVectors.size()) {
    throw std::runtime_error("source points and vectors must be same size");
  }
  if (sourcePoints.size() == 0) {
    throw std::runtime_error("must have at least one source");
  }

  ensureHaveVectorHeatSolver();

  // Consturct rhs
  size_t N = cloud.nPoints();
  Vector<std::complex<double>> dirRHS = Vector<std::complex<double>>::Zero(N);
  bool normsAllSame = true;
  double firstNorm = sourceVectors[0].norm();
  for (size_t i = 0; i < sourcePoints.size(); i++) {
    size_t ind = sourcePoints[i].getIndex();
    Vector2 vec = sourceVectors[i];
    dirRHS(ind) += std::complex<double>(vec);

    // Check if all norms same
    double thisNorm = sourceVectors[i].norm();
    if (std::abs(firstNorm - thisNorm) > std::fmax(firstNorm, thisNorm) * 1e-10) {
      normsAllSame = false;
    }
  }


  // Transport
  Vector<std::complex<double>> dirInterp = vectorHeatSolver->solve(dirRHS);

  // Normalize
  for (size_t i = 0; i < N; i++) {
    dirInterp(i) = normalizeCutoff(Vector2::fromComplex(dirInterp(i)));
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

    Vector<std::complex<double>> rhsHorizontal = Vector<std::complex<double>>::Zero(N);
    Vector<std::complex<double>> rhsOutward = Vector<std::complex<double>>::Zero(N);

    rhsHorizontal[sourcePoint.getIndex()] = 1.;

    std::vector<Point>& neighbors = geom.neighbors->neighbors[sourcePoint];
    for (size_t iN = 0; iN < neighbors.size(); iN++) {
      Point neigh = neighbors[iN];
      Vector2 neighOutward = geom.tangentCoordinates[sourcePoint][iN];
      Vector2 outward = geom.tangentTransport[sourcePoint][iN] * neighOutward;
      rhsOutward[neigh.getIndex()] = outward;
    }

    // Transport
    Vector<std::complex<double>> dirHorizontal = vectorHeatSolver->solve(rhsHorizontal);
    Vector<std::complex<double>> dirOutward = vectorHeatSolver->solve(rhsOutward);

    // Store directional component of logmap
    for (size_t i = 0; i < N; i++) {
      logmapResult[i] = normalizeCutoff(Vector2::fromComplex(dirOutward(i) / dirHorizontal(i)));
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
