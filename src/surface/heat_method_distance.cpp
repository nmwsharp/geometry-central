#include "geometrycentral/surface/heat_method_distance.h"

#include "geometrycentral/surface/intrinsic_mollification.h"
#include "geometrycentral/surface/simple_idt.h"
#include "geometrycentral/surface/tufted_laplacian.h"


namespace geometrycentral {
namespace surface {

VertexData<double> heatMethodDistance(IntrinsicGeometryInterface& geom, Vertex v) {
  return HeatMethodDistanceSolver(geom).computeDistance(v);
}

HeatMethodDistanceSolver::HeatMethodDistanceSolver(IntrinsicGeometryInterface& geom_, double tCoef_,
                                                   bool useRobustLaplacian_)
    : tCoef(tCoef_), useRobustLaplacian(useRobustLaplacian_), mesh(geom_.mesh), geom(geom_) {

  // === Build & factor the linear systems
  if (useRobustLaplacian) {
    geom.requireEdgeLengths();

    // Build operators using robust Laplacian (see [Sharp & Crane. A Laplacian for Nonmanifold Triangle Meshes. SGP
    // 2020]) NOTE: we build explicitly here rather than just calling buildTuftedLaplacian() in order to use an
    // intrinsic geometry.

    // Create a copy of the mesh / geometry to operate on, and build the mollified (tufted if nonmanifold) cover
    EdgeData<double> tuftedEdgeLengths;
    if (mesh.usesImplicitTwin()) {
      tuftedMesh = mesh.copy();
      tuftedEdgeLengths = geom.edgeLengths.reinterpretTo(*tuftedMesh);
    } else {
      tuftedMesh = mesh.copyToSurfaceMesh();
      tuftedEdgeLengths = geom.edgeLengths.reinterpretTo(*tuftedMesh);
      buildIntrinsicTuftedCover(*tuftedMesh, tuftedEdgeLengths);
    }
    mollifyIntrinsic(*tuftedMesh, tuftedEdgeLengths, 1e-5);
    size_t nFlips = flipToDelaunay(*tuftedMesh, tuftedEdgeLengths);
    tuftedIntrinsicGeom.reset(new EdgeLengthGeometry(*tuftedMesh, tuftedEdgeLengths));
  }

  // Compute mean edge length and set shortTime
  getGeom().requireEdgeLengths();
  double meanEdgeLength = 0.;
  for (Edge e : getMesh().edges()) {
    meanEdgeLength += getGeom().edgeLengths[e];
  }
  meanEdgeLength /= getMesh().nEdges();
  shortTime = tCoef * meanEdgeLength * meanEdgeLength;


  // Mass matrix
  getGeom().requireVertexLumpedMassMatrix();
  SparseMatrix<double>& M = getGeom().vertexLumpedMassMatrix;

  // Laplacian
  getGeom().requireCotanLaplacian();
  SparseMatrix<double>& L = getGeom().cotanLaplacian;

  // Heat operator
  SparseMatrix<double> heatOp = M + shortTime * L;
  heatSolver.reset(new PositiveDefiniteSolver<double>(heatOp));

  // Poisson solver
  // NOTE: In theory, it should not be necessary to shift the Laplacian: cotan-Laplace is always PSD. However, when the
  // matrix is only positive SEMIdefinite, some solvers may not work (ie Eigen's Cholesky solver doesn't work, but
  // Suitesparse does).
  SparseMatrix<double> Ls = L + 1e-6 * identityMatrix<double>(mesh.nVertices());
  poissonSolver.reset(new PositiveDefiniteSolver<double>(Ls));

  getGeom().unrequireEdgeLengths();
  getGeom().unrequireCotanLaplacian();
  getGeom().unrequireVertexLumpedMassMatrix();
}

SurfaceMesh& HeatMethodDistanceSolver::getMesh() { return useRobustLaplacian ? *tuftedMesh : mesh; }
IntrinsicGeometryInterface& HeatMethodDistanceSolver::getGeom() {
  return useRobustLaplacian ? *tuftedIntrinsicGeom : geom;
}

VertexData<double> HeatMethodDistanceSolver::computeDistance(const Vertex& sourceVert) {
  // call general version
  return computeDistance({SurfacePoint(sourceVert)});
}
VertexData<double> HeatMethodDistanceSolver::computeDistance(const std::vector<Vertex>& sourceVerts) {
  std::vector<SurfacePoint> surfacePoints;
  for (Vertex v : sourceVerts) {
    surfacePoints.emplace_back(v);
  }

  // call general version
  return computeDistance(surfacePoints);
}

VertexData<double> HeatMethodDistanceSolver::computeDistance(const SurfacePoint& sourcePoint) {
  // call general version
  return computeDistance(std::vector<SurfacePoint>{sourcePoint});
}

VertexData<double> HeatMethodDistanceSolver::computeDistance(const std::vector<SurfacePoint>& sourcePoints) {
  getGeom().requireHalfedgeCotanWeights();
  getGeom().requireHalfedgeVectorsInFace();
  getGeom().requireEdgeLengths();
  getGeom().requireVertexIndices();
  getGeom().requireVertexDualAreas();
  geom.requireEdgeLengths();
  geom.requireVertexIndices();

  // === Build RHS
  VertexData<double> rhs(mesh, 0.);
  for (const SurfacePoint& p : sourcePoints) {
    SurfacePoint faceP = p.inSomeFace();

    // Set initial values at the three adjacent vertices
    Halfedge he = faceP.face.halfedge();
    rhs[he.vertex()] += faceP.faceCoords.x;
    rhs[he.next().vertex()] += faceP.faceCoords.y;
    rhs[he.next().next().vertex()] += faceP.faceCoords.z;
  }
  Vector<double> rhsVec = rhs.toVector();

  Vector<double> distVec = computeDistanceRHS(rhsVec);

  // ===  Shift distance to put zero at the source set

  // Helper to measure distance between two points, given their barycentric coordinates
  auto baryDist = [&](Vector3 b1, Vector3 b2, const std::array<double, 3>& edgeLengths) {
    // Shindler & Chen 2012, Barycentric Coordinates in Olympiad Geometry, Section 3.2
    Vector3 bVec = b2 - b1;
    double d2 = 0;
    for (int i = 0; i < 3; i++) {
      d2 += edgeLengths[i] * edgeLengths[i] * bVec[i] * bVec[(i + 1) % 3];
    }
    if (!(d2 <= 0)) d2 = 0.; // ensure it's negative so the sqrt below succeed
    return std::sqrt(-d2);
  };

  double distDiffAtSource = 0;
  double weightSum = 0;
  for (const SurfacePoint& p : sourcePoints) {
    SurfacePoint faceP = p.inSomeFace();

    Halfedge he0 = faceP.face.halfedge();
    std::array<double, 3> edgeLengths{geom.edgeLengths[he0.edge()], geom.edgeLengths[he0.next().edge()],
                                      geom.edgeLengths[he0.next().next().edge()]};


    int i = 0;
    for (Halfedge he : faceP.face.adjacentHalfedges()) {

      Vector3 targetP = Vector3::zero();
      targetP[i] = 1.;

      double expectedDistAtVert = baryDist(faceP.faceCoords, targetP, edgeLengths);
      double actDistAtVert = distVec[geom.vertexIndices[he.vertex()]];

      double w = faceP.faceCoords[i];
      distDiffAtSource += (actDistAtVert - expectedDistAtVert) * w;
      weightSum += w;

      i++;
    }
  }
  distDiffAtSource /= weightSum;

  double shift = -distDiffAtSource;
  distVec = distVec.array() + shift;

  getGeom().unrequireHalfedgeCotanWeights();
  getGeom().unrequireHalfedgeVectorsInFace();
  getGeom().unrequireEdgeLengths();
  getGeom().unrequireVertexIndices();
  getGeom().unrequireVertexDualAreas();
  geom.unrequireEdgeLengths();
  geom.unrequireVertexIndices();

  return VertexData<double>(mesh, distVec);
}

Vector<double> HeatMethodDistanceSolver::computeDistanceRHS(const Vector<double>& rhsVec) {
  getGeom().requireHalfedgeCotanWeights();
  getGeom().requireHalfedgeVectorsInFace();
  getGeom().requireEdgeLengths();
  getGeom().requireVertexIndices();
  getGeom().requireVertexDualAreas();

  // === Solve heat
  Vector<double> heatVec = heatSolver->solve(rhsVec);

  // === Normalize in each face and evaluate divergence
  Vector<double> divergenceVec = Vector<double>::Zero(mesh.nVertices());
  for (Face f : getMesh().faces()) {

    Vector2 gradUDir = Vector2::zero(); // warning, wrong magnitude because we don't care
    for (Halfedge he : f.adjacentHalfedges()) {
      Vector2 ePerp = getGeom().halfedgeVectorsInFace[he.next()].rotate90();
      gradUDir += ePerp * heatVec(getGeom().vertexIndices[he.vertex()]);
    }

    gradUDir = gradUDir.normalizeCutoff();

    for (Halfedge he : f.adjacentHalfedges()) {
      double val = getGeom().halfedgeCotanWeights[he] * dot(getGeom().halfedgeVectorsInFace[he], gradUDir);
      divergenceVec[getGeom().vertexIndices[he.tailVertex()]] += val;
      divergenceVec[getGeom().vertexIndices[he.tipVertex()]] += -val;
    }
  }

  // === Integrate divergence to get distance
  Vector<double> distVec = poissonSolver->solve(divergenceVec);

  getGeom().unrequireHalfedgeVectorsInFace();
  getGeom().unrequireHalfedgeCotanWeights();
  getGeom().unrequireEdgeLengths();
  getGeom().unrequireVertexIndices();
  getGeom().unrequireVertexDualAreas();

  return distVec;
}


} // namespace surface
} // namespace geometrycentral
