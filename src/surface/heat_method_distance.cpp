#include "geometrycentral/surface/heat_method_distance.h"


namespace geometrycentral {
namespace surface {

VertexData<double> heatMethodDistance(IntrinsicGeometryInterface& geom, Vertex v) {
	return HeatMethodDistanceSolver(geom).computeDistance(v);
}

HeatMethodDistanceSolver::HeatMethodDistanceSolver(IntrinsicGeometryInterface& geom_, double tCoef_)
    : tCoef(tCoef_), mesh(geom_.mesh), geom(geom_)

{

  // Compute mean edge length and set shortTime
  geom.requireEdgeLengths();
  double meanEdgeLength = 0.;
  for (Edge e : mesh.edges()) {
    meanEdgeLength += geom.edgeLengths[e];
  }
  meanEdgeLength /= mesh.nEdges();
  shortTime = tCoef * meanEdgeLength * meanEdgeLength;


  // === Build & factor the linear systems

  // Mass matrix
  geom.requireVertexLumpedMassMatrix();
  SparseMatrix<double>& M = geom.vertexLumpedMassMatrix;

  // Laplacian
  geom.requireCotanLaplacian();
  SparseMatrix<double>& L = geom.cotanLaplacian;

  // Heat operator
  SparseMatrix<double> heatOp = M + shortTime * L;
  heatSolver.reset(new PositiveDefiniteSolver<double>(heatOp));

  // Poisson solver
  poissonSolver.reset(new PositiveDefiniteSolver<double>(L));


  geom.unrequireEdgeLengths();
  geom.unrequireCotanLaplacian();
  geom.unrequireVertexLumpedMassMatrix();
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
  geom.requireHalfedgeCotanWeights();
  geom.requireHalfedgeVectorsInFace();
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


  // === Solve heat
  Vector<double> heatVec = heatSolver->solve(rhsVec);


  // === Normalize in each face and evaluate divergence
  Vector<double> divergenceVec = Vector<double>::Zero(mesh.nVertices());
  for (Face f : mesh.faces()) {

    Vector2 gradUDir = Vector2::zero(); // warning, wrong magnitude because we don't care
    for (Halfedge he : f.adjacentHalfedges()) {
      Vector2 ePerp = geom.halfedgeVectorsInFace[he.next()].rotate90();
      gradUDir += ePerp * heatVec(geom.vertexIndices[he.vertex()]);
    }

    gradUDir = gradUDir.normalize();

    for (Halfedge he : f.adjacentHalfedges()) {
      double val = geom.halfedgeCotanWeights[he] * dot(geom.halfedgeVectorsInFace[he], gradUDir);
      divergenceVec[geom.vertexIndices[he.vertex()]] += val;
      divergenceVec[geom.vertexIndices[he.twin().vertex()]] += -val;
    }
  }

  // === Integrate divergence to get distance
  Vector<double> distVec = poissonSolver->solve(divergenceVec);


  // ===  Shift distance to put zero at the source set

  // Helper to measure distance between two points, given their barycentric coordinates
  auto baryDist = [&](Vector3 b1, Vector3 b2, const std::array<double, 3>& edgeLengths) {
    // Shindler & Chen 2012, Barycentric Coordinates in Olympiad Geometry, Section 3.2
    Vector3 bVec = b2 - b1;
    double d2 = 0;
    for (int i = 0; i < 3; i++) {
      d2 += edgeLengths[i] * edgeLengths[i] * bVec[i] * bVec[(i + 1) % 3];
    }
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

  geom.unrequireHalfedgeVectorsInFace();
  geom.unrequireHalfedgeCotanWeights();
  geom.requireEdgeLengths();
  geom.unrequireVertexIndices();

  return VertexData<double>(mesh, distVec);
}


} // namespace surface
} // namespace geometrycentral
