#include "geometrycentral/surface/vector_heat_method.h"
#include "geometrycentral/numerical/linear_algebra_utilities.h"
#include "geometrycentral/utilities/vector2.h"

namespace geometrycentral {
namespace surface {

VectorHeatMethodSolver::VectorHeatMethodSolver(IntrinsicGeometryInterface& geom_, double tCoef_)
    : tCoef(tCoef_), mesh(geom_.mesh), geom(geom_)

{
  geom.requireEdgeLengths();
  geom.requireVertexLumpedMassMatrix();

  // Compute mean edge length and set shortTime
  double meanEdgeLength = 0.;
  for (Edge e : mesh.edges()) {
    meanEdgeLength += geom.edgeLengths[e];
  }
  meanEdgeLength /= mesh.nEdges();
  shortTime = tCoef * meanEdgeLength * meanEdgeLength;

  // We always want the mass matrix
  massMat = geom.vertexLumpedMassMatrix;

  geom.unrequireVertexLumpedMassMatrix();
  geom.unrequireEdgeLengths();
}


void VectorHeatMethodSolver::ensureHaveScalarHeatSolver() {
  if (scalarHeatSolver != nullptr) return;

  // Get the ingredients
  geom.requireCotanLaplacian();
  SparseMatrix<double>& L = geom.cotanLaplacian;

  // Build the operator
  SparseMatrix<double> heatOp = massMat + shortTime * L;
  scalarHeatSolver.reset(new PositiveDefiniteSolver<double>(heatOp));

  geom.unrequireCotanLaplacian();
}

void VectorHeatMethodSolver::ensureHaveVectorHeatSolver() {
  if (vectorHeatSolver != nullptr) return;

  // Get the ingredients
  geom.requireVertexConnectionLaplacian();
  SparseMatrix<std::complex<double>>& Lconn = geom.vertexConnectionLaplacian;

  // Build the operator
  SparseMatrix<std::complex<double>> vectorOp = massMat.cast<std::complex<double>>() + shortTime * Lconn;

  // Check the Delaunay condition. If the mesh is Delaunay, then vectorOp is SPD, and we can use a
  // PositiveDefiniteSolver. Otherwise, we must use a SquareSolver
  geom.requireEdgeCotanWeights();
  bool isDelaunay = true;
  double minCotanWeight = std::numeric_limits<double>::max();
  for (Edge e : mesh.edges()) {
    minCotanWeight = std::min(minCotanWeight, geom.edgeCotanWeights[e]);
    if (geom.edgeCotanWeights[e] < -1e-6) {
      isDelaunay = false;
      break;
    }
  }
  geom.unrequireEdgeCotanWeights();

  if (isDelaunay) {
    // TODO we said the matrix should be SPD if Delaunay, but SuiteSparse is failing with non-SPD error
    // (observed by nsharp on intrinsic delaunay RamSkull100k.obj)
    // workaround with a try-catch for now
    try {
      vectorHeatSolver.reset(new PositiveDefiniteSolver<std::complex<double>>(vectorOp));
    } catch (const std::runtime_error&) {
      vectorHeatSolver.reset(new SquareSolver<std::complex<double>>(vectorOp));
    }
  } else {
    vectorHeatSolver.reset(new SquareSolver<std::complex<double>>(vectorOp)); // not necessarily SPD without Delaunay
  }

  geom.unrequireVertexConnectionLaplacian();
}

void VectorHeatMethodSolver::ensureHaveAffineHeatSolver() {
  if (affineHeatSolver != nullptr) return;

  geom.requireVertexIndices();
  geom.requireEdgeCotanWeights();
  geom.requireTransportVectorsAlongHalfedge();
  geom.requireHalfedgeVectorsInVertex();
  geom.requireVertexDualAreas();

  // Get the ingredients
  std::vector<Eigen::Triplet<double>> triplets;
  for (Halfedge he : mesh.halfedges()) {
    size_t iTail = geom.vertexIndices[he.vertex()];
    size_t iTip = geom.vertexIndices[he.twin().vertex()];

    double weight = geom.edgeCotanWeights[he.edge()];

    Vector2 eij = geom.halfedgeVectorsInVertex[he];
    Vector2 eji = geom.halfedgeVectorsInVertex[he.twin()];

    Vector2 rot = -eij / eji;
    Vector2 t = -eij;

    // Rotation from TjM -> TiM, followed by a translation TjM -> TjM
    // The order matters here!
    DenseMatrix<double> transport(3, 3);
    transport << rot.x, -rot.y, t.x, rot.y, rot.x, t.y, 0, 0, 1.;

    for (size_t p = 0; p < 3; p++) {
      triplets.emplace_back(3 * iTail + p, 3 * iTail + p, weight);
      for (size_t q = 0; q < 3; q++) {
        triplets.emplace_back(3 * iTail + p, 3 * iTip + q, -weight * transport(p, q));
      }
    }
  }

  // Build the affine connection Laplacian
  SparseMatrix<double> L_affineConnection = Eigen::SparseMatrix<double>(3 * mesh.nVertices(), 3 * mesh.nVertices());
  L_affineConnection.setFromTriplets(triplets.begin(), triplets.end());

  // Build the lumped mass matrix (with 3x identical entries for each 3x3 block)
  Vector<double> massVec(3 * mesh.nVertices());
  for (Vertex v : mesh.vertices()) {
    for (size_t p = 0; p < 3; p++) {
      massVec(3 * geom.vertexIndices[v] + p) = geom.vertexDualAreas[v];
    }
  }
  SparseMatrix<double> massMat(massVec.asDiagonal());

  // Build (and factor) the short-time heat flow operator
  SparseMatrix<double> L_affineHeatOp = massMat + shortTime * L_affineConnection;
  affineHeatSolver.reset(
      new SquareSolver<double>(L_affineHeatOp)); // not necessarily SPD in the scalar sense, even with Delaunay

  // Bookkeeping
  geom.unrequireVertexIndices();
  geom.unrequireEdgeCotanWeights();
  geom.unrequireTransportVectorsAlongHalfedge();
  geom.unrequireHalfedgeVectorsInVertex();
  geom.unrequireVertexDualAreas();
}

void VectorHeatMethodSolver::ensureHavePoissonSolver() {
  if (poissonSolver != nullptr) return;

  // Get the ingredients
  geom.requireCotanLaplacian();
  SparseMatrix<double>& L = geom.cotanLaplacian;

  // Build the operator
  poissonSolver.reset(new PositiveDefiniteSolver<double>(L));

  geom.unrequireCotanLaplacian();
}

VertexData<double> VectorHeatMethodSolver::extendScalar(const std::vector<std::tuple<Vertex, double>>& sources) {

  std::vector<std::tuple<SurfacePoint, double>> sourcePoints;
  for (auto tup : sources) {
    sourcePoints.emplace_back(SurfacePoint(std::get<0>(tup)), std::get<1>(tup));
  }

  // call general version
  return extendScalar(sourcePoints);
}

VertexData<double> VectorHeatMethodSolver::extendScalar(const std::vector<std::tuple<SurfacePoint, double>>& sources) {
  if (sources.size() == 0) {
    return VertexData<double>(mesh, std::numeric_limits<double>::quiet_NaN());
  }

  ensureHaveScalarHeatSolver();
  geom.requireVertexIndices();

  // === Build the RHS
  Vector<double> dataRHS = Vector<double>::Zero(mesh.nVertices());
  Vector<double> indicatorRHS = Vector<double>::Zero(mesh.nVertices());

  for (auto tup : sources) {
    SurfacePoint point = std::get<0>(tup);
    double value = std::get<1>(tup);

    SurfacePoint facePoint = point.inSomeFace();
    Halfedge he = facePoint.face.halfedge();

    { // First adjacent vertex
      size_t vInd = geom.vertexIndices[he.vertex()];
      double w = facePoint.faceCoords.x;
      dataRHS[vInd] += w * value;
      indicatorRHS[vInd] += w;
    }
    he = he.next();

    { // Second adjacent vertex
      size_t vInd = geom.vertexIndices[he.vertex()];
      double w = facePoint.faceCoords.y;
      dataRHS[vInd] += w * value;
      indicatorRHS[vInd] += w;
    }
    he = he.next();

    { // Third adjacent vertex
      size_t vInd = geom.vertexIndices[he.vertex()];
      double w = facePoint.faceCoords.z;
      dataRHS[vInd] += w * value;
      indicatorRHS[vInd] += w;
    }
  }


  // == Solve the systems
  Vector<double> dataSol = scalarHeatSolver->solve(dataRHS);
  Vector<double> indicatorSol = scalarHeatSolver->solve(indicatorRHS);


  // == Combine results
  Vector<double> interpResult = dataSol.array() / indicatorSol.array();
  VertexData<double> result(mesh, interpResult);

  geom.unrequireVertexIndices();

  return result;
}

VertexData<Vector2> VectorHeatMethodSolver::transportTangentVector(Vertex sourceVert, Vector2 sourceVec) {
  std::vector<std::tuple<Vertex, Vector2>> sources{std::tuple<Vertex, Vector2>{sourceVert, sourceVec}};
  return transportTangentVectors(sources);
}

VertexData<Vector2>
VectorHeatMethodSolver::transportTangentVectors(const std::vector<std::tuple<Vertex, Vector2>>& sources) {
  std::vector<std::tuple<SurfacePoint, Vector2>> sourcesSurf;
  for (const auto& tup : sources) {
    sourcesSurf.emplace_back(SurfacePoint(std::get<0>(tup)), std::get<1>(tup));
  }
  return transportTangentVectors(sourcesSurf);
}

VertexData<Vector2>
VectorHeatMethodSolver::transportTangentVectors(const std::vector<std::tuple<SurfacePoint, Vector2>>& sources) {
  if (sources.size() == 0) {
    return VertexData<Vector2>(mesh, Vector2::undefined());
  }

  geom.requireVertexIndices();


  // === Setup work

  // Don't need to do magnitude solve with a single source
  bool singleVec = sources.size() == 1;

  // Make sure systems have been built and factored
  ensureHaveVectorHeatSolver();
  if (!singleVec) {
    ensureHaveScalarHeatSolver();
  }


  // === Build the RHS

  Vector<std::complex<double>> dirRHS = Vector<std::complex<double>>::Zero(mesh.nVertices());

  // Accumulate magnitude data for scalar problem
  std::vector<std::tuple<SurfacePoint, double>> magnitudeSources;

  for (auto tup : sources) {
    SurfacePoint point = std::get<0>(tup);
    Vector2 vec = std::get<1>(tup);
    std::complex<double> unitVec = Vector2::fromComplex(vec).normalize();

    // Add to the list of magnitudes for magnitude interpolation
    magnitudeSources.emplace_back(point, vec.norm());

    SurfacePoint facePoint = point.inSomeFace();
    Halfedge he = facePoint.face.halfedge();

    { // First adjacent vertex
      size_t vInd = geom.vertexIndices[he.vertex()];
      double w = facePoint.faceCoords.x;
      dirRHS[vInd] += w * unitVec;
    }
    he = he.next();


    { // Second adjacent vertex
      size_t vInd = geom.vertexIndices[he.vertex()];
      double w = facePoint.faceCoords.y;
      dirRHS[vInd] += w * unitVec;
    }
    he = he.next();


    { // Third adjacent vertex
      size_t vInd = geom.vertexIndices[he.vertex()];
      double w = facePoint.faceCoords.z;
      dirRHS[vInd] += w * unitVec;
    }
  }


  // == Solve the system

  Vector<std::complex<double>> vecSolution = vectorHeatSolver->solve(dirRHS);


  // == Get the magnitude right

  VertexData<Vector2> result(mesh);
  if (singleVec) {
    // For one sources, can just normalize and project
    double targetNorm = std::get<1>(sources[0]).norm();

    vecSolution = (vecSolution.array() / vecSolution.array().abs()) * targetNorm;

    // Copy to output vector
    for (Vertex v : mesh.vertices()) {
      result[v] = Vector2::fromComplex(vecSolution[geom.vertexIndices[v]]);
    }
  } else {
    // For multiple sources, need to interpolate magnitudes

    // === Perform scalar interpolation
    VertexData<double> interpMags = extendScalar(magnitudeSources);

    // Scale and copy to result
    for (Vertex v : mesh.vertices()) {
      Vector2 dir = Vector2::fromComplex(vecSolution[geom.vertexIndices[v]]).normalize();
      result[v] = dir * interpMags[v];
    }
  }


  geom.unrequireVertexIndices();
  return result;
}

VertexData<Vector2> VectorHeatMethodSolver::computeLogMap(const Vertex& sourceVert, LogMapStrategy strategy) {
  switch (strategy) {
  case LogMapStrategy::VectorHeat:
    return computeLogMap_VectorHeat(sourceVert, 0.0);
  case LogMapStrategy::AffineLocal:
    return computeLogMap_AffineLocal(sourceVert);
  case LogMapStrategy::AffineAdaptive:
    return computeLogMap_AffineAdaptive(sourceVert);
  }
}

VertexData<Vector2> VectorHeatMethodSolver::computeLogMap_AffineLocal(const Vertex& sourceVert) {
  ensureHaveVectorHeatSolver();
  ensureHaveAffineHeatSolver();

  geom.requireVertexIndices();

  // Construct the parallel section
  VertexData<Vector2> parallel_section = transportTangentVector(sourceVert, Vector2{1, 0});

  // Build rhs
  Vector<double> horizontalRHS = Vector<double>::Zero(3 * mesh.nVertices());
  horizontalRHS[3 * geom.vertexIndices[sourceVert] + 2] += 1.0; // (0,0), in homogeneous coordinates

  // Short-time heat flow on of hte affine values
  Vector<double> horizontalSol = affineHeatSolver->solve(horizontalRHS);

  // Construct the map from the heat-flow'd values
  VertexData<Vector2> logmap(mesh);
  for (Vertex v : mesh.vertices()) {
    size_t vInd = geom.vertexIndices[v];

    Vector2 logCoord = Vector2{horizontalSol(3 * vInd + 0), horizontalSol(3 * vInd + 1)};
    logCoord /= horizontalSol(3 * vInd + 2); // homogenous division

    // Change of basis to coordinate frame at source
    logmap[v] = logCoord / parallel_section[v];
  }

  // Bookkeeping
  geom.unrequireVertexIndices();

  return logmap;
}

VertexData<Vector2> VectorHeatMethodSolver::computeLogMap_AffineAdaptive(const Vertex& sourceVert) {

  geom.requireVertexIndices();
  geom.requireFaceAreas();
  geom.requireVertexDualAreas();
  geom.requireEdgeLengths();
  geom.requireEdgeCotanWeights();
  geom.requireHalfedgeVectorsInVertex();

  // Construct a parallel section using the coordinate from the source
  VertexData<Vector2> parallel_section = transportTangentVector(sourceVert, Vector2{1., 0.});

  // Build the adaptive connection Laplacian in the coordinate frame
  // (this depends on the source, so it cannot be prefactored)
  std::vector<Eigen::Triplet<double>> triplets;
  for (Halfedge he : mesh.halfedges()) {
    size_t iTail = geom.vertexIndices[he.vertex()];
    size_t iTip = geom.vertexIndices[he.twin().vertex()];

    double weight = geom.edgeCotanWeights[he.edge()];

    // Translation in the edge direction in the trivialization induced by the geodesic frame X
    // We have two options on how to do it. Turns out it doesn't really matter which we pick.

    // Do the translation using the values at i
    // (we could equivalently do it at j)
    Vector2 eij = geom.halfedgeVectorsInVertex[he];
    Vector2 Xi = parallel_section[iTail];
    Vector2 t = -eij / Xi;

    // The developing connection
    DenseMatrix<double> transport(3, 3);
    transport << 1, 0, t.x, 0, 1, t.y, 0, 0, 1.;

    for (size_t p = 0; p < 3; p++) {
      triplets.emplace_back(3 * iTail + p, 3 * iTail + p, weight);
      for (size_t q = 0; q < 3; q++) {
        triplets.emplace_back(3 * iTail + p, 3 * iTip + q, -weight * transport(p, q));
      }
    }
  }
  SparseMatrix<double> L_affineAdaptive = Eigen::SparseMatrix<double>(3 * mesh.nVertices(), 3 * mesh.nVertices());
  L_affineAdaptive.setFromTriplets(triplets.begin(), triplets.end());

  // Build the corresponding (lumped) mass matrix, with 3x3 vertex area diagonal blocks
  Vector<double> massVec(3 * mesh.nVertices());
  for (Vertex v : mesh.vertices())
    for (size_t p = 0; p < 3; p++) massVec(3 * geom.vertexIndices[v] + p) = geom.vertexDualAreas[v];
  SparseMatrix<double> massMat(massVec.asDiagonal());

  // Form the short-time heat flow operator and solver
  SparseMatrix<double> affineOp = massMat + shortTime * L_affineAdaptive;
  SquareSolver<double> affineAdaptiveSolver(affineOp);

  // Build rhs, with a homogenous 0 at the source
  Vector<double> horizontalRHS = Vector<double>::Zero(3 * mesh.nVertices());
  horizontalRHS[3 * geom.vertexIndices[sourceVert] + 2] += 1.0;

  // Solve
  Vector<double> horizontalSol = affineAdaptiveSolver.solve(horizontalRHS);

  // Compose the result
  VertexData<Vector2> logmap(mesh);
  for (Vertex v : mesh.vertices()) {
    size_t vInd = geom.vertexIndices[v];

    Vector2 logCoord = Vector2{horizontalSol(3 * vInd + 0), horizontalSol(3 * vInd + 1)};
    logCoord /= horizontalSol(3 * vInd + 2);

    logmap[v] = logCoord;
  }

  // Bookkeeping
  geom.unrequireVertexIndices();
  geom.unrequireFaceAreas();
  geom.unrequireVertexDualAreas();
  geom.unrequireEdgeLengths();
  geom.unrequireEdgeCotanWeights();
  geom.unrequireHalfedgeVectorsInVertex();

  return logmap;
}

VertexData<Vector2> VectorHeatMethodSolver::computeLogMap_VectorHeat(const Vertex& sourceVert,
                                                                     double vertexDistanceShift) {
  geom.requireFaceAreas();
  geom.requireEdgeLengths();
  geom.requireCornerAngles();
  geom.requireEdgeCotanWeights();
  geom.requireHalfedgeVectorsInVertex();
  geom.requireTransportVectorsAlongHalfedge();
  geom.requireVertexIndices();


  // Make sure systems have been built and factored
  ensureHaveVectorHeatSolver();
  ensureHavePoissonSolver();

  // === Solve for "radial" field

  // Build rhs
  Vector<std::complex<double>> radialRHS = Vector<std::complex<double>>::Zero(mesh.nVertices());
  addVertexOutwardBall(sourceVert, radialRHS);

  // Solve
  Vector<std::complex<double>> radialSol = vectorHeatSolver->solve(radialRHS);

  // Normalize
  radialSol = (radialSol.array() / radialSol.array().abs());
  radialSol[geom.vertexIndices[sourceVert]] = 0.;


  // === Solve for "horizontal" field

  // Build rhs
  Vector<std::complex<double>> horizontalRHS = Vector<std::complex<double>>::Zero(mesh.nVertices());
  horizontalRHS[geom.vertexIndices[sourceVert]] += 1.0;

  // Solve
  Vector<std::complex<double>> horizontalSol = vectorHeatSolver->solve(horizontalRHS);

  // Normalize
  horizontalSol = (horizontalSol.array() / horizontalSol.array().abs());


  // === Integrate radial field to get distance

  // Build the right hand sign (divergence term)
  Vector<double> divergenceVec = Vector<double>::Zero(mesh.nVertices());
  for (Halfedge he : mesh.halfedges()) {

    // Build the vector which is the average vector along the edge, in the basis of the tail vertex
    Vector2 radAtTail = Vector2::fromComplex(radialSol[geom.vertexIndices[he.vertex()]]);
    Vector2 radAtTip = Vector2::fromComplex(radialSol[geom.vertexIndices[he.twin().vertex()]]);
    Vector2 radTipAtTail = geom.transportVectorsAlongHalfedge[he.twin()] * radAtTip;

    // Integrate the edge vector along the edge
    Vector2 vectAtEdge = 0.5 * (radAtTail + radTipAtTail);
    double fieldAlongEdge = dot(vectAtEdge, geom.halfedgeVectorsInVertex[he]);

    // Contrbution to divergence is cotan times that integral
    // (negative since we want negative divergence due to Laplacian sign)
    double weight = geom.edgeCotanWeights[he.edge()];
    divergenceVec[geom.vertexIndices[he.vertex()]] += -weight * fieldAlongEdge;
  }

  // Integrate to get distance
  Vector<double> distance = poissonSolver->solve(divergenceVec);

  // Shift distance to be zero at the source
  distance = distance.array() + (vertexDistanceShift - distance[geom.vertexIndices[sourceVert]]);


  // Combine distance and angle to get cartesian result
  VertexData<Vector2> result(mesh);
  for (Vertex v : mesh.vertices()) {
    size_t vInd = geom.vertexIndices[v];

    std::complex<double> logDir = radialSol[vInd] / horizontalSol[vInd];
    Vector2 logCoord = Vector2::fromComplex(logDir) * distance[vInd];
    result[v] = logCoord;
  }

  return result;
}

VertexData<double> VectorHeatMethodSolver::scalarDiffuse(const VertexData<double>& rhs) {
  ensureHaveScalarHeatSolver();
  return VertexData<double>(mesh, scalarHeatSolver->solve(rhs.toVector()));
}

VertexData<std::complex<double>> VectorHeatMethodSolver::vectorDiffuse(const VertexData<std::complex<double>>& rhs) {
  ensureHaveVectorHeatSolver();
  return VertexData<std::complex<double>>(mesh, vectorHeatSolver->solve(rhs.toVector()));
}

VertexData<double> VectorHeatMethodSolver::poissonSolve(const VertexData<double>& rhs) {
  ensureHavePoissonSolver();
  return VertexData<double>(mesh, poissonSolver->solve(rhs.toVector()));
}

void VectorHeatMethodSolver::addVertexOutwardBall(Vertex vert, Vector<std::complex<double>>& distGradRHS) {

  // see Vector Heat Method, Appendix A

  // Height of triangle with tip at he.vertex()
  auto heightInTriangle = [&](Halfedge he) {
    double area = geom.faceAreas[he.face()];
    double base = geom.edgeLengths[he.next().edge()];
    return 2.0 * area / base;
  };

  // Contribution to distance gradient right hand side

  size_t vInd = geom.vertexIndices[vert];
  for (Halfedge he : vert.outgoingHalfedges()) {
    Vertex vn = he.twin().vertex();
    size_t vnInd = geom.vertexIndices[vn];


    // he side
    if (he.isInterior()) {
      double h = heightInTriangle(he.next());
      // double theta = halfedgeOppositeAngles[he.next()];
      double theta = geom.cornerAngles[he.corner()];

      Vector2 valInEdgeBasis{-theta * std::sin(theta) / (2.0 * h),
                             (theta * std::cos(theta) - std::sin(theta)) / (2.0 * h)};


      distGradRHS[vnInd] +=
          static_cast<std::complex<double>>(valInEdgeBasis * geom.halfedgeVectorsInVertex[he.twin()].normalize());
    }

    // he.twin() side
    if (he.twin().isInterior()) {
      double h = heightInTriangle(he.twin());
      // double theta = halfedgeOppositeAngles[he.twin().next().next()];
      double theta = geom.cornerAngles[he.twin().next().corner()];

      Vector2 valInEdgeBasis{-theta * std::sin(theta) / (2.0 * h),
                             -(theta * std::cos(theta) - std::sin(theta)) / (2.0 * h)};

      distGradRHS[vnInd] +=
          static_cast<std::complex<double>>(valInEdgeBasis * geom.halfedgeVectorsInVertex[he.twin()].normalize());
    }

    // Contribution  to center vert
    if (he.isInterior()) {
      double h = heightInTriangle(he);
      double theta = geom.cornerAngles[he.corner()];
      double gamma = geom.cornerAngles[he.next().corner()];
      double alpha = M_PI / 2.0 - gamma;


      Vector2 valInEdgeBasis{-(theta * std::cos(alpha) + std::cos(alpha - theta) * std::sin(theta)) / (2.0 * h),
                             -(std::cos(alpha) - std::cos(alpha - 2 * theta) + 2 * theta * std::sin(alpha)) /
                                 (4.0 * h)};

      distGradRHS[vInd] +=
          static_cast<std::complex<double>>(valInEdgeBasis * geom.halfedgeVectorsInVertex[he].normalize());
    }
  }
}

VertexData<Vector2> VectorHeatMethodSolver::computeLogMap(const SurfacePoint& sourceP, LogMapStrategy strategy) {
  geom.requireHalfedgeVectorsInVertex();
  geom.requireHalfedgeVectorsInFace();

  switch (sourceP.type) {
  case SurfacePointType::Vertex: {

    return computeLogMap(sourceP.vertex, strategy);
    break;
  }
  case SurfacePointType::Edge: {
    geom.requireHalfedgeVectorsInVertex();

    // Compute logmaps at both adjacent vertices
    Halfedge he = sourceP.edge.halfedge();
    VertexData<Vector2> logmapTail = computeLogMap(he.vertex(), strategy);
    VertexData<Vector2> logmapTip = computeLogMap(he.twin().vertex(), strategy);

    // Changes of basis
    Vector2 tailRot = geom.halfedgeVectorsInVertex[he].inv().normalize();
    Vector2 tipRot = -geom.halfedgeVectorsInVertex[he.twin()].inv().normalize();

    // Blend result and store in edge basis
    VertexData<Vector2> resultMap(mesh, Vector2::zero());
    double tBlend = sourceP.tEdge;
    for (Vertex v : mesh.vertices()) {
      resultMap[v] = (1. - tBlend) * logmapTail[v] * tailRot + tBlend * logmapTip[v] * tipRot;
    }

    geom.unrequireHalfedgeVectorsInVertex();
    return resultMap;
    break;
  }
  case SurfacePointType::Face: {
    geom.requireHalfedgeVectorsInVertex();
    geom.requireHalfedgeVectorsInFace();

    // Accumulate result from adjcent halfedges
    VertexData<Vector2> resultMap(mesh, Vector2::zero());
    int iC = 0;
    for (Halfedge he : sourceP.face.adjacentHalfedges()) {

      // Commpute logmap at vertex
      VertexData<Vector2> logmapVert = computeLogMap(he.vertex(), strategy);

      // Compute change of basis to bring it back to the face
      Vector2 rot = (geom.halfedgeVectorsInFace[he] / geom.halfedgeVectorsInVertex[he]).normalize();

      // Accumulate in face fesult
      for (Vertex v : mesh.vertices()) {
        resultMap[v] += sourceP.faceCoords[iC] * rot * logmapVert[v];
      }
      iC++;
    }


    geom.unrequireHalfedgeVectorsInVertex();
    geom.unrequireHalfedgeVectorsInFace();
    return resultMap;
    break;
  }
  }

  throw std::logic_error("bad switch");
  return VertexData<Vector2>();
}

} // namespace surface
} // namespace geometrycentral
