#include "geometrycentral/surface/signed_heat_method.h"

namespace geometrycentral {
namespace surface {

SignedHeatMethodSolver::SignedHeatMethodSolver(IntrinsicGeometryInterface& geom_, double tCoef_)
    : tCoef(tCoef_), mesh(geom_.mesh), geom(geom_)

{
  geom.requireEdgeLengths();
  geom.requireCrouzeixRaviartMassMatrix();

  // Compute mean edge length and set shortTime
  double meanEdgeLength = 0.;
  for (Edge e : mesh.edges()) {
    meanEdgeLength += geom.edgeLengths[e];
  }
  meanEdgeLength /= mesh.nEdges();
  meanNodeDistance = 0.5 * meanEdgeLength;
  shortTime = tCoef * meanNodeDistance * meanNodeDistance;

  massMat = geom.crouzeixRaviartMassMatrix;
  doubleVectorOp = crouzeixRaviartDoubleMassMatrix + shortTime * crouzeixRaviartDoubleConnectionLaplacian();

  geom.unrequireCrouzeixRaviartMassMatrix();
  geom.unrequireEdgeLengths();
}

VertexData<double> SignedHeatMethodSolver::computeDistance(const std::vector<Curve>& curves,
                                                           const std::vector<SurfacePoint>& points,
                                                           const SignedHeatOptions& options) {

  Vector<std::complex<double>> Xt = integrateVectorHeatFlow(curves, points, options);
  return integrateVectorField(Xt, curves, points, options);
}

VertexData<double> SignedHeatMethodSolver::computeDistance(const std::vector<Curve>& curves,
                                                           const SignedHeatOptions& options) {
  // call generic version
  return computeDistance(curves, std::vector<SurfacePoint>(), options);
}

VertexData<double> SignedHeatMethodSolver::computeDistance(const std::vector<SurfacePoint>& points) {

  // call generic version
  return computeDistance(std::vector<Curve>(), points, options);
}

Vector<std::complex<double>> SignedHeatMethodSolver::integrateVectorHeatFlow(const std::vector<Curve>& curves,
                                                                             const std::vector<SurfacePoint>& points,
                                                                             const SignedHeatOptions& options) {

  ensureHaveVectorHeatSolver();

  geom.requireEdgeIndices();

  // Build source term.
  size_t E = mesh.nEdges();
  Vector<std::complex<double>> X0 = Vector<std::complex<double>>::Zero(E);
  for (const Curve& curve : curves) {
    if (curve.isSigned) {
      buildSignedCurveSource(curve, X0)
    } else {
      buildUnsignedCurveSource(curve, X0)
    }
  }
  geom.requireCornerAngles();
  for (const SurfacePoint& point : points) buildUnsignedPointSource(point, X0);
  geom.unrequireCornerAngles();

  Vector<std::complex<double>> Xt;
  if (!options.preserveSourceNormals) {
    Xt = vectorHeatSolver->solve(X0);
  } else {
    // Convert X0 to double vector
    Vector<double> rhs = Vector<double>::Zero(2 * E);
    for (size_t i = 0; i < E; i++) {
      rhs[i] = std::real(X0[i]);
      rhs[E + i] = std::imag(X0[i]);
    }

    // Build constraint matrix
    std::vector<Eigen::Triplet<double>> triplets;
    SparseMatrix<double> C;
    size_t m = 0;
    for (const auto& curve : curves) {
      size_t nNodes = curve.nodes.size();
      for (size_t i = 0; i < nNodes - 1; i++) {
        const SurfacePoint& pA = curve.nodes[i];
        const SurfacePoint& pB = curve.nodes[i + 1];
        SurfacePoint b = midSegmentSurfacePoint(pA, pB);
        Edge commonEdge = sharedEdge(pA, pB);
        if (commonEdge != Edge()) {
          size_t eIdx = geom.edgeIndices[commonEdge];
          // need parallel component to be zero
          triplets.emplace_back(m, eIdx, 1);
          m++;
          continue;
        } else {
          Face f = sharedFace(pA, pB);
          BarycentricVector segTangent(pA, pB);
          for (Edge e : f.adjacentEdges()) {
            size_t eIdx = geom.edgeIndices[e];
            // double w = scalarCrouzeixRaviart(b, e);
            double w = 1.; // TODO
            BarycentricVector heVec = barycentricVectorInFace(e.halfedge(), f);
            heVec /= heVec.norm(geom);
            BarycentricVector heVecN = heVec.rotate90(geom);
            triplets.emplace_back(m, eIdx, w * dot(geom, heVec, segTangent));
            triplets.emplace_back(m, E + eIdx, w * dot(geom, heVecN, segTangent));
          }
          m++;
        }
      }
    }
    C.resize(m, 2 * E);
    C.setFromTriplets(triplets.begin(), triplets.end());

    SparseMatrix<double> Z(m, m);
    SparseMatrix<double> LHS1 = horizontalStack<double>({doubleVectorOp, C.transpose()});
    SparseMatrix<double> LHS2 = horizontalStack<double>({C, Z});
    SparseMatrix<double> LHS = verticalStack<double>({LHS1, LHS2});
    Vector<double> RHS = Vector<double>::Zero(2 * E + m);
    RHS.head(2 * E) = rhs;
    Vector<double> soln = solveSquare(LHS, RHS);
    Vector<double> X = soln.head(2 * E);

    // Convert back to complex vector.
    Xt = Vector<std::complex<double>>::Zero(E);
    for (size_t i = 0; i < E; i++) Xt[i] = std::complex<double>(X[i], X[E + i]);
  }

  geom.unrequireEdgeIndices();

  return Xt;
}

VertexData<double> SignedHeatMethodSolver::integrateVectorField(const Vector<std::complex<double>>& Xt,
                                                                const std::vector<Curve>& curves,
                                                                const std::vector<SurfacePoint>& points,
                                                                const SignedHeatOptions& options) {
  geom.requireHalfedgeCotanWeights();

  size_t V = mesh.nVertices();
  FaceData<BarycentricVector> Y = sampleAtFaceBarycenters(Xt); // sample and normalizes
  Vector<double> div = Vector<double>::Zero(V);
  for (Vertex v : mesh.vertices()) {
    size_t vIdx = geom.vertexIndices[v];
    for (Corner c : v.adjacentCorners()) {
      Face f = c.face();
      Halfedge heA = c.halfedge();
      Halfedge heB = heA.next().next();
      BarycentricVector Yj = Y[f];
      BarycentricVector e1 = barycentricVectorInFace(heA, f);
      BarycentricVector e2 = -barycentricVectorInFace(heB, f);
      double cotTheta1 = geom.halfedgeCotanWeights[heA];
      double cotTheta2 = geom.halfedgeCotanWeights[heB];
      double w1 = cotTheta1 * dot(geom, e1, Yj);
      double w2 = cotTheta2 * dot(geom, e2, Yj);
      div[vIdx] += w1;
      div[vIdx] += w2;
    }
  }

  Vector<double> phi;
  switch (options.levelSetConstraint) {
  case (LevelSetConstraint::None): {
    phi = poissonSolver->solve(div);
    break;
  }
  case (LevelSetConstraint::ZeroSet): {
    phi = integrateWithZeroSetConstraint(div, curves, options);
    break;
  }
  case (LevelSetConstraint::Multiple): {
    phi = integrateWithLevelSetConstraints(div, curves, options);
    break;
  }
  }

  geom.unrequireHalfedgeCotanWeights();

  double shift = computeAverageValueOnSource(phi, curves);
  phi -= Vector<double>::Ones(V) * shift;
  return -phi; // since our Laplacian is positive
}

void SignedHeatMethodSolver::buildSignedCurveSource(const Curve& curve, Vector<std::complex<double>>& X0) const {

  size_t nNodes = curve.nodes.size();
  for (size_t i = 0; i < nNodes - 1; i++) {
    const SurfacePoint& pA = curve.nodes[i];
    const SurfacePoint& pB = curve.nodes[i + 1];
    Edge commonEdge = sharedEdge(pA, pB);
    if (commonEdge != Edge()) {
      size_t eIdx = geom.edgeIndices[commonEdge];
      double scalar = 1.;
      std::complex<double> innerProd(0, 1);
      if (pA.vertex == commonEdge.secondVertex()) innerProd.imag(-1);
      double length = lengthOfSegment(pA, pB);
      std::complex<double> contrib = length * scalar * innerProd;
      X0[eIdx] += contrib;
      continue;
    }
    Face commonFace = sharedFace(pA, pB);
    for (Edge e : commonFace.adjacentEdges()) {
      size_t eIdx = geom.edgeIndices[e];
      std::complex<double> innerProd = projectedNormal(pA, pB, e); // already incorporates length
      X0[eIdx] += innerProd;
    }
  }
}

void SignedHeatMethodSolver::buildUnsignedPointSource(const SurfacePoint& point,
                                                      Vector<std::complex<double>>& X0) const {

  switch (point.p.type) {
  case (SurfacePointType::Vertex): {
    Vertex v = point.p.vertex;
    buildUnsignedVertexSource(v, X0);
    break;
  }
  case (SurfacePointType::Edge): {
    Edge e = point.p.edge;
    double t = point.p.tEdge;
    // TODO
    buildUnsignedVertexSource(e.firstVertex(), X0, 1. - t);
    buildUnsignedVertexSource(e.secondVertex(), X0, t);
    break;
  }
  case (SurfacePointType::Face): {
    Face f = point.p.face;
    int idx = 0;
    for (Vertex v : f.adjacentVertices()) {
      buildUnsignedVertexSource(v, X0, p.faceCoords[idx]);
      idx++;
    }
    break;
  }
  default: {
    throw std::logic_error("buildUnsignedPointSource(): bad switch");
  }
  }
}

void SignedHeatMethodSolver::buildUnsignedVertexSource(const Vertex& v, Vector<std::complex<double>>& X0,
                                                       double weight) const {

  double angleSum = 0.;
  for (Corner c : v.adjacentCorners()) angleSum += geom.cornerAngles[c];
  // Note: On non-orientable or non-manifold meshes, corner/face iterator may not be reliable.
  for (Corner c : v.adjacentCorners()) {
    Face f = c.face();
    Halfedge ij = c.halfedge();
    Halfedge jk = ij.next();
    Halfedge ki = jk.next();
    double theta_i = geom.cornerAngles[c];
    double theta_j = geom.cornerAngles[jk.corner()];
    Edge e_ij = ij.edge();
    Edge e_jk = jk.edge();
    Edge e_ki = ki.edge();
    size_t eIdx_ij = geom.edgeIndices[e_ij];
    size_t eIdx_jk = geom.edgeIndices[e_jk];
    size_t eIdx_ki = geom.edgeIndices[e_ki];
    std::complex<double> exp_i = Vector2::fromAngle(theta_i);
    std::complex<double> im(0., 1.);
    double s_ij = ij.orientation() ? 1. : -1.;
    double s_jk = jk.orientation() ? 1. : -1.;
    double s_ki = ki.orientation() ? 1. : -1.;
    std::complex<double> r_ki_ij = -Vector2::fromAngle(-theta_i) * s_ki;
    std::complex<double> r_jk_ij = -Vector2::fromAngle(theta_j) * s_jk;

    std::complex<double> w_ij = im * (1. - exp_i) / angleSum;

    X0[eIdx_ij] += s_ij * w_ij;
    X0[eIdx_jk] += r_jk_ij * w_ij;
    X0[eIdx_ki] += r_ki_ij * w_ij;
  }
}

Vector<double> SignedHeatMethodSolver::integrateWithZeroSetConstraint(const Vector<double>& rhs,
                                                                      const std::vector<Curve>& curves,
                                                                      const SignedHeatOptions& options) {

  geom.requireCotanLaplacian();
  geom.requireVertexIndices();

  // If curve is constrained to mesh edges, solve the system by omitting DOFs.
  bool allVertices = true;
  for (const Curve& curve : curves) {
    for (const SurfacePoint& p : curve.nodes) {
      if (p.type != SurfacePointType::Vertex) {
        allVertices = false;
        continue;
      }
    }
  }

  size_t V = mesh.nVertices();
  if (allVertices && options.softLevelSetWeight < 0.) {
    Vector<bool> setAMembership = Vector<bool>::Ones(V); // true if interior
    for (const Curve& curve : curves) {
      for (const SurfacePoint& p : curve.nodes) {
        Vertex v = p.vertex;
        setAMembership[geom.vertexIndices[v]] = false;
      }
    }
    int nB = V - setAMembership.cast<int>().sum(); // # of boundary vertices
    Vector<double> bcVals = Vector<double>::Zero(nB);
    BlockDecompositionResult<double> decomp = blockDecomposeSquare(geom.cotanLaplacian, setAMembership, true);
    Vector<double> rhsValsA, rhsValsB;
    decomposeVector(decomp, rhs, rhsValsA, rhsValsB);
    Vector<double> combinedRHS = rhsValsA;
    PositiveDefiniteSolver constrainedSolver(decomp.AA);
    Vector<double> Aresult = constrainedSolver.solve(combinedRHS);
    phi = reassembleVector(decomp, Aresult, bcVals);
  } else {
    Curve allCurves;
    for (const Curve& curve : curves) {
      allCurves.nodes.insert(allCurves.nodes.end(), curve.nodes.begin(), curve.nodes.end());
    }
    phi = integrateWithLevelSetConstraints(rhs, {allCurves}, options);
  }

  geom.unrequireCotanLaplacian();
  geom.unrequireVertexIndices();

  return phi;
}

Vector<double> SignedHeatMethodSolver::integrateWithLevelSetConstraints(const Vector<double>& rhs,
                                                                        const std::vector<Curve>& curves,
                                                                        const SignedHeatOptions& options) {

  geom.requireCotanLaplacian();
  geom.requireVertexIndices();

  // Build constraint matrix.
  size_t V = mesh.nVertices();
  std::vector<Eigen::Triplet<double>> triplets;
  SparseMatrix<double> A; // constraint matrix
  size_t m = 0;
  for (const auto& curve : curves) {
    size_t nNodes = curve.nodes.size();
    // Technically we should check the entire curve for duplicate nodes.
    // Right now we assume the only possibly duplicated nodes are at curve endpoints.
    bool isClosed = (curve.nodes[0] == curve.nodes[nNodes - 1]);
    size_t uB = isClosed ? nNodes - 1 : nNodes;
    // Compute first node.
    const SurfacePoint& p0 = curve.nodes[0];
    std::vector<Eigen::Triplet<double>> p0_data;
    if (p0.type == SurfacePointType::Vertex) {
      size_t vIdx = geom.vertexIndices[p0.vertex];
      p0_data.emplace_back(m, vIdx, 1);
    } else if (p0.type == SurfacePointType::Edge) {
      size_t vA = geom.vertexIndices[p0.edge.firstVertex()];
      size_t vB = geom.vertexIndices[p0.edge.secondVertex()];
      p0_data.emplace_back(m, vA, 1. - p0.tEdge);
      p0_data.emplace_back(m, vB, p0.tEdge);
    }
    // Add constraints for the rest of the nodes.
    for (size_t i = 1; i < uB; i++) {
      const SurfacePoint& p = curve.nodes[i];
      if (p.type == SurfacePointType::Vertex) {
        size_t vIdx = geom.vertexIndices[p.vertex];
        triplets.emplace_back(m, vIdx, -1);
      } else if (p.type == SurfacePointType::Edge) {
        size_t vA = geom.vertexIndices[p.edge.firstVertex()];
        size_t vB = geom.vertexIndices[p.edge.secondVertex()];
        triplets.emplace_back(m, vA, -(1. - p.tEdge));
        triplets.emplace_back(m, vB, -p.tEdge);
      }
      for (auto& T : p0_data) {
        triplets.emplace_back(m, T.col(), T.value());
      }
      m++;
    }
  }
  A.resize(m, V);
  A.setFromTriplets(triplets.begin(), triplets.end());

  // Assemble and solve system.
  Vector<double> phi;
  if (options.softLevelSetWeight < 0.) {
    SparseMatrix<double> Z(m, m);
    SparseMatrix<double> LHS1 = horizontalStack<double>({geom.cotanLaplacian, A.transpose()});
    SparseMatrix<double> LHS2 = horizontalStack<double>({A, Z});
    SparseMatrix<double> LHS = verticalStack<double>({LHS1, LHS2});
    Vector<double> RHS = Vector<double>::Zero(V + m);
    RHS.head(V) = rhs;
    Vector<double> soln = solveSquare(LHS, RHS);
    phi = soln.head(V);
  } else {
    SparseMatrix<double> LHS = geom.cotanLaplacian + options.softLevelSetWeight * A.transpose() * A;
    phi = solvePositiveDefinite(LHS, rhs);
  }

  geom.unrequireCotanLaplacian();
  geom.unrequireVertexIndices();

  return phi;
}

double SignedHeatMethodSolver::lengthOfSegment(const SurfacePoint& pA, const SurfacePoint& pB) const {

  BarycentricVector segment(pA, pB);
  return segment.norm(geom);
}

SurfacePoint SignedHeatMethodSolver::midSegmentSurfacePoint(const SurfacePoint& pA, const SurfacePoint& pB) const {

  Face commonFace = sharedFace(pA, pB);
  SurfacePoint pAInFace = pA.inFace(commonFace);
  SurfacePoint pBInFace = pB.inFace(commonFace);
  SurfacePoint midpoint(commonFace, 0.5 * (pAInFace.faceCoords + pBInFace.faceCoords));
  return midpoint;
}

/*
 * Compute the normal to the segment as a complex number expressed in the local basis of the given edge. The magnitude
 * encodes the length of the segment.
 *
 * It is assumed that the curve segment lies within a face, not along an edge -- so there is a unique segment.face.
 */
std::complex<double> SignedHeatMethodSolver::projectedNormal(const SurfacePoint& pA, const SurfacePoint& pB,
                                                             const Edge& e) const {
  BarycentricVector segment(pA, pB);
  BarycentricVector segNormal = segment.rotate90(geom);

  Face f = segment.face;
  Vector3 edgeVecCoords = {0, 0, 0};
  int vIdx = 0;
  bool orientation; // relative orientation of e in face f
  for (Halfedge he : f.adjacentHalfedges()) {
    if (he.edge() == e) {
      edgeVecCoords[vIdx] = -1;
      edgeVecCoords[(vIdx + 1) % 3] = 1;
      orientation = (he.vertex() == e.firstVertex());
      break;
    }
    vIdx++;
  }
  if (!orientation) edgeVecCoords *= -1;
  BarycentricVector edgeVec(f, edgeVecCoords);
  edgeVec /= geom.edgeLengths[e];

  // <tangent, edgeVec> = sin(theta)
  // <J tangent, edgeVec> = cos(theta)
  // where normal := J tangent, theta = angle of normal with real axis.
  double sinTheta = dot(geom, segment, edgeVec);
  double cosTheta = dot(geom, segNormal, edgeVec);
  std::complex<double> normal = std::complex<double>(cosTheta, sinTheta);
  return normal;
}

BarycentricVector SignedHeatMethodSolver::barycentricVectorInFace(const Halfedge& he_, const Face& f) const {

  int eIdx = 0;
  double sign = 0.;
  for (Halfedge he : f.adjacentHalfedges()) {
    if (he.edge() == he_.edge()) {
      sign = (he.tailVertex() == he_.tailVertex() && he.tipVertex() == he_.tipVertex()) ? 1. : -1.;
      break;
    }
    eIdx++;
  }
  Vector3 faceCoords = {0, 0, 0};
  faceCoords[mod(eIdx + 1, 3)] = 1;
  faceCoords[eIdx] = -1;
  return BarycentricVector(f, sign * faceCoords);
}

FaceData<BarycentricVector> SignedHeatMethodSolver::sampleAtFaceBarycenters(const Vector<std::complex<double>>& Xt) {

  geom.requireEdgeIndices();

  FaceData<BarycentricVector> X(mesh);
  for (Face f : mesh.faces()) {
    Vector3 faceCoords = {0, 0, 0};
    for (Halfedge he : f.adjacentHalfedges()) {
      size_t eIdx = geom.edgeIndices[he.edge()];
      BarycentricVector e1 = barycentricVectorInFace(he, f);
      if (!he.orientation()) e1 *= -1;
      BarycentricVector e2 = e1.rotate90(geom);
      e1 /= e1.norm(geom);
      e2 /= e2.norm(geom);
      faceCoords += std::real(Xt[eIdx]) * e1.faceCoords;
      faceCoords += std::imag(Xt[eIdx]) * e2.faceCoords;
    }
    BarycentricVector vec(f, faceCoords);
    X[f] = vec / vec.norm(geom); // normalize
  }

  geom.unrequireEdgeIndices();

  return X;
}

double SignedHeatMethodSolver::computeAverageValueOnSource(const Vector<double>& phi,
                                                           const std::vector<Curve>& curves) const {

  geom.requireVertexIndices();

  // Compute the average value of phi along the input (integrate phi along the input geometry.)
  double shift = 0.;
  double normalization = 0.;
  for (const auto& curve : curves) {
    size_t nNodes = curve.nodes.size();
    for (size_t i = 0; i < nNodes - 1; i++) {
      const SurfacePoint& pA = curve.nodes[i];
      const SurfacePoint& pB = curve.nodes[i + 1];
      double length = lengthOfSegment(pA, pB);
      // Linearly interpolate distance from vertices to center of segment; integrate over segment using
      // one-point (midpoint) quadrature, which is exact for linear functions.
      SurfacePoint midpoint = midSegmentSurfacePoint(pA, pB);
      Face f = midpoint.face;
      int idx = 0;
      for (Vertex v : f.adjacentVertices()) {
        double w = midpoint.faceCoords[idx];
        size_t vIdx = geom.vertexIndices[v];
        shift += w * length * phi[vIdx];
        idx += 1;
      }
      normalization += length;
    }
  }
  shift /= normalization;

  geom.unrequireVertexIndices();

  return shift;
}

void SignedHeatMethodSolver::ensureHaveVectorHeatSolver() {

  if (vectorHeatSolver != nullptr) return;

  geom.requireCrouzeixRaviartConnectionLaplacian();

  SparseMatrix<std::complex<double>>& Lconn = geom.crouzeixRaviartConnectionLaplacian;
  SparseMatrix<std::complex<double>> vectorOp = massMat.cast<std::complex<double>>() + shortTime * Lconn;

  // Check the Delaunay condition. If the mesh is Delaunay, then vectorOp is SPD, and we can use a
  // PositiveDefiniteSolver. Otherwise, we must use a SquareSolver
  geom.requireEdgeCotanWeights();
  bool isDelaunay = true;
  for (Edge e : mesh.edges()) {
    if (geom.edgeCotanWeights[e] < -1e-6) {
      isDelaunay = false;
      break;
    }
  }
  geom.unrequireEdgeCotanWeights();

  if (isDelaunay) {
    vectorHeatSolver.reset(new PositiveDefiniteSolver<std::complex<double>>(vectorOp));
  } else {
    vectorHeatSolver.reset(new SquareSolver<std::complex<double>>(vectorOp));
  }

  geom.unrequireVertexConnectionLaplacian();
}

void SignedHeatMethodSolver::ensureHavePoissonSolver() {

  if (poissonSolver != nullptr) return;

  geom.requireCotanLaplacian();

  SparseMatrix<double>& L = geom.cotanLaplacian;
  poissonSolver.reset(new PositiveDefiniteSolver<double>(L));

  geom.unrequireCotanLaplacian();
}

/*
 * Build the connection Laplacian derived in Stein et al. 2020.
 * First are all the parallel elements, then all the perpendicular ones.
 */
SparseMatrix<double> SignedHeatMethodSolver::crouzeixRaviartDoubleConnectionLaplacian() const {

  geom.requireEdgeIndices();
  geom.requireEdgeLengths();
  geom.requireHalfedgeCotanWeights();

  size_t E = mesh.nEdges();
  SparseMatrix<double> L(2 * E, 2 * E);
  std::vector<Eigen::Triplet<double>> tripletList;
  for (Face f : mesh.faces()) {
    for (Halfedge he : f.adjacentHalfedges()) {
      Halfedge heA = he.next();
      Halfedge heB = heA.next();

      Corner c = heB.corner();
      size_t iE_i = geom.edgeIndices[heA.edge()];
      size_t iE_j = geom.edgeIndices[heB.edge()];

      double weight = 4.0 * geom.halfedgeCotanWeights[he];
      double s_ij = (heA.orientation() == heB.orientation()) ? 1. : -1.;
      double lOpp = geom.edgeLengths[he.edge()];
      double lA = geom.edgeLengths[heB.edge()];
      double lB = geom.edgeLengths[heA.edge()];
      double cosTheta = (lA * lA + lB * lB - lOpp * lOpp) / (2. * lA * lB);
      double sinTheta = 2. * geom.faceAreas[f] / (lA * lB);
      double a = weight * s_ij * cosTheta;
      double b = weight * s_ij * sinTheta;

      tripletList.emplace_back(iE_i, iE_i, weight);
      tripletList.emplace_back(iE_j, iE_j, weight);
      tripletList.emplace_back(E + iE_i, E + iE_i, weight);
      tripletList.emplace_back(E + iE_j, E + iE_j, weight);

      tripletList.emplace_back(iE_i, iE_j, a);
      tripletList.emplace_back(iE_j, iE_i, a);
      tripletList.emplace_back(E + iE_i, E + iE_j, a);
      tripletList.emplace_back(E + iE_j, E + iE_i, a);

      tripletList.emplace_back(iE_i, E + iE_j, b);
      tripletList.emplace_back(E + iE_j, iE_i, b);
      tripletList.emplace_back(iE_j, E + iE_i, -b);
      tripletList.emplace_back(E + iE_i, iE_j, -b);
    }
  }
  L.setFromTriplets(tripletList.begin(), tripletList.end());

  geom.unrequireEdgeIndices();
  geom.unrequireEdgeLengths();
  geom.unrequireHalfedgeCotanWeights();

  return L;
}

/*
 * Build the mass matrix derived in Stein et al. 2020.
 * First are all the parallel elements, then all the perpendicular ones.
 */
SparseMatrix<double> SignedHeatMethodSolver::crouzeixRaviartDoubleMassMatrix() const {

  geom.requireEdgeIndices();
  geom.requireFaceAreas();

  size_t E = mesh.nEdges();
  SparseMatrix<double> mat(2 * E, 2 * E);
  std::vector<Eigen::Triplet<double>> tripletList;
  for (Edge e : mesh.edges()) {
    size_t eIdx = geom.edgeIndices[e];
    double A = 0;
    for (Face f : e.adjacentFaces()) {
      A += geom.faceAreas[f];
    }
    double mass = A / 3.0;
    tripletList.emplace_back(eIdx, eIdx, mass);
    tripletList.emplace_back(E + eIdx, E + eIdx, mass);
  }
  mat.setFromTriplets(tripletList.begin(), tripletList.end());

  geom.unrequireEdgeIndices();
  geom.unrequireFaceAreas();

  return mat;
}