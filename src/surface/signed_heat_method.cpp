#include "geometrycentral/surface/signed_heat_method.h"

namespace geometrycentral {
namespace surface {

SignedHeatSolver::SignedHeatSolver(IntrinsicGeometryInterface& geom_, double tCoef_)
    : mesh(geom_.mesh), geom(geom_)

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
  shortTime = tCoef_ * meanNodeDistance * meanNodeDistance;

  massMat = geom.crouzeixRaviartMassMatrix;
  doubleMassMat = crouzeixRaviartDoubleMassMatrix();
  doubleConnectionLaplacian = crouzeixRaviartDoubleConnectionLaplacian();
  doubleVectorOp = doubleMassMat + shortTime * doubleConnectionLaplacian;

  geom.unrequireCrouzeixRaviartMassMatrix();
  geom.unrequireEdgeLengths();
}

VertexData<double> SignedHeatSolver::computeDistance(const std::vector<Curve>& curves,
                                                     const std::vector<SurfacePoint>& points,
                                                     const SignedHeatOptions& options) {

  Vector<std::complex<double>> Xt = integrateVectorHeatFlow(curves, points, options);
  return integrateVectorField(Xt, curves, points, options);
}

VertexData<double> SignedHeatSolver::computeDistance(const std::vector<Curve>& curves,
                                                     const SignedHeatOptions& options) {
  // call generic version
  return computeDistance(curves, std::vector<SurfacePoint>(), options);
}

VertexData<double> SignedHeatSolver::computeDistance(const std::vector<SurfacePoint>& points,
                                                     const SignedHeatOptions& options) {

  // call generic version
  return computeDistance(std::vector<Curve>(), points, options);
}

void SignedHeatSolver::setDiffusionTimeCoefficient(double tCoef_) {

  timeUpdated = true;
  shortTime = tCoef_ * meanNodeDistance * meanNodeDistance;
  doubleVectorOp = doubleMassMat + shortTime * doubleConnectionLaplacian;
}

Vector<std::complex<double>> SignedHeatSolver::integrateVectorHeatFlow(const std::vector<Curve>& curves,
                                                                       const std::vector<SurfacePoint>& points,
                                                                       const SignedHeatOptions& options) {

  ensureHaveVectorHeatSolver();

  geom.requireEdgeIndices();

  // Build source term.
  size_t E = mesh.nEdges();
  Vector<std::complex<double>> X0 = Vector<std::complex<double>>::Zero(E);
  for (const Curve& curve : curves) {
    if (curve.isSigned) {
      buildSignedCurveSource(curve, X0);
    } else {
      buildUnsignedCurveSource(curve, X0);
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
            double w = scalarCrouzeixRaviart(b, e);
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

VertexData<double> SignedHeatSolver::integrateVectorField(const Vector<std::complex<double>>& Xt,
                                                          const std::vector<Curve>& curves,
                                                          const std::vector<SurfacePoint>& points,
                                                          const SignedHeatOptions& options) {
  geom.requireHalfedgeCotanWeights();
  geom.requireVertexIndices();

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
  geom.unrequireHalfedgeCotanWeights();

  Vector<double> phi;
  switch (options.levelSetConstraint) {
  case (LevelSetConstraint::None): {
    ensureHavePoissonSolver();
    phi = poissonSolver->solve(div);
    break;
  }
  case (LevelSetConstraint::ZeroSet): {
    phi = integrateWithZeroSetConstraint(div, curves, points, options);
    break;
  }
  case (LevelSetConstraint::Multiple): {
    phi = integrateWithLevelSetConstraints(div, curves, options);
    break;
  }
  }
  phi *= -1; // since our Laplacian is positive
  if (curves.size() > 0) {
    double shift = computeAverageValueOnSource(phi, curves);
    phi -= Vector<double>::Ones(V) * shift;
  } else {
    phi -= Vector<double>::Ones(V) * phi.minCoeff();
  }

  geom.unrequireVertexIndices();

  return VertexData<double>(mesh, phi);
}

void SignedHeatSolver::buildSignedCurveSource(const Curve& curve, Vector<std::complex<double>>& X0) const {

  size_t nNodes = curve.nodes.size();
  for (size_t i = 0; i < nNodes - 1; i++) {
    const SurfacePoint& pA = curve.nodes[i];
    const SurfacePoint& pB = curve.nodes[i + 1];
    Edge commonEdge = sharedEdge(pA, pB);
    if (commonEdge != Edge()) {
      size_t eIdx = geom.edgeIndices[commonEdge];
      std::complex<double> innerProd(0, 1);
      if (pA.vertex == commonEdge.secondVertex()) innerProd.imag(-1);
      double length = lengthOfSegment(pA, pB);
      std::complex<double> contrib = length * innerProd;
      X0[eIdx] += contrib;
      continue;
    }
    Face commonFace = sharedFace(pA, pB);
    if (commonFace == Face()) {
      throw std::logic_error("Each curve segment must be contained within a single face.");
    }
    for (Edge e : commonFace.adjacentEdges()) {
      size_t eIdx = geom.edgeIndices[e];
      std::complex<double> innerProd = projectedNormal(pA, pB, e); // already incorporates length
      X0[eIdx] += innerProd;
    }
  }
}

int mod(int a, int b) { return (b + (a % b)) % b; }

void SignedHeatSolver::buildUnsignedCurveSource(const Curve& curve, Vector<std::complex<double>>& X0) {

  size_t nNodes = curve.nodes.size();

  // Add contributions from interior edges of the curve.
  for (size_t i = 0; i < nNodes - 1; i++) {
    const SurfacePoint& pA = curve.nodes[i];
    const SurfacePoint& pB = curve.nodes[i + 1];
    if (pA.type != SurfacePointType::Vertex || pB.type != SurfacePointType::Vertex) {
      throw std::logic_error("Unsigned curves not constrained to edges are not supported.");
    }
    Vertex vA = pA.vertex;
    Vertex vB = pB.vertex;
    Halfedge commonHalfedge = determineHalfedgeFromVertices(vA, vB);
    Edge commonEdge = commonHalfedge.edge();
    for (Halfedge he : commonEdge.adjacentHalfedges()) {
      // Ordinarily, "double-sided" vector information would cancel out along the edge. So in each face, "smear"
      // the info out a bit by parallel-transporting the initial vectors to other edges in adjacent faces.
      Face f = he.face();
      BarycentricVector segment = barycentricVectorInFace(he, f);
      BarycentricVector segNormal = segment.rotate90(geom);
      for (Edge e : f.adjacentEdges()) {
        if (e == commonEdge) continue;
        BarycentricVector edgeVec = barycentricVectorInFace(e.halfedge(), f);
        edgeVec /= edgeVec.norm(geom);
        double sinTheta = dot(geom, segment, edgeVec);
        double cosTheta = dot(geom, segNormal, edgeVec);
        std::complex<double> normal = std::complex<double>(cosTheta, sinTheta); // already length-weighted
        size_t eIdx = geom.edgeIndices[e];
        X0[eIdx] += normal;
      }
    }
  }

  // At curve endpoints,integrate only along the arc perpendicular to the adjacent curve segment (see Appendix.)
  Face f1_start, f2_start, f1_end, f2_end;
  Halfedge heStart_src, heEnd_src, heStart_dst, heEnd_dst;
  std::vector<std::vector<double>> angleBounds;

  // Determine starting/ending faces around curve endpoints.
  ensureHaveHalfedgeVectorsInVertex();
  geom.requireCornerScaledAngles();
  for (int i = 0; i < 2; i++) {
    Vertex vSrc, vDst;
    Vertex v;
    if (i == 0) {
      v = curve.nodes[0].vertex;
      vSrc = curve.nodes[0].vertex;
      vDst = curve.nodes[1].vertex;
    } else {
      v = curve.nodes[nNodes - 1].vertex;
      vSrc = curve.nodes[nNodes - 2].vertex;
      vDst = curve.nodes[nNodes - 1].vertex;
    }
    Halfedge segHalfedge = determineHalfedgeFromVertices(vSrc, vDst);
    Vector2 tangent;
    if (vSrc == v) {
      tangent = Vector2::fromComplex(halfedgeVectorsInVertex[segHalfedge]);
    } else {
      tangent = -Vector2::fromComplex(halfedgeVectorsInVertex[segHalfedge.twin()]);
    }
    Vector2 n1 = tangent.rotate90().normalize();
    Vector2 n2 = -n1;
    Halfedge heA = vertexTangentVectorHalfedge(v, n1);
    Halfedge heB = vertexTangentVectorHalfedge(v, n2);

    if (i == 0) {
      heStart_src = heB;
      heEnd_src = heA;

      std::complex<double> e1 = halfedgeVectorsInVertex[heStart_src];
      std::complex<double> e2 = halfedgeVectorsInVertex[heEnd_src];
      e1 /= std::abs(e1);
      e2 /= std::abs(e2);
      std::complex<double> angle1 = std::complex<double>(n1) / e1;
      std::complex<double> angle2 = std::complex<double>(n2) / e2;
      angleBounds.push_back({std::arg(angle1), std::arg(angle2)});
    } else {
      heStart_dst = heA;
      heEnd_dst = heB;

      std::complex<double> e1 = halfedgeVectorsInVertex[heStart_dst];
      std::complex<double> e2 = halfedgeVectorsInVertex[heEnd_dst];
      e1 /= std::abs(e1);
      e2 /= std::abs(e2);
      std::complex<double> angle1 = std::complex<double>(n2) / e1;
      std::complex<double> angle2 = std::complex<double>(n1) / e2;
      angleBounds.push_back({std::arg(angle1), std::arg(angle2)});
    }
  }
  geom.unrequireCornerScaledAngles();

  if (curve.nodes[0].type != SurfacePointType::Vertex || curve.nodes[nNodes - 1].type != SurfacePointType::Vertex)
    throw std::logic_error("Unsigned curves not constrained to edges are not supported.");
  std::vector<Vertex> endpoints = {curve.nodes[0].vertex, curve.nodes[nNodes - 1].vertex};
  std::vector<std::vector<Halfedge>> heBounds = {{heStart_src, heEnd_src}, {heStart_dst, heEnd_dst}};
  std::vector<double> angleSums(2);
  VertexData<Vector2> arcStarts(mesh);
  VertexData<Vector2> arcEnds(mesh);
  std::vector<double> endLengths = {
      geom.edgeLengths[determineHalfedgeFromVertices(curve.nodes[0].vertex, curve.nodes[1].vertex).edge()],

      geom.edgeLengths[determineHalfedgeFromVertices(curve.nodes[nNodes - 2].vertex, curve.nodes[nNodes - 1].vertex)
                           .edge()]};

  for (int i = 0; i < 2; i++) {

    // ========= Integrate over partial arcs

    // Determine angle sum
    double angleSum = 0.;
    Halfedge startHe = heBounds[i][0];
    Halfedge currHe = startHe.next().next().twin();
    while (currHe != heBounds[i][1]) {
      angleSum += geom.cornerAngles[currHe.corner()];
      currHe = currHe.next().next().twin();
    }
    angleSum += geom.cornerAngles[startHe.corner()] - angleBounds[i][0];
    angleSum += angleBounds[i][1];
    angleSums[i] = angleSum;

    Halfedge heA, heB, heC;
    Edge eA, eB, eOpp;
    Vertex v = heBounds[i][0].vertex();
    std::complex<double> im(0., 1.);
    std::complex<double> expA, expB, w, rot_AB, rot_AC;
    double sA, sB, sC, theta, thetaC;

    // Integrate over 2nd half of arc in starting face
    heA = heBounds[i][0];
    heB = heA.next().next();
    heC = heA.next();
    eA = heA.edge();
    eB = heB.edge();
    eOpp = heC.edge();
    theta = geom.cornerAngles[heA.corner()];
    thetaC = geom.cornerAngles[heC.corner()];
    expA = Vector2::fromAngle(angleBounds[i][0]);
    expB = Vector2::fromAngle(theta);
    w = -im * (expB - expA) / angleSum * endLengths[i];
    sA = heA.orientation() ? 1. : -1;
    sB = heB.orientation() ? 1. : -1.;
    sC = heC.orientation() ? 1. : -1.;
    rot_AB = -Vector2::fromAngle(-theta) * sB;
    rot_AC = -Vector2::fromAngle(thetaC) * sC;

    X0[geom.edgeIndices[eA]] += sA * w;
    X0[geom.edgeIndices[eB]] += rot_AB * w;
    X0[geom.edgeIndices[eOpp]] += rot_AC * w;

    // Integrate over 1st half of arc in ending face
    heA = heBounds[i][1];
    heB = heA.next().next();
    heC = heA.next();
    eA = heA.edge();
    eB = heB.edge();
    eOpp = heC.edge();
    theta = geom.cornerAngles[heA.corner()];
    thetaC = geom.cornerAngles[heC.corner()];
    expA = Vector2::fromAngle(0.);
    expB = Vector2::fromAngle(angleBounds[i][1]);
    sA = heA.orientation() ? 1. : -1;
    sB = heB.orientation() ? 1. : -1.;
    sC = heC.orientation() ? 1. : -1.;
    rot_AB = -Vector2::fromAngle(-theta) * sB;
    rot_AC = -Vector2::fromAngle(thetaC) * sC;
    w = -im * (expB - expA) / angleSum * endLengths[i];

    X0[geom.edgeIndices[eA]] += sA * w;
    X0[geom.edgeIndices[eB]] += rot_AB * w;
    X0[geom.edgeIndices[eOpp]] += rot_AC * w;

    // ========= Go over middle corners

    startHe = heBounds[i][0];
    // Go counterclockwise.
    currHe = startHe.next().next().twin();
    while (currHe != heBounds[i][1]) {
      Halfedge ij = currHe;
      Halfedge jk = ij.next();
      Halfedge ki = jk.next();
      size_t eIdx_ij = geom.edgeIndices[ij.edge()];
      size_t eIdx_jk = geom.edgeIndices[jk.edge()];
      size_t eIdx_ki = geom.edgeIndices[ki.edge()];
      Corner c = currHe.corner();
      double theta_i = geom.cornerAngles[c];
      double theta_j = geom.cornerAngles[jk.corner()];
      std::complex<double> exp_i = Vector2::fromAngle(theta_i);
      std::complex<double> im(0., 1.);
      double s_ij = ij.orientation() ? 1. : -1.;
      double s_jk = jk.orientation() ? 1. : -1.;
      double s_ki = ki.orientation() ? 1. : -1.;
      std::complex<double> r_jk_ij = -Vector2::fromAngle(theta_j) * s_jk;
      std::complex<double> r_ki_ij = -Vector2::fromAngle(-theta_i) * s_ki;
      std::complex<double> w_ij = s_ij * im * (1. - exp_i) / angleSum * endLengths[i];
      std::complex<double> w_jk = r_jk_ij * im * (1. - exp_i) / angleSum * endLengths[i];
      std::complex<double> w_ki = r_ki_ij * im * (1. - exp_i) / angleSum * endLengths[i];

      X0[eIdx_ij] += w_ij;
      X0[eIdx_jk] += w_jk;
      X0[eIdx_ki] += w_ki;

      currHe = currHe.next().next().twin();
    }
  }
}

void SignedHeatSolver::buildUnsignedPointSource(const SurfacePoint& p, Vector<std::complex<double>>& X0) const {

  switch (p.type) {
  case (SurfacePointType::Vertex): {
    Vertex v = p.vertex;
    buildUnsignedVertexSource(v, X0);
    break;
  }
  case (SurfacePointType::Edge): {
    // Ideally we'd blend post-diffusion results...
    throw std::logic_error("Point sources within edges are not supported.");
    break;
  }
  case (SurfacePointType::Face): {
    // Ideally we'd blend post-diffusion results...
    throw std::logic_error("Point sources within faces are not supported.");
    break;
  }
  default: {
    throw std::logic_error("buildUnsignedPointSource(): bad switch");
  }
  }
}

void SignedHeatSolver::buildUnsignedVertexSource(const Vertex& v, Vector<std::complex<double>>& X0,
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

Vector<double> SignedHeatSolver::integrateWithZeroSetConstraint(const Vector<double>& rhs,
                                                                const std::vector<Curve>& curves,
                                                                const std::vector<SurfacePoint>& points,
                                                                const SignedHeatOptions& options) {

  geom.requireCotanLaplacian();

  // If curve is constrained to mesh edges, solve the system by omitting DOFs.
  bool allVertices = true;
  for (const Curve& curve : curves) {
    for (const SurfacePoint& p : curve.nodes) {
      if (p.type != SurfacePointType::Vertex) {
        allVertices = false;
        break;
      }
    }
  }
  for (const SurfacePoint& p : points) {
    if (p.type != SurfacePointType::Vertex) {
      allVertices = false;
      break;
    }
  }

  size_t V = mesh.nVertices();
  Vector<double> phi;
  if (allVertices && options.softLevelSetWeight < 0.) {
    Vector<bool> setAMembership = Vector<bool>::Ones(V); // true if interior
    for (const Curve& curve : curves) {
      for (const SurfacePoint& p : curve.nodes) {
        Vertex v = p.vertex;
        setAMembership[geom.vertexIndices[v]] = false;
      }
    }
    for (const SurfacePoint& p : points) {
      Vertex v = p.vertex;
      setAMembership[geom.vertexIndices[v]] = false;
    }
    int nB = V - setAMembership.cast<int>().sum(); // # of boundary vertices
    Vector<double> bcVals = Vector<double>::Zero(nB);
    BlockDecompositionResult<double> decomp = blockDecomposeSquare(geom.cotanLaplacian, setAMembership, true);
    Vector<double> rhsValsA, rhsValsB;
    decomposeVector(decomp, rhs, rhsValsA, rhsValsB);
    Vector<double> combinedRHS = rhsValsA;
    shiftDiagonal(decomp.AA, 1e-16);
    PositiveDefiniteSolver<double> constrainedSolver(decomp.AA);
    Vector<double> Aresult = constrainedSolver.solve(combinedRHS);
    phi = reassembleVector(decomp, Aresult, bcVals);
  } else {
    Curve allCurves;
    for (const Curve& curve : curves) {
      allCurves.nodes.insert(allCurves.nodes.end(), curve.nodes.begin(), curve.nodes.end());
    }
    for (const SurfacePoint& p : points) allCurves.nodes.push_back(p);
    phi = integrateWithLevelSetConstraints(rhs, {allCurves}, options);
  }

  geom.unrequireCotanLaplacian();

  return phi;
}

Vector<double> SignedHeatSolver::integrateWithLevelSetConstraints(const Vector<double>& rhs,
                                                                  const std::vector<Curve>& curves,
                                                                  const SignedHeatOptions& options) {

  geom.requireCotanLaplacian();

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
    } else if (p0.type == SurfacePointType::Face) {
      int idx = 0;
      for (Vertex v : p0.face.adjacentVertices()) {
        size_t vIdx = geom.vertexIndices[v];
        p0_data.emplace_back(m, vIdx, p0.faceCoords[idx]);
        idx++;
      }
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
      } else if (p.type == SurfacePointType::Face) {
        int idx = 0;
        for (Vertex v : p.face.adjacentVertices()) {
          size_t vIdx = geom.vertexIndices[v];
          triplets.emplace_back(m, vIdx, -p.faceCoords[idx]);
          idx++;
        }
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

  return phi;
}

double SignedHeatSolver::lengthOfSegment(const SurfacePoint& pA, const SurfacePoint& pB) const {

  BarycentricVector segment(pA, pB);
  return segment.norm(geom);
}

SurfacePoint SignedHeatSolver::midSegmentSurfacePoint(const SurfacePoint& pA, const SurfacePoint& pB) const {

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
std::complex<double> SignedHeatSolver::projectedNormal(const SurfacePoint& pA, const SurfacePoint& pB,
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

BarycentricVector SignedHeatSolver::barycentricVectorInFace(const Halfedge& he_, const Face& f) const {

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
  faceCoords[(eIdx + 1) % 3] = 1;
  faceCoords[eIdx] = -1;
  return BarycentricVector(f, sign * faceCoords);
}

FaceData<BarycentricVector> SignedHeatSolver::sampleAtFaceBarycenters(const Vector<std::complex<double>>& Xt) {

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
    double magn = vec.norm(geom);
    X[f] = (magn == 0) ? BarycentricVector(f, {0, 0, 0}) : vec / magn; // normalize
  }

  geom.unrequireEdgeIndices();

  return X;
}

double SignedHeatSolver::computeAverageValueOnSource(const Vector<double>& phi,
                                                     const std::vector<Curve>& curves) const {

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

  return shift;
}

void SignedHeatSolver::ensureHaveVectorHeatSolver() {

  if (vectorHeatSolver != nullptr && !timeUpdated) return;
  timeUpdated = false;

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
}

void SignedHeatSolver::ensureHavePoissonSolver() {

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
SparseMatrix<double> SignedHeatSolver::crouzeixRaviartDoubleConnectionLaplacian() const {

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
SparseMatrix<double> SignedHeatSolver::crouzeixRaviartDoubleMassMatrix() const {

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

Halfedge SignedHeatSolver::determineHalfedgeFromVertices(const Vertex& vA, const Vertex& vB) const {

  for (Halfedge he : vA.outgoingHalfedges()) {
    if (he.tipVertex() == vB) return he;
  }
  return Halfedge();
}

/*
 * Populate halfedgeVectorsInVertex() for a general SurfaceMesh. May not be valid at non-manifold vertices, or on
 * non-orientable surfaces.
 */
void SignedHeatSolver::ensureHaveHalfedgeVectorsInVertex() {

  if (halfedgeVectorsInVertex.size() > 0) return;

  geom.requireEdgeLengths();
  geom.requireCornerScaledAngles();

  halfedgeVectorsInVertex = HalfedgeData<std::complex<double>>(mesh);

  for (Vertex v : mesh.vertices()) {
    double coordSum = 0.0;

    // orbit CCW
    Halfedge firstHe = v.halfedge();
    Halfedge currHe = firstHe;
    do {
      halfedgeVectorsInVertex[currHe] =
          std::complex<double>(Vector2::fromAngle(coordSum)) * geom.edgeLengths[currHe.edge()];
      coordSum += geom.cornerScaledAngles[currHe.corner()];
      if (!currHe.isInterior()) break;
      currHe = currHe.next().next().twin();
    } while (currHe != firstHe);
  }

  geom.unrequireEdgeLengths();
  geom.unrequireCornerScaledAngles();
}

Halfedge SignedHeatSolver::vertexTangentVectorHalfedge(const Vertex& v, const Vector2& vec) const {

  double theta = vec.arg() + M_PI; // in [-pi, pi]
  double angle = 0.;
  Halfedge startHe = v.halfedge();
  Halfedge currHe = startHe;
  do {
    Corner c = currHe.corner();
    angle += geom.cornerScaledAngles[c];
    if (angle >= theta) return currHe;
    currHe = currHe.next().next().twin();
  } while (currHe != startHe);
  throw std::logic_error("vertexTangentVectorHalfedge(): something went wrong");
}

double SignedHeatSolver::scalarCrouzeixRaviart(const SurfacePoint& p, const Edge& e) const {
  // Here I used casework. But we could also find the barycentric coordinates of `p` in the face shared by `p` and
  // `e`, and compute 1 - 2b_i, where b_i is the barycentric coordinate of `p` assoc. with the vertex across from `e`.
  Face sharedFace = Face();
  switch (p.type) {
  case (SurfacePointType::Vertex): {
    for (Face f : e.adjacentFaces()) {
      for (Vertex v : f.adjacentVertices()) {
        if (v == p.vertex) {
          sharedFace = f;
          break;
        }
      }
    }
    if (sharedFace == Face()) return 0.;
    if (p.vertex == e.firstVertex() || p.vertex == e.secondVertex()) return 1.;
    return -1;
    break;
  }
  case (SurfacePointType::Edge): {
    if (p.edge == e) return 1.;
    for (Face f : e.adjacentFaces()) {
      for (Edge ep : f.adjacentEdges()) {
        if (ep == p.edge) {
          sharedFace = f;
          break;
        }
      }
    }
    if (sharedFace == Face()) return 0.;
    SurfacePoint pInFace = p.inFace(sharedFace);
    int idx = 0;
    for (Halfedge he : sharedFace.adjacentHalfedges()) {
      if (he.next().edge() == e) return 1. - 2. * pInFace.faceCoords[idx];
      idx++;
    }
    break;
  }
  case (SurfacePointType::Face): {
    for (Face f : e.adjacentFaces()) {
      if (f == p.face) {
        sharedFace = f;
        break;
      }
    }
    if (sharedFace == Face()) return 0.;
    size_t idx = 0;
    double val = 0.;
    for (Halfedge he : sharedFace.adjacentHalfedges()) {
      double vertexValue = ((he.vertex() == e.firstVertex()) || (he.vertex() == e.secondVertex())) ? 1 : -1;
      double vertexWeight = p.faceCoords[idx];
      val += vertexValue * vertexWeight;
      idx++;
    }
    return val;
    break;
  }
  }
  throw std::logic_error("scalarCrouzeixRaviart(): bad switch");
  return 0.;
}

} // namespace surface
} // namespace geometrycentral