#include "geometrycentral/surface/embedded_geometry_interface.h"

#include <limits>

using std::cout;
using std::endl;

namespace geometrycentral {
namespace surface {

// clang-format off
EmbeddedGeometryInterface::EmbeddedGeometryInterface(SurfaceMesh& mesh_) : 
  ExtrinsicGeometryInterface(mesh_),

  vertexPositionsQ                (&vertexPositions,                std::bind(&EmbeddedGeometryInterface::computeVertexPositions, this),                quantities),
  faceNormalsQ                    (&faceNormals,                    std::bind(&EmbeddedGeometryInterface::computeFaceNormals, this),                    quantities),
  vertexNormalsQ                  (&vertexNormals,                  std::bind(&EmbeddedGeometryInterface::computeVertexNormals, this),                  quantities),
  faceTangentBasisQ               (&faceTangentBasis,               std::bind(&EmbeddedGeometryInterface::computeFaceTangentBasis, this),               quantities),
  vertexTangentBasisQ             (&vertexTangentBasis,             std::bind(&EmbeddedGeometryInterface::computeVertexTangentBasis, this),             quantities),
  vertexDualMeanCurvatureNormalsQ (&vertexDualMeanCurvatureNormals, std::bind(&EmbeddedGeometryInterface::computeVertexDualMeanCurvatureNormals, this), quantities),

  virtualRefinementAreaWeightsQ              (&virtualRefinementAreaWeights,              std::bind(&EmbeddedGeometryInterface::computeVirtualRefinementAreaWeights, this),              quantities),
  virtualRefinementLaplacianQ                (&virtualRefinementLaplacian,                std::bind(&EmbeddedGeometryInterface::computeVirtualRefinementLaplacian, this),                quantities),
  virtualRefinementVertexGalerkinMassMatrixQ (&virtualRefinementVertexGalerkinMassMatrix, std::bind(&EmbeddedGeometryInterface::computeVirtualRefinementVertexGalerkinMassMatrix, this), quantities),
  virtualRefinementVertexLumpedMassMatrixQ   (&virtualRefinementVertexLumpedMassMatrix,   std::bind(&EmbeddedGeometryInterface::computeVirtualRefinementVertexLumpedMassMatrix, this),   quantities),
 
  virtualElementLaplacianQ                 (&virtualElementLaplacian,                 std::bind(&EmbeddedGeometryInterface::computeVirtualElementLaplacian, this),                 quantities),
  virtualElementVertexGalerkinMassMatrixQ  (&virtualElementVertexGalerkinMassMatrix,  std::bind(&EmbeddedGeometryInterface::computeVirtualElementVertexGalerkinMassMatrix, this),  quantities),
  virtualElementVertexLumpedMassMatrixQ    (&virtualElementVertexLumpedMassMatrix,    std::bind(&EmbeddedGeometryInterface::computeVirtualElementVertexLumpedMassMatrix, this),    quantities),
  virtualElementVertexConnectionLaplacianQ (&virtualElementVertexConnectionLaplacian, std::bind(&EmbeddedGeometryInterface::computeVirtualElementVertexConnectionLaplacian, this), quantities),

  virtualRefinementDECOperatorArray{&virtualRefinementHodge0, &virtualRefinementHodge0Inverse, &virtualRefinementHodge1, 
                                    &virtualRefinementHodge1Inverse, &virtualRefinementHodge2, 
                                    &virtualRefinementHodge2Inverse, &virtualRefinementD0, &virtualRefinementD1},
  virtualRefinementDECOperatorsQ(&virtualRefinementDECOperatorArray, std::bind(&EmbeddedGeometryInterface::computeVirtualRefinementDECOperators, this), quantities),
  
  virtualElementDECOperatorArray{&virtualElementHodge0, &virtualElementHodge0Inverse, &virtualElementHodge1, 
                                 &virtualElementHodge1Inverse, &virtualElementHodge2, &virtualElementHodge2Inverse, 
                                 &virtualElementD0, &virtualElementD1},
  virtualElementDECOperatorsQ(&virtualElementDECOperatorArray, std::bind(&EmbeddedGeometryInterface::computeVirtualElementDECOperators, this), quantities)
  
  {}
// clang-format on

// === Overrides

// Edge lengths
void EmbeddedGeometryInterface::computeEdgeLengths() {
  vertexPositionsQ.ensureHave();

  edgeLengths = EdgeData<double>(mesh);
  for (Edge e : mesh.edges()) {
    edgeLengths[e] = norm(vertexPositions[e.halfedge().vertex()] - vertexPositions[e.halfedge().next().vertex()]);
  }
}

// Edge dihedral angles
void EmbeddedGeometryInterface::computeEdgeDihedralAngles() {
  vertexPositionsQ.ensureHave();
  faceNormalsQ.ensureHave();

  edgeDihedralAngles = EdgeData<double>(mesh, 0.);
  for (Edge e : mesh.edges()) {
    if (e.isBoundary()) continue;

    if (!e.isManifold()) {
      continue;
    }

    Vector3 N1 = faceNormals[e.halfedge().face()];
    Vector3 N2 = faceNormals[e.halfedge().sibling().face()];
    Vector3 pTail = vertexPositions[e.halfedge().vertex()];
    Vector3 pTip = vertexPositions[e.halfedge().next().vertex()];
    Vector3 edgeDir = unit(pTip - pTail);

    edgeDihedralAngles[e] = atan2(dot(edgeDir, cross(N1, N2)), dot(N1, N2));
  }
}

// === Quantities

void EmbeddedGeometryInterface::requireVertexPositions() { vertexPositionsQ.require(); }
void EmbeddedGeometryInterface::unrequireVertexPositions() { vertexPositionsQ.unrequire(); }


void EmbeddedGeometryInterface::computeFaceNormals() {
  vertexPositionsQ.ensureHave();

  faceNormals = FaceData<Vector3>(mesh);

  for (Face f : mesh.faces()) {

    // For general polygons, take the sum of the cross products at each corner
    Vector3 normalSum = Vector3::zero();
    for (Halfedge heF : f.adjacentHalfedges()) {

      // Gather vertex positions for next three vertices
      Halfedge he = heF;
      Vector3 pA = vertexPositions[he.vertex()];
      he = he.next();
      Vector3 pB = vertexPositions[he.vertex()];
      he = he.next();
      Vector3 pC = vertexPositions[he.vertex()];

      normalSum += cross(pB - pA, pC - pA);

      // In the special case of a triangle, there is no need to to repeat at all three corners; the result will be the
      // same
      if (he.next() == heF) break;
    }

    Vector3 normal = unit(normalSum);
    faceNormals[f] = normal;
  }
}
void EmbeddedGeometryInterface::requireFaceNormals() { faceNormalsQ.require(); }
void EmbeddedGeometryInterface::unrequireFaceNormals() { faceNormalsQ.unrequire(); }

// Vertex normal
void EmbeddedGeometryInterface::computeVertexNormals() {
  faceNormalsQ.ensureHave();
  cornerAnglesQ.ensureHave();

  vertexNormals = VertexData<Vector3>(mesh);

  for (Vertex v : mesh.vertices()) {
    Vector3 normalSum = Vector3::zero();

    for (Corner c : v.adjacentCorners()) {
      Vector3 normal = faceNormals[c.face()];
      double weight = cornerAngles[c];

      normalSum += weight * normal;
    }

    vertexNormals[v] = unit(normalSum);
  }
}
void EmbeddedGeometryInterface::requireVertexNormals() { vertexNormalsQ.require(); }
void EmbeddedGeometryInterface::unrequireVertexNormals() { vertexNormalsQ.unrequire(); }

// Face tangent basis
void EmbeddedGeometryInterface::computeFaceTangentBasis() {
  vertexPositionsQ.ensureHave();
  faceNormalsQ.ensureHave();

  faceTangentBasis = FaceData<std::array<Vector3, 2>>(mesh);

  if (!mesh.usesImplicitTwin()) {
    // For a nonmanifold mesh, just compute any extrinsic basis
    for (Face f : mesh.faces()) {
      Vector3 normal = faceNormals[f];
      faceTangentBasis[f] = normal.buildTangentBasis();
    }
    return;
  }

  halfedgeVectorsInFaceQ.ensureHave();

  for (Face f : mesh.faces()) {
    // TODO this implementation seems a bit silly...

    // For general polygons, take the average of each edge vector projected to tangent plane
    Vector3 basisXSum = Vector3::zero();
    Vector3 N = faceNormals[f];
    bool isTriangular = f.isTriangle();
    for (Halfedge heF : f.adjacentHalfedges()) {

      Vector3 eVec = vertexPositions[heF.next().vertex()] - vertexPositions[heF.vertex()];
      eVec = eVec.removeComponent(N);

      double angle = halfedgeVectorsInFace[heF].arg();
      Vector3 eVecX = eVec.rotateAround(N, -angle);

      basisXSum += eVecX;

      // In the special case of a triangle, there is no need to to repeat at all three corners; the result will be the
      // same
      if (isTriangular) break;
    }

    Vector3 basisX = unit(basisXSum);
    Vector3 basisY = cross(N, basisX);
    faceTangentBasis[f] = {{basisX, basisY}};
  }
}
void EmbeddedGeometryInterface::requireFaceTangentBasis() { faceTangentBasisQ.require(); }
void EmbeddedGeometryInterface::unrequireFaceTangentBasis() { faceTangentBasisQ.unrequire(); }

// Vertex tangent basis
void EmbeddedGeometryInterface::computeVertexTangentBasis() {
  vertexPositionsQ.ensureHave();
  vertexNormalsQ.ensureHave();

  vertexTangentBasis = VertexData<std::array<Vector3, 2>>(mesh);

  if (!mesh.usesImplicitTwin()) {
    // For a nonmanifold mesh, just compute any extrinsic basis
    for (Vertex v : mesh.vertices()) {
      Vector3 normal = vertexNormals[v];
      vertexTangentBasis[v] = normal.buildTangentBasis();
    }
    return;
  }

  halfedgeVectorsInVertexQ.ensureHave();

  for (Vertex v : mesh.vertices()) {

    // For general polygons, take the average of each edge vector projected to tangent plane
    Vector3 basisXSum = Vector3::zero();
    Vector3 N = vertexNormals[v];
    for (Halfedge he : v.outgoingHalfedges()) {

      Vector3 eVec = vertexPositions[he.next().vertex()] - vertexPositions[he.vertex()];
      eVec = eVec.removeComponent(N);

      // TODO can surely do this with less trig
      double angle = halfedgeVectorsInVertex[he].arg();
      Vector3 eVecX = eVec.rotateAround(N, -angle);

      basisXSum += eVecX;
    }

    Vector3 basisX = unit(basisXSum);
    Vector3 basisY = cross(N, basisX);
    vertexTangentBasis[v] = {{basisX, basisY}};
  }
}
void EmbeddedGeometryInterface::requireVertexTangentBasis() { vertexTangentBasisQ.require(); }
void EmbeddedGeometryInterface::unrequireVertexTangentBasis() { vertexTangentBasisQ.unrequire(); }

void EmbeddedGeometryInterface::computeVertexDualMeanCurvatureNormals() {
  edgeCotanWeightsQ.ensureHave();
  vertexPositionsQ.ensureHave();

  vertexDualMeanCurvatureNormals = VertexData<Vector3>(mesh, Vector3::zero());

  // These are defined by the property that the mean curvature normals are the (half) the laplacian of the vertex
  // positions
  // WARNING: this means that vertexMeanCurvatures != vertexMeanCurvatureNormals.norm()

  // Rather than building the whole cotan laplacian, we evaluate 0.5 * L * positions using our edge cotan weights
  for (Edge e : mesh.edges()) {
    double w = edgeCotanWeights[e];

    Vertex vTail = e.halfedge().tailVertex();
    Vertex vTip = e.halfedge().tipVertex();

    Vector3 pTail = vertexPositions[vTail];
    Vector3 pTip = vertexPositions[vTip];

    vertexDualMeanCurvatureNormals[vTail] += w * (pTail - pTip) / 2;
    vertexDualMeanCurvatureNormals[vTip] += w * (pTip - pTail) / 2;
  }
}
void EmbeddedGeometryInterface::requireVertexDualMeanCurvatureNormals() { vertexDualMeanCurvatureNormalsQ.require(); }
void EmbeddedGeometryInterface::unrequireVertexDualMeanCurvatureNormals() {
  vertexDualMeanCurvatureNormalsQ.unrequire();
}


// == Overrides to compute things better using vertex positions

// Override to compute directly from vertex positions
void EmbeddedGeometryInterface::computeFaceAreas() {
  vertexPositionsQ.ensureHave();

  faceAreas = FaceData<double>(mesh);

  for (Face f : mesh.faces()) {

    // WARNING: Logic duplicated between cached and immediate version
    Halfedge he = f.halfedge();
    Vector3 pA = vertexPositions[he.vertex()];
    he = he.next();
    Vector3 pB = vertexPositions[he.vertex()];
    he = he.next();
    Vector3 pC = vertexPositions[he.vertex()];

    GC_SAFETY_ASSERT(he.next() == f.halfedge(), "faces must be triangular");

    double area = 0.5 * norm(cross(pB - pA, pC - pA));
    faceAreas[f] = area;
  }
}

// Override to compute directly from vertex positions
void EmbeddedGeometryInterface::computeCornerAngles() {
  vertexPositionsQ.ensureHave();

  cornerAngles = CornerData<double>(mesh);

  for (Corner c : mesh.corners()) {

    // WARNING: Logic duplicated between cached and immediate version
    Halfedge he = c.halfedge();
    Vector3 pA = vertexPositions[he.vertex()];
    he = he.next();
    Vector3 pB = vertexPositions[he.vertex()];
    he = he.next();
    Vector3 pC = vertexPositions[he.vertex()];

    GC_SAFETY_ASSERT(he.next() == c.halfedge(), "faces must be triangular");

    double q = dot(unit(pB - pA), unit(pC - pA));
    q = clamp(q, -1.0, 1.0);
    double angle = std::acos(q);

    cornerAngles[c] = angle;
  }
}


// Override to compute directly from vertex positions
void EmbeddedGeometryInterface::computeHalfedgeCotanWeights() {
  vertexPositionsQ.ensureHave();

  halfedgeCotanWeights = HalfedgeData<double>(mesh);

  for (Halfedge heI : mesh.interiorHalfedges()) {
    // WARNING: Logic duplicated between cached and immediate version

    Halfedge he = heI;
    Vector3 pB = vertexPositions[he.vertex()];
    he = he.next();
    Vector3 pC = vertexPositions[he.vertex()];
    he = he.next();
    Vector3 pA = vertexPositions[he.vertex()];
    GC_SAFETY_ASSERT(he.next() == heI, "faces must be triangular");

    Vector3 vecR = pB - pA;
    Vector3 vecL = pC - pA;

    double cotValue = dot(vecR, vecL) / norm(cross(vecR, vecL));

    halfedgeCotanWeights[heI] = cotValue / 2;
  }
}


// Override to compute directly from vertex positions
void EmbeddedGeometryInterface::computeEdgeCotanWeights() {
  vertexPositionsQ.ensureHave();

  edgeCotanWeights = EdgeData<double>(mesh);

  for (Edge e : mesh.edges()) {
    double cotSum = 0.;

    for (Halfedge he : e.adjacentInteriorHalfedges()) {
      // WARNING: Logic duplicated between cached and immediate version
      Halfedge heFirst = he;
      Vector3 pB = vertexPositions[he.vertex()];
      he = he.next();
      Vector3 pC = vertexPositions[he.vertex()];
      he = he.next();
      Vector3 pA = vertexPositions[he.vertex()];
      GC_SAFETY_ASSERT(he.next() == heFirst, "faces must be triangular");

      Vector3 vecR = pB - pA;
      Vector3 vecL = pC - pA;

      double cotValue = dot(vecR, vecL) / norm(cross(vecR, vecL));
      cotSum += cotValue / 2;
    }

    edgeCotanWeights[e] = cotSum;
  }
}

// === Polygon Operators

// = Bunge et al. "Polygon Laplacian Made Simple" (2020), based on virtual refinement (virtual node method).

// Laplacian
void EmbeddedGeometryInterface::computeVirtualRefinementLaplacian() {
  vertexIndicesQ.ensureHave();

  size_t V = mesh.nVertices();
  virtualRefinementLaplacian = Eigen::SparseMatrix<double>(V, V);
  std::vector<Eigen::Triplet<double>> triplets;
  std::vector<size_t> vIndices; // indices of vertices of polygon face
  Eigen::MatrixXd Si;           // local per-polygon matrix
  for (Face f : mesh.faces()) {
    vIndices.clear();
    for (Vertex v : f.adjacentVertices()) vIndices.push_back(vertexIndices[v]);
    int n = f.degree();
    // Get local stiffness matrix.
    Si = buildPolygonStiffnessMatrix(f);
    // Add contribution to global mass matrix.
    for (int j = 0; j < n; j++) {
      for (int i = 0; i < n; i++) {
        triplets.emplace_back(vIndices[i], vIndices[j], Si(i, j));
      }
    }
  }
  virtualRefinementLaplacian.setFromTriplets(triplets.begin(), triplets.end());
}
void EmbeddedGeometryInterface::requireVirtualRefinementLaplacian() { virtualRefinementLaplacianQ.require(); }
void EmbeddedGeometryInterface::unrequireVirtualRefinementLaplacian() { virtualRefinementLaplacianQ.unrequire(); }


// Vertex Galerkin mass matrix (unlumped)
void EmbeddedGeometryInterface::computeVirtualRefinementVertexGalerkinMassMatrix() {
  vertexIndicesQ.ensureHave();

  size_t V = mesh.nVertices();
  virtualRefinementVertexGalerkinMassMatrix = Eigen::SparseMatrix<double>(V, V);
  std::vector<Eigen::Triplet<double>> triplets;
  std::vector<size_t> vIndices; // indices of vertices of polygon face
  Eigen::MatrixXd Mi;           // local per-polygon matrix
  for (Face f : mesh.faces()) {
    vIndices.clear();
    for (Vertex v : f.adjacentVertices()) vIndices.push_back(vertexIndices[v]);
    int n = f.degree();
    // Get local mass matrix.
    Mi = buildPolygonMassMatrix(f);
    // Add contribution to global mass matrix.
    for (int j = 0; j < n; j++) {
      for (int i = 0; i < n; i++) {
        triplets.emplace_back(vIndices[i], vIndices[j], Si(i, j));
      }
    }
  }
  virtualRefinementVertexGalerkinMassMatrix.setFromTriplets(triplets.begin(), triplets.end());
}
void EmbeddedGeometryInterface::requireVirtualRefinementVertexGalerkinMassMatrix() {
  virtualRefinementVertexGalerkinMassMatrixQ.require();
}
void EmbeddedGeometryInterface::unrequireVirtualRefinementVertexGalerkinMassMatrix() {
  virtualRefinementVertexGalerkinMassMatrixQ.unrequire();
}


// Vertex mass matrix (lumped)
void EmbeddedGeometryInterface::computeVirtualRefinementVertexLumpedMassMatrix() {
  vertexIndicesQ.ensureHave();

  size_t V = mesh.nVertices();
  virtualRefinementVertexGalerkinMassMatrix = Eigen::SparseMatrix<double>(V, V);
  std::vector<Eigen::Triplet<double>> triplets;
  std::vector<size_t> vIndices; // indices of vertices of polygon face
  Eigen::MatrixXd Mi;           // local per-polygon matrix
  for (Face f : mesh.faces()) {
    vIndices.clear();
    for (Vertex v : f.adjacentVertices()) vIndices.push_back(vertexIndices[v]);
    int n = f.degree();
    // Get local mass matrix.
    Mi = buildPolygonMassMatrix(f);
    // Add contribution to global mass matrix.
    for (int j = 0; j < n; j++) {
      for (int i = 0; i < n; i++) {
        triplets.emplace_back(vIndices[i], vIndices[i], Si(i, j));
      }
    }
  }
  virtualRefinementVertexGalerkinMassMatrix.setFromTriplets(triplets.begin(), triplets.end());
}
void EmbeddedGeometryInterface::requireVirtualRefinementVertexLumpedMassMatrix() {
  virtualRefinementVertexLumpedMassMatrixQ.require();
}
void EmbeddedGeometryInterface::unrequireVirtualRefinementVertexLumpedMassMatrix() {
  virtualRefinementVertexLumpedMassMatrixQ.unrequire();
}


// DEC Operators
void EmbeddedGeometryInterface::computeVirtualRefinementDECOperators() {
  // TODO
}
void EmbeddedGeometryInterface::requireVirtualRefinementDECOperators() { virtualRefinementDECOperatorsQ.require(); }
void EmbeddedGeometryInterface::unrequireVirtualRefinementDECOperators() { virtualRefinementDECOperatorsQ.unrequire(); }


// Helper functions

Eigen::MatrixXd EmbeddedGeometryInterface::buildPolygonMassMatrix(const Face& f) const {
  virtualRefinementAreaWeightsQ.ensureHave();

  int n = f.degree();
  Eigen::MatrixXd poly = getPolygonPositionMatrix(f);
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(n, n);
  const Eigen::VectorXd& weights = virtualRefinementAreaWeights[f];
  Eigen::Vector3d virtualVertex = poly.transpose() * weights;
  Eigen::VectorXd ln = Eigen::VectorXd::Zero(n + 1);
  double l[3], l2[3]; // lengths, lengths squared
  // Build triangle fan mass and cotan matrices
  for (int i = 0; i < n; i++) {
    const int i1 = (i + 1) % n;
    l2[2] = (poly.row(i) - poly.row(i1)).squaredNorm();
    l2[0] = (poly.row(i1) - virtualVertex.transpose()).squaredNorm();
    l2[1] = (poly.row(i) - virtualVertex.transpose()).squaredNorm();
    l[0] = std::sqrt(l2[0]);
    l[1] = std::sqrt(l2[1]);
    l[2] = std::sqrt(l2[2]);
    const double arg =
        (l[0] + (l[1] + l[2])) * (l[2] - (l[0] - l[1])) * (l[2] + (l[0] - l[1])) * (l[0] + (l[1] - l[2]));
    const double area = 0.25 * std::sqrt(arg);
    l[0] = 1.0 / 6.0 * area;
    l[1] = 1.0 / 12.0 * area;
    M(i1, i1) += 1.0 / 6.0 * area;
    M(i, i) += 1.0 / 6.0 * area;
    M(i1, i) += 1.0 / 12.0 * area;
    M(i, i1) += 1.0 / 12.0 * area;
    ln(i1) += l[1];
    ln(i) += l[1];
    ln(n) += l[0];
  }
  // Apply prolongation
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i) M(i, j) += weights(i) * ln(j) + weights(j) * ln(i) + weights(i) * weights(j) * ln(n);

  return M;
}

Eigen::MatrixXd EmbeddedGeometryInterface::buildPolygonStiffnessMatrix(const Face& f) const {
  virtualRefinementAreaWeightsQ.ensureHave();

  int n = f.degree();
  Eigen::MatrixXd poly = getPolygonPositionMatrix(f);
  Eigen::MatrixXd S = Eigen::MatrixXd::Zero(n, n);
  const Eigen::VectorXd& weights = virtualRefinementAreaWeights[f];
  Eigen::Vector3d virtualVertex = poly.transpose() * weights;
  Eigen::VectorXd ln = Eigen::VectorXd::Zero(n + 1);
  double l[3], l2[3]; // lengths, lengths squared
  // Build triangle fan mass and cotan matrices
  for (int i = 0; i < n; i++) {
    const int i1 = (i + 1) % n;
    l2[2] = (poly.row(i) - poly.row(i1)).squaredNorm();
    l2[0] = (poly.row(i1) - virtualVertex.transpose()).squaredNorm();
    l2[1] = (poly.row(i) - virtualVertex.transpose()).squaredNorm();
    l[0] = std::sqrt(l2[0]);
    l[1] = std::sqrt(l2[1]);
    l[2] = std::sqrt(l2[2]);
    const double arg =
        (l[0] + (l[1] + l[2])) * (l[2] - (l[0] - l[1])) * (l[2] + (l[0] - l[1])) * (l[0] + (l[1] - l[2]));
    const double area = 0.5 * std::sqrt(arg);
    if (area > 1e-7) {
      l[0] = 0.25 * (l2[1] + l2[2] - l2[0]) / area;
      l[1] = 0.25 * (l2[2] + l2[0] - l2[1]) / area;
      l[2] = 0.25 * (l2[0] + l2[1] - l2[2]) / area;

      S(i1, i1) += l[0];
      S(i, i) += l[1];
      S(i1, i) -= l[2];
      S(i, i1) -= l[2];
      S(i, i) += l[2];
      S(i1, i1) += l[2];

      ln(i1) -= l[0];
      ln(i) -= l[1];
      ln(n) += l[0] + l[1];
    }
  }
  // Apply prolongation
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i) S(i, j) += weights(i) * ln(j) + weights(j) * ln(i) + weights(i) * weights(j) * ln(n);

  return S;
}

SparseMatrix<double> EmbeddedGeometryInterface::buildDivergenceMatrix() const {
  SparseMatrix<double> G = buildGradientMatrix();
  SparseMatrix<double> M = buildGradientMassMatrix();
  return -G.transpose() * M; // take convention that outflow is positive
}

SparseMatrix<double> EmbeddedGeometryInterface::buildGradientMatrix() const {
  virtualRefinementAreaWeightsQ.ensureHave();
  vertexIndicesQ.ensureHave();
  vertexPositionsQ.ensureHave();
  faceIndicesQ.ensureHave();

  // Builds a matrix G ∈ (R^3)^{|T^f| x |V|}, which gets built as a 3|T^f| x |V| matrix.
  size_t V = mesh.nVertices();
  size_t F = mesh.nFaces();
  SparseMatrix<double> G;
  std::vector<Eigen::Triplet<double>> triplets;
  size_t nTriangles = 0;
  size_t k = 0;
  for (Face f : mesh.faces()) {
    size_t fIdx = faceIndices[f];
    nTriangles += f.degree();
    Eigen::MatrixXd poly = getPolygonPositionMatrix(f);
    const Eigen::VectorXd& weights = virtualRefinementAreaWeights[f];
    Eigen::Vector3d p = poly.transpose() * weights;
    for (Halfedge he : f.adjacentHalfedges()) {
      size_t v0 = vertexIndices[he.tailVertex()];
      size_t v1 = vertexIndices[he.tipVertex()];
      Vector3 p0 = vertexPositions[v0];
      Vector3 p1 = vertexPositions[v1];
      Vector3 grad_p = gradientHatFunction(p, p0, p1);
      Vector3 grad_p0 = gradientHatFunction(p0, p1, p);
      Vector3 grad_p1 = gradientHatFunction(p1, p, p0);
      for (size_t j = 0; j < 3; j++) {
        triplets.emplace_back(3 * k + j, V + fIdx, grad_p[j]);
        triplets.emplace_back(3 * k + j, v0, grad_p0[j]);
        triplets.emplace_back(3 * k + j, v1, grad_p1[j]);
      }
      k++;
    }
  }
  G.resize(3 * nTriangles, V + F);
  G.setFromTriplets(triplets.begin(), triplets.end());
  SparseMatrix<double> P = buildProlongationMatrix();
  G = G * P;
  return G;
}

SparseMatrix<double> EmbeddedGeometryInterface::buildGradientMassMatrix() const {
  virtualRefinementAreaWeightsQ.ensureHave();
  vertexPositionsQ.ensureHave();

  // Block diagonal matrix whose i-th block consists of the 3×3 identity matrix multiplied by the area of the i-th
  // triangle.
  SparseMatrix<double> M;
  std::vector<Eigen::Triplet<double>> triplets;
  size_t c = 0;
  for (Face f : mesh.faces()) {
    Eigen::MatrixXd poly = getPolygonPositionMatrix(f);
    const Eigen::VectorXd& weights = virtualRefinementAreaWeights[f];
    Eigen::Vector3d areaPoint = poly.transpose() * weights;
    Vector3 ap = {areaPoint[0], areaPoint[1], areaPoint[2]};
    size_t i = 0;
    for (Halfedge he : f.adjacentHalfedges()) {
      Vector3 p0 = vertexPositions[he.tailVertex()];
      Vector3 p1 = vertexPositions[he.tipVertex()];
      double area = 0.5 * cross(p0 - ap, p1 - ap).norm();
      for (size_t j = 0; j < 3; j++) {
        size_t idx = c + 3 * i + j;
        triplets.emplace_back(idx, idx, area);
      }
      i++;
    }
    c += f.degree() * 3;
  }
  M.resize(c, c);
  M.setFromTriplets(triplets.begin(), triplets.end());
  return M;
}

SparseMatrix<double> EmbeddedGeometryInterface::buildProlongationMatrix() const {
  virtualRefinementAreaWeightsQ.ensureHave();
  vertexIndicesQ.ensureHave();

  size_t V = mesh.nVertices();
  size_t F = mesh.nFaces();
  std::vector<Eigen::Triplet<double>> triplets;
  SparseMatrix<double> P(V + F, V);
  for (size_t i = 0; i < V; i++) triplets.emplace_back(i, i, 1);
  int j = 0;
  for (Face f : mesh.faces()) {
    Eigen::VectorXd weights = virtualRefinementAreaWeights[f];
    int i = 0;
    for (Vertex v : f.adjacentVertices()) {
      size_t vIdx = vertexIndices[v];
      triplets.emplace_back(V + j, vIdx, weights[i]);
      i++;
    }
    j++;
  }
  P.setFromTriplets(triplets.begin(), triplets.end());
  return P;
}

void EmbeddedGeometryInterface::computeVirtualRefinementAreaWeights() {
  vertexPositionsQ.ensureHave();

  virtualRefinementAreaWeights = FaceData<Eigen::VectorXd>(mesh);

  for (Face f : mesh.faces()) {
    Eigen::MatrixXd poly = getPolygonPositionMatrix(f);
    Eigen::VectorXd weights = computeVirtualVertex(poly);
    virtualRefinementAreaWeights[f] = weights;
  }
}

Eigen::MatrixXd EmbeddedGeometryInterface::getPolygonPositionMatrix(const Face& f) const {
  vertexPositionsQ.ensureHave();

  Eigen::MatrixXd poly(f.degree(), 3);
  int i = 0;
  for (Vertex v : f.adjacentVertices()) {
    for (int j = 0; j < 3; j++) {
      poly(i, j) = vertexPositions[v][j];
    }
    i++;
  }
  return poly;
}

Eigen::VectorXd EmbeddedGeometryInterface::computeVirtualVertex(const Eigen::MatrixXd& poly) const {

  // Given a polygon face, computes the affine weights that determine the position of the virtual vertex that minimizes
  // the sum of the squared areas of the triangles in the induced triangle fan. While the location of this vertex (the
  // minimizer) is unique, its expression as an affine combination of the polygon verties may not be -- regularize by
  // picking the weights with minimum L_2 norm, which encourages the weights to be as uniform as possible.

  int n = poly.rows();
  Eigen::VectorXd weights(n);
  Eigen::MatrixXd J(n, n);
  Eigen::VectorXd b(n);
  for (int i = 0; i < n; i++) {
    Eigen::Vector3d pk = poly.row(i);

    double Bk1_d2 = 0.0;
    double Bk1_d1 = 0.0;

    double Bk2_d0 = 0.0;
    double Bk2_d2 = 0.0;

    double Bk3_d0 = 0.0;
    double Bk3_d1 = 0.0;

    double CBk = 0.0;
    Eigen::Vector3d d = Eigen::MatrixXd::Zero(3, 1);

    for (int j = 0; j < n; j++) {
      Eigen::Vector3d pi = poly.row(j);
      Eigen::Vector3d pj = poly.row((j + 1) % n);
      d = pi - pj;

      double Bik1 = d(1) * pk(2) - d(2) * pk(1);
      double Bik2 = d(2) * pk(0) - d(0) * pk(2);
      double Bik3 = d(0) * pk(1) - d(1) * pk(0);

      double Ci1 = d(1) * pi(2) - d(2) * pi(1);
      double Ci2 = d(2) * pi(0) - d(0) * pi(2);
      double Ci3 = d(0) * pi(1) - d(1) * pi(0);

      Bk1_d1 += d(1) * Bik1;
      Bk1_d2 += d(2) * Bik1;

      Bk2_d0 += d(0) * Bik2;
      Bk2_d2 += d(2) * Bik2;

      Bk3_d0 += d(0) * Bik3;
      Bk3_d1 += d(1) * Bik3;

      CBk += Ci1 * Bik1 + Ci2 * Bik2 + Ci3 * Bik3;
    }
    for (int k = 0; k < n; k++) {
      Eigen::Vector3d xj = poly.row(k);
      J(i, k) =
          0.5 * (xj(2) * Bk1_d1 - xj(1) * Bk1_d2 + xj(0) * Bk2_d2 - xj(2) * Bk2_d0 + xj(1) * Bk3_d0 - xj(0) * Bk3_d1);
    }
    b(i) = 0.5 * CBk;
  }

  Eigen::MatrixXd M(n + 1, n);
  M.block(0, 0, n, n) = 4 * J;
  M.block(n, 0, 1, n).setOnes();

  Eigen::VectorXd b_(n + 1);
  b_.block(0, 0, n, 1) = 4 * b;

  b_(n) = 1.;
  weights = M.completeOrthogonalDecomposition().solve(b_).topRows(n);

  return weights;
}

Vector3 EmbeddedGeometryInterface::gradientHatFunction(const Vector3& a, const Vector3& b, const Vector3& c) const {
  Vector3 gradient;
  Vector3 site = a - b;
  Vector3 base = c - b;
  double area = 0.5 * (cross(site, base)).norm();
  double baseNorm = base.norm();
  Vector3 grad = site - (dot(site, base) / baseNorm) * base / baseNorm;
  if (area < 1e-10) {
    gradient = Vector3::Zero;
  } else {
    grad = baseNorm * grad / grad.norm();
    gradient = grad / (2.0 * area);
  }
  return gradient;
}


// = de Goes et al. "Discrete Differential Operators on Polygonal Meshes" (2020), based on the virtual element method.

// Laplacian
void EmbeddedGeometryInterface::computeVirtualElementLaplacian() {
  // TODO
}
void EmbeddedGeometryInterface::requireVirtualElementLaplacian() { virtualElementLaplacianQ.require(); }
void EmbeddedGeometryInterface::unrequireVirtualElementLaplacian() { virtualElementLaplacianQ.unrequire(); }


// Vertex Galerkin mass matrix (unlumped)
void EmbeddedGeometryInterface::computeVirtualElementVertexGalerkinMassMatrix() {
  // TODO
}
void EmbeddedGeometryInterface::requireVirtualElementVertexGalerkinMassMatrix() {
  virtualElementVertexGalerkinMassMatrixQ.require();
}
void EmbeddedGeometryInterface::unrequireVirtualElementVertexGalerkinMassMatrix() {
  virtualElementVertexGalerkinMassMatrixQ.unrequire();
}


// Vertex mass matrix (lumped)
void EmbeddedGeometryInterface::computeVirtualElementVertexLumpedMassMatrix() {
  // TODO
}
void EmbeddedGeometryInterface::requireVirtualElementVertexLumpedMassMatrix() {
  virtualElementVertexLumpedMassMatrixQ.require();
}
void EmbeddedGeometryInterface::unrequireVirtualElementVertexLumpedMassMatrix() {
  virtualElementVertexLumpedMassMatrixQ.unrequire();
}


// vertex connection Laplacian
void EmbeddedGeometryInterface::computeVirtualElementVertexConnectionLaplacian() {
  // TODO
}
void EmbeddedGeometryInterface::requireVirtualElementVertexConnectionLaplacian() {
  virtualElementVertexConnectionLaplacianQ.require();
}
void EmbeddedGeometryInterface::unrequireVirtualElementVertexConnectionLaplacian() {
  virtualElementVertexConnectionLaplacianQ.unrequire();
}


// DEC Operators
void EmbeddedGeometryInterface::computeVirtualElementDECOperators() {
  // TODO
}
void EmbeddedGeometryInterface::requireVirtualElementDECOperators() { virtualElementDECOperatorsQ.require(); }
void EmbeddedGeometryInterface::unrequireVirtualElementDECOperators() { virtualElementDECOperatorsQ.unrequire(); }

} // namespace surface
} // namespace geometrycentral
