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

  simplePolygonLaplacianQ                (&simplePolygonLaplacian,                std::bind(&EmbeddedGeometryInterface::computeSimplePolygonLaplacian, this),                quantities),
  simplePolygonVertexGalerkinMassMatrixQ (&simplePolygonVertexGalerkinMassMatrix, std::bind(&EmbeddedGeometryInterface::computeSimplePolygonVertexGalerkinMassMatrix, this), quantities),
  simplePolygonVertexLumpedMassMatrixQ   (&simplePolygonVertexLumpedMassMatrix,   std::bind(&EmbeddedGeometryInterface::computeSimplePolygonVertexLumpedMassMatrix, this),   quantities),
  virtualRefinementAreaWeightsQ          (&virtualRefinementAreaWeights,          std::bind(&EmbeddedGeometryInterface::computeVirtualRefinementAreaWeights, this),          quantities),
 
  polygonLaplacianQ                 (&polygonLaplacian,                 std::bind(&EmbeddedGeometryInterface::computePolygonLaplacian, this),                 quantities),
  polygonGradientMatrixQ            (&polygonGradientMatrix,            std::bind(&EmbeddedGeometryInterface::computePolygonGradientMatrix, this),            quantities),
  polygonDivergenceMatrixQ          (&polygonDivergenceMatrix,          std::bind(&EmbeddedGeometryInterface::computePolygonDivergenceMatrix, this),          quantities),
  polygonVertexLumpedMassMatrixQ    (&polygonVertexLumpedMassMatrix,    std::bind(&EmbeddedGeometryInterface::computePolygonVertexLumpedMassMatrix, this),    quantities),
  polygonVertexConnectionLaplacianQ (&polygonVertexConnectionLaplacian, std::bind(&EmbeddedGeometryInterface::computePolygonVertexConnectionLaplacian, this), quantities),
  
  polygonDECOperatorArray{&polygonHodge0, &polygonHodge0Inverse, &polygonHodge1, &polygonHodge2, &polygonHodge2Inverse, &polygonD0, &polygonD1},
  polygonDECOperatorsQ(&polygonDECOperatorArray, std::bind(&EmbeddedGeometryInterface::computePolygonDECOperators, this), quantities),
  polygonVertexNormalsQ             (&polygonVertexNormals,             std::bind(&EmbeddedGeometryInterface::computePolygonVertexNormals, this),             quantities)
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
void EmbeddedGeometryInterface::computeSimplePolygonLaplacian() {
  vertexIndicesQ.ensureHave();

  size_t V = mesh.nVertices();
  simplePolygonLaplacian = Eigen::SparseMatrix<double>(V, V);
  std::vector<Eigen::Triplet<double>> triplets;
  std::vector<size_t> vIndices; // indices of vertices of polygon face
  Eigen::MatrixXd Si;           // local per-polygon matrix
  for (Face f : mesh.faces()) {
    vIndices.clear();
    for (Vertex v : f.adjacentVertices()) vIndices.push_back(vertexIndices[v]);
    size_t n = f.degree();
    // Get local stiffness matrix.
    Si = simplePolygonStiffnessMatrix(f);
    // Add contribution to global mass matrix.
    for (size_t j = 0; j < n; j++) {
      for (size_t i = 0; i < n; i++) {
        triplets.emplace_back(vIndices[i], vIndices[j], Si(i, j));
      }
    }
  }
  simplePolygonLaplacian.setFromTriplets(triplets.begin(), triplets.end());
}
void EmbeddedGeometryInterface::requireSimplePolygonLaplacian() { simplePolygonLaplacianQ.require(); }
void EmbeddedGeometryInterface::unrequireSimplePolygonLaplacian() { simplePolygonLaplacianQ.unrequire(); }


// Vertex Galerkin mass matrix (unlumped)
void EmbeddedGeometryInterface::computeSimplePolygonVertexGalerkinMassMatrix() {
  vertexIndicesQ.ensureHave();

  size_t V = mesh.nVertices();
  simplePolygonVertexGalerkinMassMatrix = Eigen::SparseMatrix<double>(V, V);
  std::vector<Eigen::Triplet<double>> triplets;
  std::vector<size_t> vIndices; // indices of vertices of polygon face
  Eigen::MatrixXd Mi;           // local per-polygon matrix
  for (Face f : mesh.faces()) {
    vIndices.clear();
    for (Vertex v : f.adjacentVertices()) vIndices.push_back(vertexIndices[v]);
    size_t n = f.degree();
    // Get local mass matrix.
    Mi = simplePolygonMassMatrix(f);
    // Add contribution to global mass matrix.
    for (size_t j = 0; j < n; j++) {
      for (size_t i = 0; i < n; i++) {
        triplets.emplace_back(vIndices[i], vIndices[j], Mi(i, j));
      }
    }
  }
  simplePolygonVertexGalerkinMassMatrix.setFromTriplets(triplets.begin(), triplets.end());
}
void EmbeddedGeometryInterface::requireSimplePolygonVertexGalerkinMassMatrix() {
  simplePolygonVertexGalerkinMassMatrixQ.require();
}
void EmbeddedGeometryInterface::unrequireSimplePolygonVertexGalerkinMassMatrix() {
  simplePolygonVertexGalerkinMassMatrixQ.unrequire();
}


// Vertex mass matrix (lumped)
void EmbeddedGeometryInterface::computeSimplePolygonVertexLumpedMassMatrix() {
  vertexIndicesQ.ensureHave();

  size_t V = mesh.nVertices();
  simplePolygonVertexLumpedMassMatrix = Eigen::SparseMatrix<double>(V, V);
  std::vector<Eigen::Triplet<double>> triplets;
  std::vector<size_t> vIndices; // indices of vertices of polygon face
  Eigen::MatrixXd Mi;           // local per-polygon matrix
  for (Face f : mesh.faces()) {
    vIndices.clear();
    for (Vertex v : f.adjacentVertices()) vIndices.push_back(vertexIndices[v]);
    size_t n = f.degree();
    // Get local mass matrix.
    Mi = simplePolygonMassMatrix(f);
    // Add contribution to global mass matrix.
    for (size_t j = 0; j < n; j++) {
      for (size_t i = 0; i < n; i++) {
        triplets.emplace_back(vIndices[i], vIndices[i], Mi(i, j));
      }
    }
  }
  simplePolygonVertexLumpedMassMatrix.setFromTriplets(triplets.begin(), triplets.end());
}
void EmbeddedGeometryInterface::requireSimplePolygonVertexLumpedMassMatrix() {
  simplePolygonVertexLumpedMassMatrixQ.require();
}
void EmbeddedGeometryInterface::unrequireSimplePolygonVertexLumpedMassMatrix() {
  simplePolygonVertexLumpedMassMatrixQ.unrequire();
}

// Helper functions

Eigen::MatrixXd EmbeddedGeometryInterface::simplePolygonMassMatrix(const Face& f) {
  virtualRefinementAreaWeightsQ.ensureHave();

  size_t n = f.degree();
  Eigen::MatrixXd poly = polygonPositionMatrix(f);
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(n, n);
  const Eigen::VectorXd& weights = virtualRefinementAreaWeights[f];
  Eigen::Vector3d virtualVertex = poly.transpose() * weights;
  Eigen::VectorXd ln = Eigen::VectorXd::Zero(n + 1);
  double l[3], l2[3]; // lengths, lengths squared
  // Build triangle fan mass and cotan matrices
  for (size_t i = 0; i < n; i++) {
    const size_t i1 = (i + 1) % n;
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
  for (size_t j = 0; j < n; ++j)
    for (size_t i = 0; i < n; ++i) M(i, j) += weights(i) * ln(j) + weights(j) * ln(i) + weights(i) * weights(j) * ln(n);

  return M;
}

Eigen::MatrixXd EmbeddedGeometryInterface::simplePolygonStiffnessMatrix(const Face& f) {
  virtualRefinementAreaWeightsQ.ensureHave();

  size_t n = f.degree();
  Eigen::MatrixXd poly = polygonPositionMatrix(f);
  Eigen::MatrixXd S = Eigen::MatrixXd::Zero(n, n);
  const Eigen::VectorXd& weights = virtualRefinementAreaWeights[f];
  Eigen::Vector3d virtualVertex = poly.transpose() * weights;
  Eigen::VectorXd ln = Eigen::VectorXd::Zero(n + 1);
  double l[3], l2[3]; // lengths, lengths squared
  // Build triangle fan mass and cotan matrices
  for (size_t i = 0; i < n; i++) {
    const size_t i1 = (i + 1) % n;
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
  for (size_t j = 0; j < n; ++j)
    for (size_t i = 0; i < n; ++i) S(i, j) += weights(i) * ln(j) + weights(j) * ln(i) + weights(i) * weights(j) * ln(n);

  return S;
}

SparseMatrix<double> EmbeddedGeometryInterface::simplePolygonProlongationMatrix() {
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
    Eigen::MatrixXd poly = polygonPositionMatrix(f);
    Eigen::VectorXd weights = simplePolygonVirtualVertex(poly);
    virtualRefinementAreaWeights[f] = weights;
  }
}

Eigen::MatrixXd EmbeddedGeometryInterface::polygonPositionMatrix(const Face& f) {
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

Eigen::VectorXd EmbeddedGeometryInterface::simplePolygonVirtualVertex(const Eigen::MatrixXd& poly) const {

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


// = de Goes et al. "Discrete Differential Operators on Polygonal Meshes" (2020), based on the virtual element method.

// Laplacian
void EmbeddedGeometryInterface::computePolygonLaplacian() {
  vertexIndicesQ.ensureHave();

  size_t V = mesh.nVertices();
  polygonLaplacian = Eigen::SparseMatrix<double>(V, V);
  std::vector<Eigen::Triplet<double>> triplets;
  // Assemble per-face operators.
  Eigen::MatrixXd Lf;           // local per-polygon matrix
  std::vector<size_t> vIndices; // indices of vertices of polygon face
  for (Face f : mesh.faces()) {
    vIndices.clear();
    for (Vertex v : f.adjacentVertices()) vIndices.push_back(vertexIndices[v]);
    size_t n = f.degree();
    // Add contribution to global matrix.
    Lf = polygonPerFaceLaplacian(f);
    for (size_t j = 0; j < n; j++) {
      for (size_t i = 0; i < n; i++) {
        triplets.emplace_back(vIndices[i], vIndices[j], Lf(i, j));
      }
    }
  }
  polygonLaplacian.setFromTriplets(triplets.begin(), triplets.end());
}
void EmbeddedGeometryInterface::requirePolygonLaplacian() { polygonLaplacianQ.require(); }
void EmbeddedGeometryInterface::unrequirePolygonLaplacian() { polygonLaplacianQ.unrequire(); }


// Gradient matrix
void EmbeddedGeometryInterface::computePolygonGradientMatrix() {
  vertexIndicesQ.ensureHave();
  faceIndicesQ.ensureHave();

  size_t V = mesh.nVertices();
  size_t F = mesh.nFaces();
  polygonGradientMatrix = Eigen::SparseMatrix<double>(3 * F, V);
  std::vector<Eigen::Triplet<double>> triplets;
  // Assemble per-face operators.
  Eigen::MatrixXd Gf;           // local per-polygon matrix
  std::vector<size_t> vIndices; // indices of vertices of polygon face
  for (Face f : mesh.faces()) {
    size_t fIdx = faceIndices[f];
    vIndices.clear();
    for (Vertex v : f.adjacentVertices()) vIndices.push_back(vertexIndices[v]);
    size_t n = f.degree();
    // Add contribution to global matrix.
    Gf = polygonPerFaceGradientMatrix(f); // 3 x nf
    for (size_t i = 0; i < 3; i++) {
      for (size_t j = 0; j < n; j++) {
        triplets.emplace_back(3 * fIdx + i, vIndices[j], Gf(i, j));
      }
    }
  }
  polygonGradientMatrix.setFromTriplets(triplets.begin(), triplets.end());
}
void EmbeddedGeometryInterface::requirePolygonGradientMatrix() { polygonGradientMatrixQ.require(); }
void EmbeddedGeometryInterface::unrequirePolygonGradientMatrix() { polygonGradientMatrixQ.unrequire(); }


// Gradient matrix
void EmbeddedGeometryInterface::computePolygonDivergenceMatrix() {
  vertexIndicesQ.ensureHave();
  faceIndicesQ.ensureHave();

  size_t V = mesh.nVertices();
  size_t F = mesh.nFaces();
  polygonDivergenceMatrix = Eigen::SparseMatrix<double>(V, 3 * F);
  std::vector<Eigen::Triplet<double>> triplets;
  // Assemble per-face operators.
  Eigen::MatrixXd Gf;           // local per-polygon matrix
  std::vector<size_t> vIndices; // indices of vertices of polygon face
  for (Face f : mesh.faces()) {
    size_t fIdx = faceIndices[f];
    double area = polygonArea(f);
    vIndices.clear();
    for (Vertex v : f.adjacentVertices()) vIndices.push_back(vertexIndices[v]);
    size_t n = f.degree();
    // Add contribution to global matrix.
    Gf = polygonPerFaceGradientMatrix(f); // 3 x nf
    for (size_t i = 0; i < 3; i++) {
      for (size_t j = 0; j < n; j++) {
        triplets.emplace_back(vIndices[j], 3 * fIdx + i, Gf(i, j) * area);
      }
    }
  }
  polygonDivergenceMatrix.setFromTriplets(triplets.begin(), triplets.end());
}
void EmbeddedGeometryInterface::requirePolygonDivergenceMatrix() { polygonDivergenceMatrixQ.require(); }
void EmbeddedGeometryInterface::unrequirePolygonDivergenceMatrix() { polygonDivergenceMatrixQ.unrequire(); }


// Vertex mass matrix (lumped)
void EmbeddedGeometryInterface::computePolygonVertexLumpedMassMatrix() {
  vertexIndicesQ.ensureHave();

  size_t V = mesh.nVertices();
  Eigen::VectorXd hodge0(V);
  for (Vertex v : mesh.vertices()) {
    size_t vIdx = vertexIndices[v];
    double area = 0.;
    for (Face f : v.adjacentFaces()) {
      area += polygonArea(f) / f.degree();
    }
    hodge0[vIdx] = area;
  }
  polygonVertexLumpedMassMatrix = hodge0.asDiagonal();
}
void EmbeddedGeometryInterface::requirePolygonVertexLumpedMassMatrix() { polygonVertexLumpedMassMatrixQ.require(); }
void EmbeddedGeometryInterface::unrequirePolygonVertexLumpedMassMatrix() { polygonVertexLumpedMassMatrixQ.unrequire(); }


// vertex connection Laplacian
void EmbeddedGeometryInterface::computePolygonVertexConnectionLaplacian() {
  vertexIndicesQ.ensureHave();

  size_t V = mesh.nVertices();
  polygonVertexConnectionLaplacian = Eigen::SparseMatrix<std::complex<double>>(V, V);
  std::vector<Eigen::Triplet<std::complex<double>>> triplets;
  Eigen::MatrixXd Lf;
  std::vector<size_t> vIndices;
  for (Face f : mesh.faces()) {
    vIndices.clear();
    for (Vertex v : f.adjacentVertices()) vIndices.push_back(vertexIndices[v]);
    size_t n = f.degree();
    Lf = polygonPerFaceConnectionLaplacian(f);
    for (size_t j = 0; j < n; j++) {
      for (size_t i = 0; i < n; i++) {
        double re = Lf(2 * i, 2 * j);     // == Lf(2 * i + 1, 2 * j + 1)
        double im = Lf(2 * i, 2 * j + 1); // == -Lf(2 * i + 1, 2 * j)
        GC_SAFETY_ASSERT(std::abs(Lf(2 * i, 2 * j) - Lf(2 * i + 1, 2 * j + 1)) < 1e-5,
                         "polygon vertex connection Laplacian error");
        GC_SAFETY_ASSERT(std::abs(Lf(2 * i, 2 * j + 1) + Lf(2 * i + 1, 2 * j)) < 1e-5,
                         "polygon vertex connection Laplacian error");
        triplets.emplace_back(vIndices[i], vIndices[j], std::complex<double>(re, im));
      }
    }
  }
  polygonVertexConnectionLaplacian.setFromTriplets(triplets.begin(), triplets.end());
}
void EmbeddedGeometryInterface::requirePolygonVertexConnectionLaplacian() {
  polygonVertexConnectionLaplacianQ.require();
}
void EmbeddedGeometryInterface::unrequirePolygonVertexConnectionLaplacian() {
  polygonVertexConnectionLaplacianQ.unrequire();
}


void EmbeddedGeometryInterface::computePolygonDECOperators() {
  vertexIndicesQ.ensureHave();
  edgeIndicesQ.ensureHave();
  faceIndicesQ.ensureHave();

  size_t V = mesh.nVertices();
  size_t E = mesh.nEdges();
  size_t F = mesh.nFaces();
  std::vector<Eigen::Triplet<double>> tripletsD0, tripletsD1, tripletsH1;

  // exterior derivatives the same as the simplicial case
  polygonD0 = Eigen::SparseMatrix<double>(E, V);
  for (Edge e : mesh.edges()) {
    size_t eIdx = edgeIndices[e];
    size_t vA = vertexIndices[e.firstVertex()];
    size_t vB = vertexIndices[e.secondVertex()];
    tripletsD0.emplace_back(eIdx, vA, -1.);
    tripletsD0.emplace_back(eIdx, vB, 1.);
  }
  polygonD0.setFromTriplets(tripletsD0.begin(), tripletsD0.end());

  polygonD1 = Eigen::SparseMatrix<double>(F, E);
  for (Face f : mesh.faces()) {
    size_t fIdx = faceIndices[f];
    for (Halfedge he : f.adjacentHalfedges()) {
      size_t eIdx = edgeIndices[he.edge()];
      double sign = he.orientation() ? 1. : -1.;
      tripletsD1.emplace_back(fIdx, eIdx, sign);
    }
  }
  polygonD1.setFromTriplets(tripletsD1.begin(), tripletsD1.end());

  // hodge0
  Eigen::VectorXd h0(V);
  for (Vertex v : mesh.vertices()) {
    size_t vIdx = vertexIndices[v];
    double area = 0.;
    for (Face f : v.adjacentFaces()) {
      area += polygonArea(f) / f.degree();
    }
    h0[vIdx] = area;
  }
  polygonHodge0 = h0.asDiagonal();
  polygonHodge0Inverse = h0.asDiagonal().inverse();

  // hodge1 (inner product matrix on 1-forms)
  polygonHodge1 = Eigen::SparseMatrix<double>(E, E);
  Eigen::MatrixXd Mf;           // local per-polygon matrix
  std::vector<size_t> eIndices; // indices of edges of polygon face
  for (Face f : mesh.faces()) {
    eIndices.clear();
    // Assumes edge iterator goes in same order as vertex iterator around a face
    for (Edge e : f.adjacentEdges()) eIndices.push_back(edgeIndices[e]);
    size_t n = f.degree();
    // Add contribution to global mass matrix.
    Mf = polygonPerFaceInnerProductMatrix(f);
    for (size_t j = 0; j < n; j++) {
      for (size_t i = 0; i < n; i++) {
        tripletsH1.emplace_back(eIndices[i], eIndices[j], Mf(i, j));
      }
    }
  }
  polygonHodge1.setFromTriplets(tripletsH1.begin(), tripletsH1.end());

  // hodge2
  Eigen::VectorXd h2(F);
  for (Face f : mesh.faces()) {
    size_t fIdx = faceIndices[f];
    h2[fIdx] = 1. / polygonArea(f);
  }
  polygonHodge2 = h2.asDiagonal();
  polygonHodge2Inverse = h2.asDiagonal().inverse();
}
void EmbeddedGeometryInterface::requirePolygonDECOperators() { polygonDECOperatorsQ.require(); }
void EmbeddedGeometryInterface::unrequirePolygonDECOperators() { polygonDECOperatorsQ.unrequire(); }


// Helper functions

Eigen::MatrixXd EmbeddedGeometryInterface::polygonPerFaceLaplacian(const Face& f) {
  Eigen::MatrixXd Df = polygonDerivativeMatrix(f);
  return Df.transpose() * polygonPerFaceInnerProductMatrix(f) * Df; // build positive-definite
}

Eigen::MatrixXd EmbeddedGeometryInterface::polygonPerFaceInnerProductMatrix(const Face& f) {
  Eigen::MatrixXd Uf = polygonSharp(f);
  Eigen::MatrixXd Pf = polygonProjectionMatrix(f);
  double A = polygonArea(f);
  return A * Uf.transpose() * Uf + polygonLambda * Pf.transpose() * Pf;
}

Eigen::MatrixXd EmbeddedGeometryInterface::polygonPerFaceConnectionLaplacian(const Face& f) {
  Eigen::MatrixXd G = polygonCovariantGradient(f);
  Eigen::MatrixXd P = polygonCovariantProjection(f);
  double A = polygonArea(f);
  Eigen::MatrixXd L = A * G.transpose() * G + polygonLambda * P.transpose() * P; // build positive-definite
  return L;
}

Eigen::MatrixXd EmbeddedGeometryInterface::polygonBlockConnection(const Face& f) {
  size_t d = f.degree();
  Eigen::MatrixXd R = Eigen::MatrixXd::Zero(2 * d, 2 * d);
  size_t cpt = 0;
  for (Vertex v : f.adjacentVertices()) {
    Eigen::Matrix2d Rv = Rvf(v, f);
    R.block<2, 2>(2 * cpt, 2 * cpt) = Rv;
    cpt++;
  }
  return R;
}

Eigen::MatrixXd EmbeddedGeometryInterface::polygonCovariantGradient(const Face& f) {
  return kroneckerWithI2(Tf(f).transpose() * polygonPerFaceGradientMatrix(f)) * polygonBlockConnection(f);
}

Eigen::MatrixXd EmbeddedGeometryInterface::polygonCovariantProjection(const Face& f) {
  return kroneckerWithI2(polygonProjectionMatrix(f) * polygonDerivativeMatrix(f)) * polygonBlockConnection(f);
}

Eigen::MatrixXd EmbeddedGeometryInterface::polygonProjectionMatrix(const Face& f) {
  size_t d = f.degree();
  Eigen::MatrixXd P = Eigen::MatrixXd::Identity(d, d) - polygonFlat(f) * polygonSharp(f);
  return P;
}

Eigen::MatrixXd EmbeddedGeometryInterface::polygonPerFaceGradientMatrix(const Face& f) {
  // equivalent to applying the (per-face) sharp operator to the (per-face) exterior derivative
  Eigen::Vector3d An = polygonVectorArea(f);
  double A = An.norm();
  Eigen::Vector3d N = An / A;
  return -1. / A * bracket(N) * polygonCoGradientMatrix(f);
}

Eigen::MatrixXd EmbeddedGeometryInterface::polygonCoGradientMatrix(const Face& f) {
  return polygonEdgeVectorMatrix(f).transpose() * polygonAveragingMatrix(f);
}

Eigen::MatrixXd EmbeddedGeometryInterface::polygonAveragingMatrix(const Face& f) const {
  size_t d = f.degree();
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(d, d);
  for (size_t i = 0; i < d; i++) {
    A(i, (i + 1) % d) = 0.5;
    A(i, i) = 0.5;
  }
  return A;
}

Eigen::MatrixXd EmbeddedGeometryInterface::polygonDerivativeMatrix(const Face& f) const {
  size_t d = f.degree();
  Eigen::MatrixXd D = Eigen::MatrixXd::Zero(d, d);
  for (size_t i = 0; i < d; i++) {
    D(i, (i + 1) % d) = 1.;
    D(i, i) = -1.;
  }
  return D;
}

Eigen::MatrixXd EmbeddedGeometryInterface::polygonEdgeVectorMatrix(const Face& f) {
  return polygonDerivativeMatrix(f) * polygonPositionMatrix(f);
}

Eigen::MatrixXd EmbeddedGeometryInterface::polygonEdgeMidpointMatrix(const Face& f) {
  return polygonAveragingMatrix(f) * polygonPositionMatrix(f);
}

Eigen::MatrixXd EmbeddedGeometryInterface::polygonFlat(const Face& f) {

  Eigen::Vector3d N = polygonNormal(f);
  return polygonEdgeVectorMatrix(f) * (Eigen::MatrixXd::Identity(3, 3) - N * N.transpose());
}

Eigen::MatrixXd EmbeddedGeometryInterface::polygonSharp(const Face& f) {

  size_t d = f.degree();
  Eigen::Vector3d An = polygonVectorArea(f);
  double A = An.norm();
  Eigen::Vector3d N = An / A;
  Eigen::Vector3d c = polygonCentroid(f);
  return 1. / A * bracket(N) * (polygonEdgeMidpointMatrix(f).transpose() - c * Eigen::VectorXd::Ones(d).transpose());
}

void EmbeddedGeometryInterface::computePolygonVertexNormals() {

  polygonVertexNormals = VertexData<Eigen::VectorXd>(mesh);
  for (Vertex v : mesh.vertices()) {
    double totalArea = 0.;
    Eigen::Vector3d vN(0., 0., 0.);
    for (Face f : v.adjacentFaces()) {
      Eigen::Vector3d N = polygonVectorArea(f);
      vN += N;
      totalArea += N.norm();
    }
    vN /= totalArea;
    polygonVertexNormals[v] = vN;
  }
}

Eigen::Vector3d EmbeddedGeometryInterface::polygonVectorArea(const Face& f) {
  vertexPositionsQ.ensureHave();

  Vector3 N = {0, 0, 0};
  for (Halfedge he : f.adjacentHalfedges()) {
    Vertex vA = he.vertex();
    Vertex vB = he.next().vertex();
    Vector3 pA = vertexPositions[vA];
    Vector3 pB = vertexPositions[vB];
    N += cross(pA, pB);
  }
  N *= 0.5;
  return Eigen::Vector3d(N[0], N[1], N[2]);
}

double EmbeddedGeometryInterface::polygonArea(const Face& f) {
  Eigen::Vector3d An = polygonVectorArea(f);
  double A = An.norm();
  return A;
}

Eigen::Vector3d EmbeddedGeometryInterface::polygonNormal(const Face& f) {
  Eigen::Vector3d N = polygonVectorArea(f);
  N /= N.norm();
  return N;
}

Eigen::Vector3d EmbeddedGeometryInterface::polygonCentroid(const Face& f) {
  vertexPositionsQ.ensureHave();

  Vector3 c = {0, 0, 0};
  for (Vertex v : f.adjacentVertices()) {
    c += vertexPositions[v];
  }
  c /= f.degree();
  return Eigen::Vector3d(c[0], c[1], c[2]);
}

Eigen::MatrixXd EmbeddedGeometryInterface::Tv(const Vertex& v) {
  polygonVertexNormalsQ.ensureHave();

  // Return 3 x 2 matrix defining the tangent space at vertex v, with basis vectors in columns.
  Eigen::Vector3d nv = polygonVertexNormals[v];
  Vector3 pA = vertexPositions[v.halfedge().vertex()];
  Vector3 pB = vertexPositions[v.halfedge().next().vertex()];
  Vector3 heVec = pB - pA;
  Eigen::Vector3d w(heVec[0], heVec[1], heVec[2]);
  Eigen::Vector3d uu = project(w, nv).normalized();
  Eigen::Vector3d vv = nv.cross(uu);
  Eigen::MatrixXd B(3, 2);
  B.col(0) = uu;
  B.col(1) = vv;
  return B;
}

Eigen::MatrixXd EmbeddedGeometryInterface::Tf(const Face& f) {
  vertexPositionsQ.ensureHave();

  // Return 3 x 2 matrix defining the tangent space at face f, with basis vectors in columns.
  Eigen::Vector3d nf = polygonNormal(f);
  Vector3 pA = vertexPositions[f.halfedge().vertex()];
  Vector3 pB = vertexPositions[f.halfedge().next().vertex()];
  Vector3 heVec = pB - pA;
  Eigen::Vector3d w(heVec[0], heVec[1], heVec[2]);
  Eigen::Vector3d uu = project(w, nf).normalized();
  Eigen::Vector3d vv = nf.cross(uu);
  Eigen::MatrixXd B(3, 2);
  B.col(0) = uu;
  B.col(1) = vv;
  return B;
}

Eigen::Matrix2d EmbeddedGeometryInterface::Rvf(const Vertex& v, const Face& f) {
  return Tf(f).transpose() * Qvf(v, f) * Tv(v);
}

Eigen::Matrix3d EmbeddedGeometryInterface::Qvf(const Vertex& v, const Face& f) {
  polygonVertexNormalsQ.ensureHave();

  // Return 3 x 3 rotation matrix to align n_v to n_f.
  Eigen::Vector3d nf = polygonNormal(f);
  Eigen::Vector3d nv = polygonVertexNormals[v];
  double c = nv.dot(nf);

  // Special case for opposite nv and nf vectors.
  if (std::abs(c + 1.0) < 1e-5) return -Eigen::Matrix3d::Identity();

  Eigen::Vector3d vv = nv.cross(nf);
  Eigen::Matrix3d skew = bracket(vv);
  return Eigen::Matrix3d::Identity() + skew + 1.0 / (1.0 + c) * skew * skew;
}

Eigen::Matrix3d EmbeddedGeometryInterface::bracket(const Eigen::Vector3d& n) const {
  Eigen::Matrix3d B;
  B << 0., -n[2], n[1], n[2], 0., -n[0], -n[1], n[0], 0.;
  return B;
}

Eigen::Vector3d EmbeddedGeometryInterface::project(const Eigen::Vector3d& u, const Eigen::Vector3d& n) const {
  // Project u on the orthgonal of n.
  return u - (u.dot(n) / n.squaredNorm()) * n;
}

Eigen::MatrixXd EmbeddedGeometryInterface::kroneckerWithI2(const Eigen::MatrixXd& M) const {
  size_t h = M.rows();
  size_t w = M.cols();
  Eigen::MatrixXd MK = Eigen::MatrixXd::Zero(h * 2, w * 2);
  for (size_t j = 0; j < h; j++)
    for (size_t i = 0; i < w; i++) {
      MK(2 * j, 2 * i) = M(j, i);
      MK(2 * j + 1, 2 * i + 1) = M(j, i);
    }
  return MK;
}

} // namespace surface
} // namespace geometrycentral
