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

  simplePolygonLaplacianQ                 (&simplePolygonLaplacian,                 std::bind(&EmbeddedGeometryInterface::computeSimplePolygonLaplacian, this),                 quantities),
  simplePolygonDivergenceMatrixQ          (&simplePolygonDivergenceMatrix,          std::bind(&EmbeddedGeometryInterface::computeSimplePolygonDivergenceMatrix, this),          quantities),
  simplePolygonGradientMatrixQ            (&simplePolygonGradientMatrix,            std::bind(&EmbeddedGeometryInterface::computeSimplePolygonGradientMatrix, this),            quantities),
  simplePolygonProlongationMatrixQ        (&simplePolygonProlongationMatrix,        std::bind(&EmbeddedGeometryInterface::computeSimplePolygonProlongationMatrix, this),        quantities),
  simplePolygonVertexConnectionLaplacianQ (&simplePolygonVertexConnectionLaplacian, std::bind(&EmbeddedGeometryInterface::computeSimplePolygonVertexConnectionLaplacian, this), quantities),
  simplePolygonVertexGalerkinMassMatrixQ  (&simplePolygonVertexGalerkinMassMatrix,  std::bind(&EmbeddedGeometryInterface::computeSimplePolygonVertexGalerkinMassMatrix, this),  quantities),
  simplePolygonVertexLumpedMassMatrixQ    (&simplePolygonVertexLumpedMassMatrix,    std::bind(&EmbeddedGeometryInterface::computeSimplePolygonVertexLumpedMassMatrix, this),    quantities),
  virtualRefinementAreaWeightsQ           (&virtualRefinementAreaWeights,           std::bind(&EmbeddedGeometryInterface::computeVirtualRefinementAreaWeights, this),           quantities),
  virtualRefinementAreaPointsQ            (&virtualRefinementAreaPoints,            std::bind(&EmbeddedGeometryInterface::computeVirtualRefinementAreaWeights, this),           quantities)

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
    Vector3 N = {0, 0, 0};
    for (Halfedge he : f.adjacentHalfedges()) {
      Vertex vA = he.vertex();
      Vertex vB = he.next().vertex();
      Vector3 pA = vertexPositions[vA];
      Vector3 pB = vertexPositions[vB];
      N += cross(pA, pB);
    }
    double area = 0.5 * norm(N);
    faceAreas[f] = area;
  }
}

// Override to compute directly from vertex positions
void EmbeddedGeometryInterface::computeCornerAngles() {
  vertexPositionsQ.ensureHave();

  cornerAngles = CornerData<double>(mesh);

  for (Face f : mesh.faces()) {
    for (Halfedge he : f.adjacentHalfedges()) {
      // WARNING: Logic duplicated between cached and immediate version
      Vector3 pA = vertexPositions[he.vertex()];
      Halfedge heNext = he.next();
      Vector3 pB = vertexPositions[heNext.vertex()];
      Vector3 pC = vertexPositions[heNext.next().vertex()];

      double q = dot(unit(pC - pB), unit(pA - pB));
      q = clamp(q, -1.0, 1.0);
      double angle = std::acos(q);

      cornerAngles[heNext.corner()] = angle;
    }
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
// Copyright (C) 2020 Astrid Bunge, Philipp Herholz, Misha Kazhdan, Mario Botsch, MIT license
// (Modified to work in geometry-central. Original code can be found here: https://github.com/mbotsch/polygon-laplacian)

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

void EmbeddedGeometryInterface::computeSimplePolygonDivergenceMatrix() {

  simplePolygonGradientMatrixQ.ensureHave();

  SparseMatrix<double> M = simplePolygonGradientMassMatrix();
  simplePolygonDivergenceMatrix = -simplePolygonGradientMatrix.transpose() * M;
}
void EmbeddedGeometryInterface::requireSimplePolygonDivergenceMatrix() { simplePolygonDivergenceMatrixQ.require(); }
void EmbeddedGeometryInterface::unrequireSimplePolygonDivergenceMatrix() { simplePolygonDivergenceMatrixQ.unrequire(); }

// 3|T^f| x |V| gradient matrix
void EmbeddedGeometryInterface::computeSimplePolygonGradientMatrix() {

  faceIndicesQ.ensureHave();
  vertexIndicesQ.ensureHave();
  vertexPositionsQ.ensureHave();
  simplePolygonProlongationMatrixQ.ensureHave();

  size_t V = mesh.nVertices();
  size_t F = mesh.nFaces();
  SparseMatrix<double> G;
  std::vector<Eigen::Triplet<double>> triplets;
  int nTriangles = 0;
  int k = 0;
  for (Face f : mesh.faces()) {
    size_t fIdx = faceIndices[f];
    nTriangles += f.degree();
    Eigen::Vector3d p = virtualRefinementAreaPoints[f];
    for (Halfedge he : f.adjacentHalfedges()) {
      size_t v0 = vertexIndices[he.tailVertex()];
      size_t v1 = vertexIndices[he.tipVertex()];
      const Vector3& pos0 = vertexPositions[v0];
      const Vector3& pos1 = vertexPositions[v1];
      Eigen::Vector3d p0 = {pos0[0], pos0[1], pos0[2]};
      Eigen::Vector3d p1 = {pos1[0], pos1[1], pos1[2]};
      Eigen::Vector3d grad_p = gradientHatFunction(p, p0, p1);
      Eigen::Vector3d grad_p0 = gradientHatFunction(p0, p1, p);
      Eigen::Vector3d grad_p1 = gradientHatFunction(p1, p, p0);
      for (int j = 0; j < 3; j++) {
        triplets.emplace_back(3 * k + j, V + fIdx, grad_p(j));
        triplets.emplace_back(3 * k + j, v0, grad_p0(j));
        triplets.emplace_back(3 * k + j, v1, grad_p1(j));
      }
      k++;
    }
  }
  G.resize(3 * nTriangles, V + F);
  G.setFromTriplets(triplets.begin(), triplets.end());
  simplePolygonGradientMatrix = G * simplePolygonProlongationMatrix;
}
void EmbeddedGeometryInterface::requireSimplePolygonGradientMatrix() { simplePolygonGradientMatrixQ.require(); }
void EmbeddedGeometryInterface::unrequireSimplePolygonGradientMatrix() { simplePolygonGradientMatrixQ.unrequire(); }

// Vertex Connection Laplacian. Formed by augmenting a scalar Laplacian with rotations between each pair of vertices
// within each polygon. The rotations are based on extrinsic (esimated) tangent planes.
void EmbeddedGeometryInterface::computeSimplePolygonVertexConnectionLaplacian() {

  simplePolygonLaplacianQ.ensureHave();
  vertexNormalsQ.ensureHave();
  vertexTangentBasisQ.ensureHave();
  vertexIndicesQ.ensureHave();

  // Compute rotation between tangent plane at vertex i and at vertex j.
  // Rotations are not just along halfedges, but between any pair of vertices sharing a face.
  std::map<std::pair<size_t, size_t>, std::complex<double>> rotations; // entry at ij = rotation from i to j
  for (Face f : mesh.faces()) {
    for (Halfedge he : f.adjacentHalfedges()) {
      Vertex vi = he.tailVertex();
      Vertex vj = he.tipVertex();
      Vector3 ni = vertexNormals[vi];
      Vector3 nj = vertexNormals[vj];
      Vector3 axis = cross(ni, nj);
      Vector3 ei = vertexTangentBasis[vi][0];
      Vector3 ej = vertexTangentBasis[vj][0];
      Vector3 R_ei = ei;
      if (axis.norm() > 1e-8) {
        double angle = angleInPlane(ni, nj, axis);
        R_ei = ei.rotateAround(axis, angle);
      }
      std::complex<double> r_ij = std::complex<double>(Vector2::fromAngle(angleInPlane(R_ei, ej, nj)));

      size_t viIdx = vertexIndices[vi];
      size_t vjIdx = vertexIndices[vj];
      if (viIdx > vjIdx) {
        viIdx = vertexIndices[vj];
        vjIdx = vertexIndices[vi];
        r_ij = std::conj(r_ij);
      }
      std::pair<size_t, size_t> key = std::make_pair(viIdx, vjIdx);
      rotations[key] = r_ij;
    }
  }

  size_t V = mesh.nVertices();
  const SparseMatrix<double>& L = simplePolygonLaplacian;
  simplePolygonVertexConnectionLaplacian = Eigen::SparseMatrix<std::complex<double>>(V, V);
  std::vector<Eigen::Triplet<std::complex<double>>> triplets;
  for (int k = 0; k < L.outerSize(); ++k) {
    for (SparseMatrix<double>::InnerIterator it(L, k); it; ++it) {
      size_t i = it.row();
      size_t j = it.col();
      if (i == j) {
        triplets.emplace_back(it.row(), it.col(), it.value());
        continue;
      }
      bool sgn = (i < j);
      std::pair<size_t, size_t> key = sgn ? std::make_pair(i, j) : std::make_pair(j, i);
      std::complex<double> rot = sgn ? rotations[key] : std::conj(rotations[key]);
      triplets.emplace_back(it.row(), it.col(), it.value() * rot);
    }
  }
  simplePolygonVertexConnectionLaplacian.setFromTriplets(triplets.begin(), triplets.end());
}
void EmbeddedGeometryInterface::requireSimplePolygonVertexConnectionLaplacian() {
  simplePolygonVertexConnectionLaplacianQ.require();
}
void EmbeddedGeometryInterface::unrequireSimplePolygonVertexConnectionLaplacian() {
  simplePolygonVertexConnectionLaplacianQ.unrequire();
}

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

Eigen::VectorXd EmbeddedGeometryInterface::simplePolygonVirtualVertex(const Eigen::MatrixXd& poly) {

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

/*
 * Block diagonal matrix whose i-th block consists of the 3×3 identity matrix multiplied by the area of the i-th
 * triangle.
 */
SparseMatrix<double> EmbeddedGeometryInterface::simplePolygonGradientMassMatrix() {

  virtualRefinementAreaWeightsQ.ensureHave();
  vertexPositionsQ.ensureHave();

  SparseMatrix<double> M;
  std::vector<Eigen::Triplet<double>> triplets;
  int c = 0;
  for (Face f : mesh.faces()) {
    Eigen::Vector3d areaPoint = virtualRefinementAreaPoints[f];
    Vector3 ap = {areaPoint[0], areaPoint[1], areaPoint[2]};
    int i = 0;
    for (Halfedge he : f.adjacentHalfedges()) {
      Vector3 p0 = vertexPositions[he.tailVertex()];
      Vector3 p1 = vertexPositions[he.tipVertex()];
      double area = 0.5 * cross(p0 - ap, p1 - ap).norm();
      for (int j = 0; j < 3; j++) {
        int idx = c + 3 * i + j;
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

Eigen::MatrixXd EmbeddedGeometryInterface::simplePolygonMassMatrix(const Face& f) {
  virtualRefinementAreaWeightsQ.ensureHave();

  size_t n = f.degree();
  Eigen::MatrixXd poly = polygonPositionMatrix(f);
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(n, n);
  const Eigen::VectorXd& weights = virtualRefinementAreaWeights[f];
  Eigen::Vector3d virtualVertex = virtualRefinementAreaPoints[f];
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
  Eigen::Vector3d virtualVertex = virtualRefinementAreaPoints[f];
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

void EmbeddedGeometryInterface::computeSimplePolygonProlongationMatrix() {
  virtualRefinementAreaWeightsQ.ensureHave();
  vertexIndicesQ.ensureHave();

  size_t V = mesh.nVertices();
  size_t F = mesh.nFaces();
  std::vector<Eigen::Triplet<double>> triplets;
  simplePolygonProlongationMatrix.resize(V + F, V);
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
  simplePolygonProlongationMatrix.setFromTriplets(triplets.begin(), triplets.end());
}
void EmbeddedGeometryInterface::requireSimplePolygonProlongationMatrix() { simplePolygonProlongationMatrixQ.require(); }
void EmbeddedGeometryInterface::unrequireSimplePolygonProlongationMatrix() {
  simplePolygonProlongationMatrixQ.unrequire();
}

void EmbeddedGeometryInterface::computeVirtualRefinementAreaWeights() {
  vertexPositionsQ.ensureHave();

  virtualRefinementAreaWeights = FaceData<Eigen::VectorXd>(mesh);
  virtualRefinementAreaPoints = FaceData<Eigen::Vector3d>(mesh);

  for (Face f : mesh.faces()) {
    Eigen::MatrixXd poly = polygonPositionMatrix(f);
    Eigen::VectorXd weights = simplePolygonVirtualVertex(poly);
    virtualRefinementAreaPoints[f] = poly.transpose() * weights;
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

Eigen::Vector3d EmbeddedGeometryInterface::gradientHatFunction(const Eigen::Vector3d& a, const Eigen::Vector3d& b,
                                                               const Eigen::Vector3d& c) const {

  Eigen::Vector3d gradient;
  Eigen::Vector3d site = a - b;
  Eigen::Vector3d base = c - b;
  double area = 0.5 * (site.cross(base)).norm();
  double baseNorm = base.norm();
  Eigen::Vector3d grad = site - (site.dot(base) / baseNorm) * base / baseNorm;
  grad = baseNorm * grad / grad.norm();
  gradient = Eigen::Vector3d(grad[0], grad[1], grad[2]) / (2.0 * area);
  return gradient;
}

} // namespace surface
} // namespace geometrycentral
