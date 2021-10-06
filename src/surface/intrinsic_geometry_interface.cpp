#include "geometrycentral/surface/intrinsic_geometry_interface.h"

//#include "geometrycentral/surface/discrete_operators.h"

#include <fstream>
#include <limits>

namespace geometrycentral {
namespace surface {

// clang-format off
IntrinsicGeometryInterface::IntrinsicGeometryInterface(SurfaceMesh& mesh_) : 
  BaseGeometryInterface(mesh_), 

  edgeLengthsQ              (&edgeLengths,                  std::bind(&IntrinsicGeometryInterface::computeEdgeLengths, this),               quantities),
  faceAreasQ                (&faceAreas,                    std::bind(&IntrinsicGeometryInterface::computeFaceAreas, this),                 quantities),
  vertexDualAreasQ          (&vertexDualAreas,              std::bind(&IntrinsicGeometryInterface::computeVertexDualAreas, this),           quantities),
  cornerAnglesQ             (&cornerAngles,                 std::bind(&IntrinsicGeometryInterface::computeCornerAngles, this),              quantities),
  vertexAngleSumsQ          (&vertexAngleSums,              std::bind(&IntrinsicGeometryInterface::computeVertexAngleSums, this),           quantities),
  cornerScaledAnglesQ       (&cornerScaledAngles,           std::bind(&IntrinsicGeometryInterface::computeCornerScaledAngles, this),        quantities),
  vertexGaussianCurvaturesQ (&vertexGaussianCurvatures,     std::bind(&IntrinsicGeometryInterface::computeVertexGaussianCurvatures, this),  quantities),
  faceGaussianCurvaturesQ   (&faceGaussianCurvatures,       std::bind(&IntrinsicGeometryInterface::computeFaceGaussianCurvatures, this),    quantities),
  halfedgeCotanWeightsQ     (&halfedgeCotanWeights,         std::bind(&IntrinsicGeometryInterface::computeHalfedgeCotanWeights, this),      quantities),
  edgeCotanWeightsQ         (&edgeCotanWeights,             std::bind(&IntrinsicGeometryInterface::computeEdgeCotanWeights, this),          quantities),
  shapeLengthScaleQ         (&shapeLengthScale,             std::bind(&IntrinsicGeometryInterface::computeShapeLengthScale, this),          quantities),
  meshLengthScaleQ          (&meshLengthScale,              std::bind(&IntrinsicGeometryInterface::computeMeshLengthScale, this),          quantities),
  
  halfedgeVectorsInFaceQ            (&halfedgeVectorsInFace,            std::bind(&IntrinsicGeometryInterface::computeHalfedgeVectorsInFace, this),             quantities),
  transportVectorsAcrossHalfedgeQ   (&transportVectorsAcrossHalfedge,   std::bind(&IntrinsicGeometryInterface::computeTransportVectorsAcrossHalfedge, this),    quantities),
  halfedgeVectorsInVertexQ          (&halfedgeVectorsInVertex,          std::bind(&IntrinsicGeometryInterface::computeHalfedgeVectorsInVertex, this),           quantities),
  transportVectorsAlongHalfedgeQ    (&transportVectorsAlongHalfedge,    std::bind(&IntrinsicGeometryInterface::computeTransportVectorsAlongHalfedge, this),     quantities),

  cotanLaplacianQ               (&cotanLaplacian,               std::bind(&IntrinsicGeometryInterface::computeCotanLaplacian, this),                quantities),
  vertexLumpedMassMatrixQ       (&vertexLumpedMassMatrix,       std::bind(&IntrinsicGeometryInterface::computeVertexLumpedMassMatrix, this),        quantities),
  vertexGalerkinMassMatrixQ     (&vertexGalerkinMassMatrix,     std::bind(&IntrinsicGeometryInterface::computeVertexGalerkinMassMatrix, this),      quantities),
  vertexConnectionLaplacianQ    (&vertexConnectionLaplacian,    std::bind(&IntrinsicGeometryInterface::computeVertexConnectionLaplacian, this),     quantities),
  faceGalerkinMassMatrixQ       (&faceGalerkinMassMatrix,       std::bind(&IntrinsicGeometryInterface::computeFaceGalerkinMassMatrix, this),        quantities),
  faceConnectionLaplacianQ      (&faceConnectionLaplacian,      std::bind(&IntrinsicGeometryInterface::computeFaceConnectionLaplacian, this),       quantities),


  // DEC operators need some extra work since 8 members are grouped under one require
  DECOperatorArray{&hodge0, &hodge0Inverse, &hodge1, &hodge1Inverse, &hodge2, &hodge2Inverse, &d0, &d1},
  DECOperatorsQ(&DECOperatorArray, std::bind(&IntrinsicGeometryInterface::computeDECOperators, this), quantities)


  { }
// clang-format on

// === Quantity implementations

// Edge lengths
void IntrinsicGeometryInterface::requireEdgeLengths() { edgeLengthsQ.require(); }
void IntrinsicGeometryInterface::unrequireEdgeLengths() { edgeLengthsQ.unrequire(); }

void IntrinsicGeometryInterface::computeFaceAreas() {
  edgeLengthsQ.ensureHave();

  // ONEDAY try these for better accuracy in near-degenerate triangles?
  // "Miscalculating Area and Angles of a Needle-like Triangle" https://www.cs.unc.edu/~snoeyink/c/c205/Triangle.pdf

  faceAreas = FaceData<double>(mesh);
  for (Face f : mesh.faces()) {
    // WARNING: Logic duplicated between cached and immediate version

    Halfedge he = f.halfedge();
    double a = edgeLengths[he.edge()];
    he = he.next();
    double b = edgeLengths[he.edge()];
    he = he.next();
    double c = edgeLengths[he.edge()];

    GC_SAFETY_ASSERT(he.next() == f.halfedge(), "faces mush be triangular");

    // Herons formula
    double s = (a + b + c) / 2.0;
    double arg = s * (s - a) * (s - b) * (s - c);
    arg = std::fmax(0., arg);
    double area = std::sqrt(arg);

    faceAreas[f] = area;
  }
}
void IntrinsicGeometryInterface::requireFaceAreas() { faceAreasQ.require(); }
void IntrinsicGeometryInterface::unrequireFaceAreas() { faceAreasQ.unrequire(); }


// Vertex dual area
void IntrinsicGeometryInterface::computeVertexDualAreas() {
  faceAreasQ.ensureHave();

  vertexDualAreas = VertexData<double>(mesh, 0.);

  for (Face f : mesh.faces()) {
    double A = faceAreas[f];
    for (Vertex v : f.adjacentVertices()) {
      vertexDualAreas[v] += A / 3.0;
    }
  }
}
void IntrinsicGeometryInterface::requireVertexDualAreas() { vertexDualAreasQ.require(); }
void IntrinsicGeometryInterface::unrequireVertexDualAreas() { vertexDualAreasQ.unrequire(); }


void IntrinsicGeometryInterface::computeCornerAngles() {
  edgeLengthsQ.ensureHave();

  cornerAngles = CornerData<double>(mesh);

  for (Corner c : mesh.corners()) {
    // WARNING: Logic duplicated between cached and immediate version
    Halfedge heA = c.halfedge();
    Halfedge heOpp = heA.next();
    Halfedge heB = heOpp.next();

    GC_SAFETY_ASSERT(heB.next() == heA, "faces mush be triangular");

    double lOpp = edgeLengths[heOpp.edge()];
    double lA = edgeLengths[heA.edge()];
    double lB = edgeLengths[heB.edge()];

    double q = (lA * lA + lB * lB - lOpp * lOpp) / (2. * lA * lB);
    q = clamp(q, -1.0, 1.0);
    double angle = std::acos(q);

    cornerAngles[c] = angle;
  }
}
void IntrinsicGeometryInterface::requireCornerAngles() { cornerAnglesQ.require(); }
void IntrinsicGeometryInterface::unrequireCornerAngles() { cornerAnglesQ.unrequire(); }


// Vertex angle sums
void IntrinsicGeometryInterface::computeVertexAngleSums() {
  cornerAnglesQ.ensureHave();

  vertexAngleSums = VertexData<double>(mesh, 0.);
  for (Corner c : mesh.corners()) {
    vertexAngleSums[c.vertex()] += cornerAngles[c];
  }
}
void IntrinsicGeometryInterface::requireVertexAngleSums() { vertexAngleSumsQ.require(); }
void IntrinsicGeometryInterface::unrequireVertexAngleSums() { vertexAngleSumsQ.unrequire(); }


// Corner scaled angles
void IntrinsicGeometryInterface::computeCornerScaledAngles() {
  cornerAnglesQ.ensureHave();
  vertexAngleSumsQ.ensureHave();

  cornerScaledAngles = CornerData<double>(mesh);

  for (Corner c : mesh.corners()) {
    if (c.vertex().isBoundary()) {
      double s = PI / vertexAngleSums[c.vertex()];
      cornerScaledAngles[c] = s * cornerAngles[c];
    } else {
      double s = 2.0 * PI / vertexAngleSums[c.vertex()];
      cornerScaledAngles[c] = s * cornerAngles[c];
    }
  }
}
void IntrinsicGeometryInterface::requireCornerScaledAngles() { cornerScaledAnglesQ.require(); }
void IntrinsicGeometryInterface::unrequireCornerScaledAngles() { cornerScaledAnglesQ.unrequire(); }


// Vertex gaussian curvatures
void IntrinsicGeometryInterface::computeVertexGaussianCurvatures() {
  vertexAngleSumsQ.ensureHave();

  vertexGaussianCurvatures = VertexData<double>(mesh, 0);

  for (Vertex v : mesh.vertices()) {
    if (!v.isBoundary()) {
      vertexGaussianCurvatures[v] = 2. * PI - vertexAngleSums[v];
    }
  }
}
void IntrinsicGeometryInterface::requireVertexGaussianCurvatures() { vertexGaussianCurvaturesQ.require(); }
void IntrinsicGeometryInterface::unrequireVertexGaussianCurvatures() { vertexGaussianCurvaturesQ.unrequire(); }

// Face gaussian curvatures
void IntrinsicGeometryInterface::computeFaceGaussianCurvatures() {
  cornerScaledAnglesQ.ensureHave();

  faceGaussianCurvatures = FaceData<double>(mesh);

  for (Face f : mesh.faces()) {

    double angleDefect = -PI;
    Halfedge he = f.halfedge();
    for (int i = 0; i < 3; i++) {
      angleDefect += cornerScaledAngles[he.corner()];
      he = he.next();
    }
    GC_SAFETY_ASSERT(he == f.halfedge(), "faces mush be triangular");

    faceGaussianCurvatures[f] = angleDefect;
  }
}
void IntrinsicGeometryInterface::requireFaceGaussianCurvatures() { faceGaussianCurvaturesQ.require(); }
void IntrinsicGeometryInterface::unrequireFaceGaussianCurvatures() { faceGaussianCurvaturesQ.unrequire(); }

// Halfedge cotan weights
void IntrinsicGeometryInterface::computeHalfedgeCotanWeights() {
  edgeLengthsQ.ensureHave();
  faceAreasQ.ensureHave();

  halfedgeCotanWeights = HalfedgeData<double>(mesh, 0.);

  for (Halfedge he : mesh.interiorHalfedges()) {

    Halfedge heF = he;
    double l_ij = edgeLengths[heF.edge()];
    heF = heF.next();
    double l_jk = edgeLengths[heF.edge()];
    heF = heF.next();
    double l_ki = edgeLengths[heF.edge()];
    heF = heF.next();

    GC_SAFETY_ASSERT(heF == he, "faces mush be triangular");

    double area = faceAreas[he.face()];
    double cotValue = (-l_ij * l_ij + l_jk * l_jk + l_ki * l_ki) / (4. * area);
    halfedgeCotanWeights[he] = cotValue / 2;
  }
}
void IntrinsicGeometryInterface::requireHalfedgeCotanWeights() { halfedgeCotanWeightsQ.require(); }
void IntrinsicGeometryInterface::unrequireHalfedgeCotanWeights() { halfedgeCotanWeightsQ.unrequire(); }

// Edge cotan weights
void IntrinsicGeometryInterface::computeEdgeCotanWeights() {
  edgeLengthsQ.ensureHave();
  faceAreasQ.ensureHave();

  edgeCotanWeights = EdgeData<double>(mesh, 0.);

  for (Edge e : mesh.edges()) {
    // WARNING: Logic duplicated between cached and immediate version
    double cotSum = 0.;
    for (Halfedge he : e.adjacentInteriorHalfedges()) {
      Halfedge heFirst = he;
      double l_ij = edgeLengths[he.edge()];
      he = he.next();
      double l_jk = edgeLengths[he.edge()];
      he = he.next();
      double l_ki = edgeLengths[he.edge()];
      he = he.next();
      GC_SAFETY_ASSERT(he == heFirst, "faces mush be triangular");
      double area = faceAreas[he.face()];
      double cotValue = (-l_ij * l_ij + l_jk * l_jk + l_ki * l_ki) / (4. * area);
      cotSum += cotValue / 2;
    }
    edgeCotanWeights[e] = cotSum;
  }
}
void IntrinsicGeometryInterface::requireEdgeCotanWeights() { edgeCotanWeightsQ.require(); }
void IntrinsicGeometryInterface::unrequireEdgeCotanWeights() { edgeCotanWeightsQ.unrequire(); }

// Shape length scale
void IntrinsicGeometryInterface::computeShapeLengthScale() {
  faceAreasQ.ensureHave();

  double totalArea = 0.;
  for (Face f : mesh.faces()) {
    totalArea += faceAreas[f];
  }

  shapeLengthScale = std::sqrt(totalArea);
}
void IntrinsicGeometryInterface::requireShapeLengthScale() { shapeLengthScaleQ.require(); }
void IntrinsicGeometryInterface::unrequireShapeLengthScale() { shapeLengthScaleQ.unrequire(); }

// Mesh length scale
void IntrinsicGeometryInterface::computeMeshLengthScale() {
  edgeLengthsQ.ensureHave();

  double totalEdgeLength = 0.;
  for (Edge e : mesh.edges()) {
    totalEdgeLength += edgeLengths[e];
  }

  meshLengthScale = totalEdgeLength / mesh.nEdges();
}
void IntrinsicGeometryInterface::requireMeshLengthScale() { meshLengthScaleQ.require(); }
void IntrinsicGeometryInterface::unrequireMeshLengthScale() { meshLengthScaleQ.unrequire(); }


// Halfedge vectors in face
void IntrinsicGeometryInterface::computeHalfedgeVectorsInFace() {
  edgeLengthsQ.ensureHave();
  faceAreasQ.ensureHave();

  halfedgeVectorsInFace = HalfedgeData<Vector2>(mesh);

  for (Face f : mesh.faces()) {

    // Gather some values
    Halfedge heAB = f.halfedge();
    Halfedge heBC = heAB.next();
    Halfedge heCA = heBC.next();
    GC_SAFETY_ASSERT(heCA.next() == heAB, "faces must be triangular");

    double lAB = edgeLengths[heAB.edge()];
    double lBC = edgeLengths[heBC.edge()];
    double lCA = edgeLengths[heCA.edge()];

    // Assign positions to all three vertices
    // Vector2 pA{0., 0.}; // used implicitly
    Vector2 pB{lAB, 0.};
    // pC is the hard one:

    double tArea = faceAreas[f];

    // Compute width and height of right triangle formed via altitude from C
    double h = 2. * tArea / lAB;
    double w = std::sqrt(std::max(0., lCA * lCA - h * h));

    // Take the closer of the positive and negative solutions
    if (lBC * lBC > (lAB * lAB + lCA * lCA)) w *= -1.0;

    // Project some vectors to get the actual position
    Vector2 pC{w, h};

    // Now, all halfedge vectors are just coordinates
    halfedgeVectorsInFace[heAB] = pB;
    halfedgeVectorsInFace[heBC] = pC - pB;
    halfedgeVectorsInFace[heCA] = -pC;
  }

  // Set all the exterior ones to NaN
  for (Halfedge he : mesh.exteriorHalfedges()) {
    halfedgeVectorsInFace[he] = Vector2::undefined();
  }
}
void IntrinsicGeometryInterface::requireHalfedgeVectorsInFace() { halfedgeVectorsInFaceQ.require(); }
void IntrinsicGeometryInterface::unrequireHalfedgeVectorsInFace() { halfedgeVectorsInFaceQ.unrequire(); }

// Transport vectors across halfedge
void IntrinsicGeometryInterface::computeTransportVectorsAcrossHalfedge() {
  halfedgeVectorsInFaceQ.ensureHave();

  transportVectorsAcrossHalfedge = HalfedgeData<Vector2>(mesh, Vector2::undefined());

  for (Edge e : mesh.edges()) {
    if (e.isBoundary()) continue;

    Halfedge heA = e.halfedge();
    Halfedge heB = heA.twin();

    Vector2 vecA = halfedgeVectorsInFace[heA];
    Vector2 vecB = halfedgeVectorsInFace[heB];
    Vector2 rot = unit(-vecB / vecA);

    transportVectorsAcrossHalfedge[heA] = rot;
    transportVectorsAcrossHalfedge[heB] = rot.inv();
  }
}
void IntrinsicGeometryInterface::requireTransportVectorsAcrossHalfedge() { transportVectorsAcrossHalfedgeQ.require(); }
void IntrinsicGeometryInterface::unrequireTransportVectorsAcrossHalfedge() {
  transportVectorsAcrossHalfedgeQ.unrequire();
}


// Halfedge vectors in vertex
void IntrinsicGeometryInterface::computeHalfedgeVectorsInVertex() {

  if (!mesh.usesImplicitTwin()) {
    throw std::runtime_error("ERROR: Tangent spaces not implemented for general SurfaceMesh, use ManifoldSurfaceMesh");
  }

  edgeLengthsQ.ensureHave();
  cornerScaledAnglesQ.ensureHave();

  halfedgeVectorsInVertex = HalfedgeData<Vector2>(mesh);

  for (Vertex v : mesh.vertices()) {
    double coordSum = 0.0;

    // Custom loop to orbit CCW
    Halfedge firstHe = v.halfedge();
    Halfedge currHe = firstHe;
    do {
      halfedgeVectorsInVertex[currHe] = Vector2::fromAngle(coordSum) * edgeLengths[currHe.edge()];
      coordSum += cornerScaledAngles[currHe.corner()];
      if (!currHe.isInterior()) break;
      currHe = currHe.next().next().twin();
    } while (currHe != firstHe);
  }
}
void IntrinsicGeometryInterface::requireHalfedgeVectorsInVertex() { halfedgeVectorsInVertexQ.require(); }
void IntrinsicGeometryInterface::unrequireHalfedgeVectorsInVertex() { halfedgeVectorsInVertexQ.unrequire(); }


// Transport vectors along halfedge
void IntrinsicGeometryInterface::computeTransportVectorsAlongHalfedge() {
  halfedgeVectorsInVertexQ.ensureHave();

  transportVectorsAlongHalfedge = HalfedgeData<Vector2>(mesh);

  for (Edge e : mesh.edges()) {

    Halfedge heA = e.halfedge();
    Halfedge heB = heA.twin();

    Vector2 vecA = halfedgeVectorsInVertex[heA];
    Vector2 vecB = halfedgeVectorsInVertex[heB];
    Vector2 rot = unit(-vecB / vecA);

    transportVectorsAlongHalfedge[heA] = rot;
    transportVectorsAlongHalfedge[heB] = rot.inv();
  }
}
void IntrinsicGeometryInterface::requireTransportVectorsAlongHalfedge() { transportVectorsAlongHalfedgeQ.require(); }
void IntrinsicGeometryInterface::unrequireTransportVectorsAlongHalfedge() {
  transportVectorsAlongHalfedgeQ.unrequire();
}


// Cotan Laplacian
void IntrinsicGeometryInterface::computeCotanLaplacian() {
  vertexIndicesQ.ensureHave();
  edgeCotanWeightsQ.ensureHave();

  cotanLaplacian = Eigen::SparseMatrix<double>(mesh.nVertices(), mesh.nVertices());
  std::vector<Eigen::Triplet<double>> tripletList;

  for (Edge e : mesh.edges()) {
    Halfedge he = e.halfedge();
    Vertex vTail = he.vertex();
    Vertex vHead = he.next().vertex();

    size_t iVHead = vertexIndices[vHead];
    size_t iVTail = vertexIndices[vTail];

    double weight = edgeCotanWeights[e];

    tripletList.emplace_back(iVTail, iVTail, weight);
    tripletList.emplace_back(iVHead, iVHead, weight);
    tripletList.emplace_back(iVTail, iVHead, -weight);
    tripletList.emplace_back(iVHead, iVTail, -weight);
  }

  cotanLaplacian.setFromTriplets(tripletList.begin(), tripletList.end());
}
void IntrinsicGeometryInterface::requireCotanLaplacian() { cotanLaplacianQ.require(); }
void IntrinsicGeometryInterface::unrequireCotanLaplacian() { cotanLaplacianQ.unrequire(); }


// Vertex lumped mass matrix
void IntrinsicGeometryInterface::computeVertexLumpedMassMatrix() {
  vertexDualAreasQ.ensureHave();

  size_t nVerts = mesh.nVertices();
  Eigen::VectorXd hodge0V(nVerts);
  size_t iV = 0;
  for (Vertex v : mesh.vertices()) {
    double dualArea = vertexDualAreas[v];
    hodge0V[iV] = dualArea;
    iV++;
  }
  vertexLumpedMassMatrix = hodge0V.asDiagonal();
}
void IntrinsicGeometryInterface::requireVertexLumpedMassMatrix() { vertexLumpedMassMatrixQ.require(); }
void IntrinsicGeometryInterface::unrequireVertexLumpedMassMatrix() { vertexLumpedMassMatrixQ.unrequire(); }


// Vertex Galerkin mass matrix
void IntrinsicGeometryInterface::computeVertexGalerkinMassMatrix() {
  vertexIndicesQ.ensureHave();
  faceAreasQ.ensureHave();

  std::vector<Eigen::Triplet<double>> triplets;
  for (Face f : mesh.faces()) {
    double area = faceAreas[f];

    // Gather indices for vertices on faces
    Halfedge he = f.halfedge();
    Vertex vA = he.vertex();
    he = he.next();
    Vertex vB = he.vertex();
    he = he.next();
    Vertex vC = he.vertex();
    GC_SAFETY_ASSERT(he.next() == f.halfedge(), "faces must be triangular");

    std::array<size_t, 3> indices{vertexIndices[vA], vertexIndices[vB], vertexIndices[vC]};

    // Set entries
    for (int root = 0; root < 3; root++) {
      size_t i = indices[root];
      size_t j = indices[(root + 1) % 3];
      size_t k = indices[(root + 2) % 3];
      triplets.emplace_back(i, i, area / 6.);
      triplets.emplace_back(i, j, area / 12.);
      triplets.emplace_back(i, k, area / 12.);
    }
  }

  vertexGalerkinMassMatrix = Eigen::SparseMatrix<double>(mesh.nVertices(), mesh.nVertices());
  vertexGalerkinMassMatrix.setFromTriplets(triplets.begin(), triplets.end());
}
void IntrinsicGeometryInterface::requireVertexGalerkinMassMatrix() { vertexGalerkinMassMatrixQ.require(); }
void IntrinsicGeometryInterface::unrequireVertexGalerkinMassMatrix() { vertexGalerkinMassMatrixQ.unrequire(); }

// Face Galerkin mass matrix
void IntrinsicGeometryInterface::computeFaceGalerkinMassMatrix() {
  faceIndicesQ.ensureHave();
  faceAreasQ.ensureHave();

  std::vector<Eigen::Triplet<double>> triplets;
  for (Face f : mesh.faces()) {
    size_t i = faceIndices[f];
    triplets.emplace_back(i, i, faceAreas[f]);
  }

  faceGalerkinMassMatrix = Eigen::SparseMatrix<double>(mesh.nFaces(), mesh.nFaces());
  faceGalerkinMassMatrix.setFromTriplets(triplets.begin(), triplets.end());
}
void IntrinsicGeometryInterface::requireFaceGalerkinMassMatrix() { faceGalerkinMassMatrixQ.require(); }
void IntrinsicGeometryInterface::unrequireFaceGalerkinMassMatrix() { faceGalerkinMassMatrixQ.unrequire(); }


// Vertex connection Laplacian
void IntrinsicGeometryInterface::computeVertexConnectionLaplacian() {
  vertexIndicesQ.ensureHave();
  edgeCotanWeightsQ.ensureHave();
  transportVectorsAlongHalfedgeQ.ensureHave();


  std::vector<Eigen::Triplet<std::complex<double>>> triplets;
  for (Halfedge he : mesh.halfedges()) {

    size_t iTail = vertexIndices[he.vertex()];
    size_t iTip = vertexIndices[he.next().vertex()];

    double weight = edgeCotanWeights[he.edge()];
    Vector2 rot = transportVectorsAlongHalfedge[he.twin()];
    triplets.emplace_back(iTail, iTail, weight);
    triplets.emplace_back(iTail, iTip, -weight * rot);
  }

  vertexConnectionLaplacian = Eigen::SparseMatrix<std::complex<double>>(mesh.nVertices(), mesh.nVertices());
  vertexConnectionLaplacian.setFromTriplets(triplets.begin(), triplets.end());
}
void IntrinsicGeometryInterface::requireVertexConnectionLaplacian() { vertexConnectionLaplacianQ.require(); }
void IntrinsicGeometryInterface::unrequireVertexConnectionLaplacian() { vertexConnectionLaplacianQ.unrequire(); }

// Face connection Laplacian
void IntrinsicGeometryInterface::computeFaceConnectionLaplacian() {
  faceIndicesQ.ensureHave();
  transportVectorsAcrossHalfedgeQ.ensureHave();

  std::vector<Eigen::Triplet<std::complex<double>>> triplets;
  for (Face f : mesh.faces()) {
    size_t i = faceIndices[f];

    std::complex<double> weightISum = 0;
    for (Halfedge he : f.adjacentHalfedges()) {

      if (!he.twin().isInterior()) {
        continue;
      }

      Face neighFace = he.twin().face();
      unsigned int j = faceIndices[neighFace];

      // LC connection between the faces
      Vector2 rot = transportVectorsAcrossHalfedge[he.twin()];
      double weight = 1; // FIXME TODO figure out weights
      triplets.emplace_back(i, j, -weight * rot);

      weightISum += weight;
    }

    triplets.emplace_back(i, i, weightISum);
  }
  faceConnectionLaplacian = Eigen::SparseMatrix<std::complex<double>>(mesh.nFaces(), mesh.nFaces());
  faceConnectionLaplacian.setFromTriplets(triplets.begin(), triplets.end());
}
void IntrinsicGeometryInterface::requireFaceConnectionLaplacian() { faceConnectionLaplacianQ.require(); }
void IntrinsicGeometryInterface::unrequireFaceConnectionLaplacian() { faceConnectionLaplacianQ.unrequire(); }


void IntrinsicGeometryInterface::computeDECOperators() {
  vertexIndicesQ.ensureHave();
  edgeIndicesQ.ensureHave();
  faceIndicesQ.ensureHave();
  vertexDualAreasQ.ensureHave();
  edgeCotanWeightsQ.ensureHave();
  faceAreasQ.ensureHave();

  size_t nVerts = mesh.nVertices();
  size_t nEdges = mesh.nEdges();
  size_t nFaces = mesh.nFaces();

  { // Hodge 0
    Eigen::VectorXd hodge0V(nVerts);
    for (Vertex v : mesh.vertices()) {
      double primalArea = 1.0;
      double dualArea = vertexDualAreas[v];
      double ratio = dualArea / primalArea;
      size_t iV = vertexIndices[v];
      hodge0V[iV] = ratio;
    }

    hodge0 = hodge0V.asDiagonal();
    hodge0Inverse = hodge0V.asDiagonal().inverse();
  }


  { // Hodge 1
    Eigen::VectorXd hodge1V(nEdges);
    for (Edge e : mesh.edges()) {
      double ratio = edgeCotanWeights[e];
      size_t iE = edgeIndices[e];
      hodge1V[iE] = ratio;
    }

    hodge1 = hodge1V.asDiagonal();
    hodge1Inverse = hodge1V.asDiagonal().inverse();
  }

  { // Hodge 2
    Eigen::VectorXd hodge2V(nFaces);
    for (Face f : mesh.faces()) {
      double primalArea = faceAreas[f];
      double dualArea = 1.0;
      double ratio = dualArea / primalArea;

      size_t iF = faceIndices[f];
      hodge2V[iF] = ratio;
    }
    hodge2 = hodge2V.asDiagonal();
    hodge2Inverse = hodge2V.asDiagonal().inverse();
  }


  { // D0
    d0 = Eigen::SparseMatrix<double>(nEdges, nVerts);
    std::vector<Eigen::Triplet<double>> tripletList;

    for (Edge e : mesh.edges()) {
      size_t iEdge = edgeIndices[e];
      Halfedge he = e.halfedge();
      Vertex vTail = he.vertex();
      Vertex vHead = he.next().vertex();

      size_t iVHead = vertexIndices[vHead];
      tripletList.emplace_back(iEdge, iVHead, 1.0);

      size_t iVTail = vertexIndices[vTail];
      tripletList.emplace_back(iEdge, iVTail, -1.0);
    }

    d0.setFromTriplets(tripletList.begin(), tripletList.end());
  }

  { // D1
    d1 = Eigen::SparseMatrix<double>(nFaces, nEdges);
    std::vector<Eigen::Triplet<double>> tripletList;

    for (Face f : mesh.faces()) {
      size_t iFace = faceIndices[f];

      for (Halfedge he : f.adjacentHalfedges()) {
        size_t iEdge = edgeIndices[he.edge()];
        double sign = (he == he.edge().halfedge()) ? (1.0) : (-1.0);
        tripletList.emplace_back(iFace, iEdge, sign);
      }
    }

    d1.setFromTriplets(tripletList.begin(), tripletList.end());
  }
}
void IntrinsicGeometryInterface::requireDECOperators() { DECOperatorsQ.require(); }
void IntrinsicGeometryInterface::unrequireDECOperators() { DECOperatorsQ.unrequire(); }

} // namespace surface
} // namespace geometrycentral
