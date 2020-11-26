#include <geometrycentral/surface/direction_fields.h>

#include "geometrycentral/numerical/linear_solvers.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>

using std::cout;
using std::endl;

namespace geometrycentral {
namespace surface {

namespace {

SparseMatrix<std::complex<double>> computeVertexConnectionLaplacian(IntrinsicGeometryInterface& geometry, int nSym) {

  SurfaceMesh& mesh = geometry.mesh;

  geometry.requireVertexIndices();
  geometry.requireEdgeCotanWeights();
  geometry.requireTransportVectorsAlongHalfedge();

  std::vector<Eigen::Triplet<std::complex<double>>> triplets;
  for (Halfedge he : mesh.halfedges()) {

    size_t iTail = geometry.vertexIndices[he.vertex()];
    size_t iTip = geometry.vertexIndices[he.next().vertex()];

    // Levi-Civita connection between vertices
    Vector2 rot = geometry.transportVectorsAlongHalfedge[he.twin()].pow(nSym);
    double weight = geometry.edgeCotanWeights[he.edge()];

    triplets.emplace_back(iTail, iTail, weight);
    triplets.emplace_back(iTail, iTip, -weight * rot);
  }
  // assemble matrix from triplets
  Eigen::SparseMatrix<std::complex<double>> vertexConnectionLaplacian(mesh.nVertices(), mesh.nVertices());
  vertexConnectionLaplacian.setFromTriplets(triplets.begin(), triplets.end());

  // Shift to avoid singularity
  Eigen::SparseMatrix<std::complex<double>> eye(mesh.nVertices(), mesh.nVertices());
  eye.setIdentity();
  vertexConnectionLaplacian += 1e-9 * eye;

  return vertexConnectionLaplacian;
}

SparseMatrix<std::complex<double>> computeFaceConnectionLaplacian(IntrinsicGeometryInterface& geometry, int nSym) {

  SurfaceMesh& mesh = geometry.mesh;

  geometry.requireFaceIndices();
  geometry.requireTransportVectorsAcrossHalfedge();

  std::vector<Eigen::Triplet<std::complex<double>>> triplets;
  for (Face f : mesh.faces()) {
    size_t i = geometry.faceIndices[f];

    std::complex<double> weightDiagSum = 0;
    for (Halfedge he : f.adjacentHalfedges()) {
      if (he.twin().isInterior()) {

        Face neighFace = he.twin().face();
        size_t j = geometry.faceIndices[neighFace];

        // Levi-Civita connection between the faces
        Vector2 rot = geometry.transportVectorsAcrossHalfedge[he.twin()].pow(nSym);

        double weight = 1;
        triplets.emplace_back(i, j, -weight * rot);

        weightDiagSum += weight;
      }
    }

    triplets.emplace_back(i, i, weightDiagSum + 1e-9); // Shift to avoid singularity
  }
  // assemble matrix from triplets
  Eigen::SparseMatrix<std::complex<double>> faceConnectionLaplacian(mesh.nFaces(), mesh.nFaces());
  faceConnectionLaplacian.setFromTriplets(triplets.begin(), triplets.end());

  return faceConnectionLaplacian;
}

} // namespace

VertexData<Vector2> computeSmoothestVertexDirectionField(IntrinsicGeometryInterface& geometry, int nSym) {


  SurfaceMesh& mesh = geometry.mesh;

  geometry.requireVertexIndices();
  geometry.requireVertexGalerkinMassMatrix();

  // Mass matrix
  SparseMatrix<std::complex<double>> massMatrix = geometry.vertexGalerkinMassMatrix.cast<std::complex<double>>();

  // Energy matrix
  SparseMatrix<std::complex<double>> energyMatrix = computeVertexConnectionLaplacian(geometry, nSym);

  // Find the smallest eigenvector
  Vector<std::complex<double>> solution = smallestEigenvectorSquare(energyMatrix, massMatrix);

  // Copy the result to a VertexData vector
  VertexData<Vector2> toReturn(mesh);
  for (Vertex v : mesh.vertices()) {
    toReturn[v] = Vector2::fromComplex(solution(geometry.vertexIndices[v]));
    toReturn[v] = unit(toReturn[v]);
  }

  return toReturn;
}

FaceData<Vector2> computeSmoothestFaceDirectionField(IntrinsicGeometryInterface& geometry, int nSym) {

  SurfaceMesh& mesh = geometry.mesh;

  geometry.requireFaceIndices();
  geometry.requireFaceGalerkinMassMatrix();

  // Mass matrix
  SparseMatrix<std::complex<double>> massMatrix = geometry.faceGalerkinMassMatrix.cast<std::complex<double>>();

  // Energy matrix
  SparseMatrix<std::complex<double>> energyMatrix = computeFaceConnectionLaplacian(geometry, nSym);

  // Find the smallest eigenvector
  Vector<std::complex<double>> solution = smallestEigenvectorSquare(energyMatrix, massMatrix);

  // Copy the result to a FaceData vector
  FaceData<Vector2> toReturn(mesh);
  for (Face f : mesh.faces()) {
    toReturn[f] = Vector2::fromComplex(solution(geometry.faceIndices[f]));
    toReturn[f] = unit(toReturn[f]);
  }

  return toReturn;
}

VertexData<Vector2> computeSmoothestBoundaryAlignedVertexDirectionField(IntrinsicGeometryInterface& geometry,
                                                                        int nSym) {
  SurfaceMesh& mesh = geometry.mesh;

  if (!mesh.hasBoundary()) {
    throw std::logic_error("tried to compute smoothest boundary aligned direction field on a mesh without boundary");
  }

  geometry.requireVertexGalerkinMassMatrix();
  geometry.requireHalfedgeVectorsInVertex();

  // Mass matrix
  SparseMatrix<std::complex<double>> massMatrix = geometry.vertexGalerkinMassMatrix.cast<std::complex<double>>();

  // Energy matrix
  SparseMatrix<std::complex<double>> energyMatrix = computeVertexConnectionLaplacian(geometry, nSym);

  // Compute the boundary values
  VertexData<std::complex<double>> boundaryValues(mesh);
  for (Vertex v : mesh.vertices()) {
    if (v.isBoundary()) {

      // Find incoming and outgoing boundary vectors as tangent
      Halfedge heBoundaryA = v.halfedge();
      Halfedge heBoundaryB = heBoundaryA.twin().next();

      Vector2 vecA = geometry.halfedgeVectorsInVertex[heBoundaryA];
      Vector2 vecB = geometry.halfedgeVectorsInVertex[heBoundaryB];

      Vector2 tangentV = unit(-vecA + vecB);
      Vector2 normalV = tangentV.rotate90();

      boundaryValues[v] = normalV.pow(nSym);
    } else {
      boundaryValues[v] = 0;
    }
  }

  // Assemble right-hand side from these boundary values
  Vector<std::complex<double>> b(mesh.nVertices());
  b.setZero();

  for (Vertex v : mesh.vertices()) {
    if (!v.isBoundary()) {

      for (Halfedge he : v.incomingHalfedges()) {
        Face neighFace = he.twin().face();
        if (he.vertex().isBoundary()) {

          size_t i = geometry.vertexIndices[v];
          size_t j = geometry.vertexIndices[he.vertex()];

          // move boundary terms to the right-hand side and remove the corresponding matrix entries
          std::complex<double> Aij = energyMatrix.coeff(i, j);
          energyMatrix.coeffRef(i, j) = 0;
          std::complex<double> bVal = boundaryValues[he.vertex()];
          b(i) += -Aij * bVal;
        }
      }
    }
  }

  Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor> LHS = energyMatrix;
  Eigen::VectorXcd RHS = massMatrix * b;
  Eigen::VectorXcd solution = solveSquare(LHS, RHS);

  // Copy the result to a VertexData vector for both the boundary and interior
  VertexData<Vector2> toReturn(mesh);
  for (Vertex v : mesh.vertices()) {
    if (v.isBoundary()) {
      toReturn[v] = Vector2::fromComplex(boundaryValues[v]);
    } else {
      toReturn[v] = Vector2::fromComplex(solution(geometry.vertexIndices[v]));
      toReturn[v] = unit(toReturn[v]);
    }
  }

  return toReturn;
}

FaceData<Vector2> computeSmoothestBoundaryAlignedFaceDirectionField(IntrinsicGeometryInterface& geometry, int nSym) {

  SurfaceMesh& mesh = geometry.mesh;

  if (!mesh.hasBoundary()) {
    throw std::logic_error("tried to compute smoothest boundary aligned direction field on a mesh without boundary");
  }

  geometry.requireFaceGalerkinMassMatrix();
  geometry.requireHalfedgeVectorsInFace();

  // Mass matrix
  SparseMatrix<std::complex<double>> massMatrix = geometry.faceGalerkinMassMatrix.cast<std::complex<double>>();

  // Energy matrix
  SparseMatrix<std::complex<double>> energyMatrix = computeFaceConnectionLaplacian(geometry, nSym);

  // Compute boundary values
  FaceData<bool> isInterior(mesh);
  FaceData<std::complex<double>> boundaryValues(mesh);
  for (Face f : mesh.faces()) {
    bool isBoundary = false;
    for (Edge e : f.adjacentEdges()) {
      isBoundary |= e.isBoundary();
    }
    isInterior[f] = !isBoundary;

    if (isInterior[f]) {
      boundaryValues[f] = 0;
    } else {
      Vector2 bC = Vector2::zero();
      for (Halfedge he : f.adjacentHalfedges()) {
        if (he.edge().isBoundary()) {
          bC -= geometry.halfedgeVectorsInFace[he].rotate90(); // negate the vector to point outwards
        }
      }
      bC = unit(bC);
      boundaryValues[f] = bC.pow(nSym);
    }
  }

  // Assemble right-hand side from these boundary values
  Vector<std::complex<double>> b(mesh.nFaces());
  b.setZero();

  for (Face f : mesh.faces()) {
    if (isInterior[f]) {

      for (Halfedge he : f.adjacentHalfedges()) {
        Face neighFace = he.twin().face();
        if (!isInterior[neighFace]) {

          size_t i = geometry.faceIndices[f];
          size_t j = geometry.faceIndices[neighFace];

          // move boundary terms to the right-hand side and remove the corresponding matrix entries
          std::complex<double> Aij = energyMatrix.coeff(i, j);
          energyMatrix.coeffRef(i, j) = 0;
          std::complex<double> bVal = boundaryValues[neighFace];
          b(i) += -Aij * bVal;
        }
      }
    }
  }

  Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor> LHS = energyMatrix;
  Eigen::VectorXcd RHS = massMatrix * b;
  Eigen::VectorXcd solution = solveSquare(LHS, RHS);

  // Copy the result to a FaceData object
  FaceData<Vector2> field(mesh);
  for (Face f : mesh.faces()) {
    if (isInterior[f]) {
      field[f] = Vector2::fromComplex(solution(geometry.faceIndices[f]));
      field[f] = unit(field[f]);
    } else {
      field[f] = Vector2::fromComplex(boundaryValues[f]);
    }
  }

  return field;
}

VertexData<Vector2> computeCurvatureAlignedVertexDirectionField(ExtrinsicGeometryInterface& geometry, int nSym) {


  SurfaceMesh& mesh = geometry.mesh;
  size_t N = mesh.nVertices();

  geometry.requireVertexIndices();
  geometry.requireVertexGalerkinMassMatrix();
  geometry.requireVertexPrincipalCurvatureDirections();

  // Mass matrix
  SparseMatrix<std::complex<double>> massMatrix = geometry.vertexGalerkinMassMatrix.cast<std::complex<double>>();

  // Energy matrix
  SparseMatrix<std::complex<double>> energyMatrix = computeVertexConnectionLaplacian(geometry, nSym);

  Vector<std::complex<double>> dirVec(N);
  if (nSym == 2) {
    for (Vertex v : mesh.vertices()) {
      dirVec[geometry.vertexIndices[v]] = geometry.vertexPrincipalCurvatureDirections[v];
    }
  } else if (nSym == 4) {
    for (Vertex v : mesh.vertices()) {
      dirVec[geometry.vertexIndices[v]] =
          std::pow(std::complex<double>(geometry.vertexPrincipalCurvatureDirections[v]), 2);
    }
  } else {
    throw std::logic_error("ERROR: It only makes sense to align with curvature when nSym = 2 or 4");
  }

  // Normalize the alignment field
  double scale = std::sqrt(std::abs((dirVec.adjoint() * massMatrix * dirVec)[0]));
  dirVec /= scale;

  // this is something of a magical constant, see "Globally Optimal Direction Fields", eqn 16
  double lambdaT = 0;

  Eigen::VectorXcd RHS = massMatrix * dirVec;
  Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor> LHS = energyMatrix - lambdaT * massMatrix;
  Eigen::VectorXcd solution = solveSquare(LHS, RHS);

  // Copy the result to a VertexData vector
  VertexData<Vector2> toReturn(mesh);
  for (Vertex v : mesh.vertices()) {
    toReturn[v] = Vector2::fromComplex(solution(geometry.vertexIndices[v]));
    toReturn[v] = unit(toReturn[v]);
  }

  return toReturn;
}

FaceData<Vector2> computeCurvatureAlignedFaceDirectionField(EmbeddedGeometryInterface& geometry, int nSym) {

  SurfaceMesh& mesh = geometry.mesh;
  const unsigned int N = mesh.nFaces();

  geometry.requireFaceIndices();
  geometry.requireFaceGalerkinMassMatrix();
  geometry.requireFacePrincipalCurvatureDirections();

  // Mass matrix
  SparseMatrix<std::complex<double>> massMatrix = geometry.faceGalerkinMassMatrix.cast<std::complex<double>>();

  // Energy matrix
  SparseMatrix<std::complex<double>> energyMatrix = computeFaceConnectionLaplacian(geometry, nSym);

  Eigen::VectorXcd dirVec(N);
  if (nSym == 2) {
    for (Face f : mesh.faces()) {
      dirVec[geometry.faceIndices[f]] = geometry.facePrincipalCurvatureDirections[f];
    }
  } else if (nSym == 4) {
    for (Face f : mesh.faces()) {
      dirVec[geometry.faceIndices[f]] = std::pow(std::complex<double>(geometry.facePrincipalCurvatureDirections[f]), 2);
    }
  } else {
    throw std::logic_error("ERROR: It only makes sense to align with curvature when nSym = 2 or 4");
  }

  // Normalize the alignment field
  double scale = std::sqrt(std::abs((dirVec.adjoint() * massMatrix * dirVec)[0]));
  dirVec /= scale;

  double lambdaT = 0.0; // this is something of a magical constant, see "Globally Optimal Direction Fields", eqn 16

  Eigen::VectorXcd RHS = massMatrix * dirVec;
  Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor> LHS = energyMatrix - lambdaT * massMatrix;
  Eigen::VectorXcd solution = solveSquare(LHS, RHS);

  // Copy the result to a FaceData object
  FaceData<Vector2> field(mesh);
  for (Face f : mesh.faces()) {
    field[f] = Vector2::fromComplex(solution[geometry.faceIndices[f]]);
    field[f] = unit(field[f]);
  }

  return field;
}


FaceData<int> computeFaceIndex(IntrinsicGeometryInterface& geometry, const VertexData<Vector2>& directionField,
                               int nSym) {

  SurfaceMesh& mesh = geometry.mesh;

  geometry.requireTransportVectorsAlongHalfedge();
  geometry.requireFaceGaussianCurvatures();

  // Store the result here
  FaceData<int> indices(mesh);

  // TODO haven't tested that this correctly reports the index when it is larger
  // than +-1

  for (Face f : mesh.faces()) {
    // Trace the direction field around the face and see how many times it
    // spins!
    double totalRot = geometry.faceGaussianCurvatures[f] * nSym;

    for (Halfedge he : f.adjacentHalfedges()) {
      if (he.twin().isInterior()) {
        // Compute the rotation along the halfedge implied by the field
        Vector2 x0 = directionField[he.vertex()];
        Vector2 x1 = directionField[he.twin().vertex()];
        Vector2 rot = geometry.transportVectorsAlongHalfedge[he].pow(nSym);

        // Find the difference in angle
        double theta0 = (rot * x0).arg();
        double theta1 = x1.arg();
        double deltaTheta = regularizeAngle(theta1 - theta0 + PI) - PI; // regularize to [-PI,PI]

        totalRot += deltaTheta; // accumulate
      }
    }

    // Compute the net rotation and corresponding index
    int index = static_cast<int>(std::round(totalRot / (2 * PI))); // should be very close to a multiple of 2PI
    indices[f] = index;
  }

  return indices;
}


VertexData<int> computeVertexIndex(IntrinsicGeometryInterface& geometry, const FaceData<Vector2>& directionField,
                                   int nSym) {

  SurfaceMesh& mesh = geometry.mesh;

  geometry.requireTransportVectorsAcrossHalfedge();
  geometry.requireVertexGaussianCurvatures();

  // Store the result here
  VertexData<int> indices(mesh);

  // TODO haven't tested that this correctly reports the index when it is larger
  // than +-1

  for (Vertex v : mesh.vertices()) {

    // Trace the direction field around the face and see how many times it
    // spins!
    double totalRot = geometry.vertexGaussianCurvatures[v] * nSym;

    if (!v.isBoundary()) {
      for (Halfedge he : v.incomingHalfedges()) {
        // Compute the rotation along the halfedge implied by the field
        Vector2 x0 = directionField[he.face()];
        Vector2 x1 = directionField[he.twin().face()];
        Vector2 rot = geometry.transportVectorsAcrossHalfedge[he].pow(nSym);

        double deltaTheta = regularizeAngle(x1.arg() - (rot * x0).arg() + PI) - PI; // regularize to [-PI,PI]

        totalRot += deltaTheta;
      }
    }

    // Compute the net rotation and corresponding index
    int index = static_cast<int>(std::round(totalRot / (2 * PI))); // should be very close to a multiple of 2PI
    indices[v] = index;
  }

  return indices;
}
/*
VertexData<Vector2> computeSmoothestBoundaryAlignedVertexDirectionField(IntrinsicGeometryInterface& geometry,
                                                                        int nSym) {
  SurfaceMesh& mesh = geometry.mesh;


  if (!mesh.hasBoundary()) {
    throw std::logic_error("tried to compute smoothest boundary aligned direction field on a mesh without boundary");
  }

  geometry.requireVertexGalerkinMassMatrix();
  geometry.requireVertexConnectionLaplacian();
  geometry.requireHalfedgeVectorsInVertex();

  // Mass matrix
  SparseMatrix<std::complex<double>> massMatrix = geometry.vertexGalerkinMassMatrix.cast<std::complex<double>>();

  // Energy matrix
  SparseMatrix<std::complex<double>> energyMatrix = geometry.vertexConnectionLaplacian;

  // Compute the boundary values
  VertexData<std::complex<double>> boundaryValues(mesh);
  VertexData<char> isBoundary(mesh, false);
  for (Vertex v : mesh.vertices()) {
    if (v.isBoundary()) {
      isBoundary[v] = true;

      // Find incoming and outgoing boundary vectors as tangent
      Halfedge heBoundaryA = v.halfedge();
      Halfedge heBoundaryB = heBoundaryA.twin().next();

      Vector2 vecA = geometry.halfedgeVectorsInVertex[heBoundaryA];
      Vector2 vecB = geometry.halfedgeVectorsInVertex[heBoundaryB];

      Vector2 tangentV = unit(-vecA + vecB);
      Vector2 normalV = tangentV.rotate90();

      boundaryValues[v] = normalV.pow(nSym);
    } else {
      boundaryValues[v] = 0;
    }
  }
  Vector<std::complex<double>> b = boundaryValues.toVector();
  Vector<char> isBoundaryVec = isBoundary.toVector();

  // Block decompose problem


  // Compute the actual solution
  std::cout << "Solving linear problem..." << std::endl;

  // Store the solution here
  Eigen::VectorXcd solution;

  // If requested, align to principal curvatures
  if (alignCurvature) {

    geometry.requirePrincipalDirections();

    Eigen::VectorXcd dirVec(nInterior);
    for (Vertex v : mesh.vertices()) {
      if (v.isBoundary()) {
        continue;
      }

      Vector2 directionVal = geometry.principalDirections[v];
      if (nSym == 4) {
        directionVal = std::pow(directionVal, 2);
      }

      // Normalize the curvature vectors. By doing so, we lose the property of adjusting the strength of the alignment
      // based on the strength of the curvature, but resolve any scaling issues between the magnitude of the normals and
      // the magnitude of the desired field.  Be careful when interpreting this as opposed to the usual direction field
      // optimization.
      dirVec[geometry.interiorVertexIndices[v]] = directionVal / std::abs(directionVal);
    }

    double t = 0.01; // this is something of a magical constant, see "Globally
                     // Optimal Direction Fields", eqn 9

    Eigen::VectorXcd RHS = massMatrix * (t * dirVec + b);
    Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor> LHS = energyMatrix;
    solution = solveSquare(LHS, RHS);
  }
  // Otherwise find the general closest solution
  else {
    std::cout << "Solving smoothest field dirichlet problem..." << std::endl;
    Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor> LHS = energyMatrix;
    Eigen::VectorXcd RHS = massMatrix * b;
    solution = solveSquare(LHS, RHS);
  }

  // Copy the result to a VertexData vector for both the boudary and interior
  VertexData<Vector2> toReturn(*mesh);
  for (Vertex v : mesh.vertices()) {
    if (v.isBoundary()) {
      toReturn[v] = boundaryValues[v];
    } else {
      toReturn[v] = unit(solution[geometry.interiorVertexIndices[v]]);
    }
  }

  return toReturn;
}

// Helpers for computing face-based direction fields
namespace {

FaceData<Vector2> computeSmoothestFaceDirectionField_noBoundary(IntrinsicGeometryInterface& geometry, int nSym = 1,
                                                                bool alignCurvature = false) {


  SurfaceMesh* mesh = geometry->getMesh();
  unsigned int N = mesh.nFaces();

  GeometryCache<Euclidean>& geometry = geometry->cache;
  geometry.requireFaceTransportCoefs();
  geometry.requireFaceNormals();
  geometry.requireFaceAreas();
  geometry.requireDihedralAngles();
  geometry.requireFaceIndices();

  // === Allocate matrices
  // Energy matrix
  Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor> energyMatrix(N, N);
  energyMatrix.reserve(Eigen::VectorXi::Constant(N, 4));

  // Mass matrix
  Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor> massMatrix(N, N);
  massMatrix.reserve(Eigen::VectorXi::Constant(N, 1));


  // === Build matrices

  // Build the mass matrix
  for (Face f : mesh.faces()) {
    size_t i = geometry.faceIndices[f];
    massMatrix.insert(i, i) = geometry.faceAreas[f];
  }

  // Build the energy matrix
  for (Face f : mesh.faces()) {
    size_t i = geometry.faceIndices[f];

    std::complex<double> weightISum = 0;
    for (Halfedge he : f.adjacentHalfedges()) {

      if (!he.twin().isInterior()) {
        continue;
      }

      Face neighFace = he.twin().face();
      unsigned int j = geometry.faceIndices[neighFace];

      // LC connection between the faces
      Vector2 rBar = std::pow(geometry.faceTransportCoefs[he.twin()], nSym);

      double weight = 1; // FIXME TODO figure out weights
      energyMatrix.insert(i, j) = -weight * rBar;
      weightISum += weight;
    }

    energyMatrix.insert(i, i) = weightISum;
  }

  // Shift to avoid singularity
  Eigen::SparseMatrix<Vector2> eye(N, N);
  eye.setIdentity();
  energyMatrix += 1e-4 * eye;

  // Store the solution here
  Eigen::VectorXcd solution;

  // If requested, align to principal curvatures
  if (alignCurvature) {

    Eigen::VectorXcd dirVec(N);
    for (Face f : mesh.faces()) {

      // Compute something like the principal directions
      double weightSum = 0;
      Vector2 sum = 0;

      for (Halfedge he : f.adjacentHalfedges()) {

        double dihedralAngle = std::abs(geometry.dihedralAngles[he.edge()]);
        double weight = norm(geometry->vector(he));
        weightSum += weight;
        double angleCoord = angleInPlane(geometry->vector(f.halfedge()), geometry->vector(he), geometry.faceNormals[f]);
        Vector2 coord = std::exp(angleCoord * IM_I *
                                 (double)nSym); // nsym should be 2 or 4, checked in the funciton which calls this

        sum += coord * weight * dihedralAngle;
      }

      sum /= weightSum;

      dirVec[geometry.faceIndices[f]] = sum;
    }

    // Normalize the alignment field
    double scale = std::sqrt(std::abs((dirVec.adjoint() * massMatrix * dirVec)[0]));
    dirVec /= scale;

    double lambdaT = 0.0; // this is something of a magical constant, see "Globally Optimal Direction Fields", eqn 16

    // Eigen::VectorXcd RHS = massMatrix * dirVec;
    Eigen::VectorXcd RHS = dirVec;
    Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor> LHS = energyMatrix - lambdaT * massMatrix;
    solution = solveSquare(LHS, RHS);

  }
  // Otherwise find the smallest eigenvector
  else {
    std::cout << "Solving smoothest field eigenvalue problem..." << std::endl;
    solution = smallestEigenvectorPositiveDefinite(energyMatrix, massMatrix);
  }


  // Copy the result to a FaceData object
  FaceData<Vector2> field(*mesh);
  for (Face f : mesh.faces()) {
    field[f] = solution[geometry.faceIndices[f]] / std::abs(solution[geometry.faceIndices[f]]);
  }

  return field;
}

FaceData<Vector2> computeSmoothestFaceDirectionField_boundary(IntrinsicGeometryInterface& geometry, int nSym = 1,
                                                              bool alignCurvature = false) {

  SurfaceMesh* mesh = geometry->getMesh();

  GeometryCache<Euclidean>& geometry = geometry->cache;
  geometry.requireFaceTransportCoefs();
  geometry.requireFaceNormals();
  geometry.requireFaceAreas();
  geometry.requireDihedralAngles();


  // Index interior faces
  size_t nInteriorFace = 0;
  FaceData<size_t> interiorFaceInd(*mesh, -77);
  FaceData<char> isInterior(*mesh);
  for (Face f : mesh.faces()) {
    bool isBoundary = false;
    for (Edge e : f.adjacentEdges()) {
      isBoundary |= e.isBoundary();
    }
    isInterior[f] = !isBoundary;
    if (!isBoundary) {
      interiorFaceInd[f] = nInteriorFace++;
    }
  }

  // Compute boundary values
  FaceData<Vector2> boundaryValues(*mesh);
  for (Face f : mesh.faces()) {
    if (isInterior[f]) {
      boundaryValues[f] = 0;
    } else {
      Vector3 bVec = Vector3::zero();
      for (Halfedge he : f.adjacentHalfedges()) {
        if (he.edge().isBoundary()) {
          bVec += geometry->vector(he).rotate_around(geometry.faceNormals[f], -PI / 2.0);
        }
      }
      Vector2 bC(dot(geometry.faceBases[f][0], bVec), dot(geometry.faceBases[f][1], bVec));
      bC = unit(bC);
      boundaryValues[f] = std::pow(bC, nSym);
    }
  }


  // === Allocate matrices
  // Energy matrix
  Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor> energyMatrix(nInteriorFace, nInteriorFace);
  energyMatrix.reserve(Eigen::VectorXi::Constant(nInteriorFace, 4));

  // Mass matrix
  Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor> massMatrix(nInteriorFace, nInteriorFace);
  massMatrix.reserve(Eigen::VectorXi::Constant(nInteriorFace, 1));

  // RHS
  Eigen::VectorXcd b(nInteriorFace);

  // === Build matrices

  // Build the mass matrix
  for (Face f : mesh.faces()) {
    if (isInterior[f]) {
      size_t i = interiorFaceInd[f];
      massMatrix.insert(i, i) = geometry.faceAreas[f];
    }
  }

  // Build the energy matrix
  for (Face f : mesh.faces()) {
    if (isInterior[f]) {
      size_t i = interiorFaceInd[f];

      std::complex<double> weightISum = 0;
      for (Halfedge he : f.adjacentHalfedges()) {

        Face neighFace = he.twin().face();
        double weight = 1; // FIXME TODO figure out weights
        Vector2 rBar = std::pow(geometry.faceTransportCoefs[he.twin()], nSym);

        if (isInterior[neighFace]) {
          size_t j = interiorFaceInd[neighFace];
          energyMatrix.insert(i, j) = -weight * rBar;
        } else {
          std::complex<double> bVal = boundaryValues[neighFace];
          b(i) += weight * rBar * bVal;
        }

        weightISum += weight;
      }

      energyMatrix.insert(i, i) = weightISum;
    }
  }

  // Shift to avoid singularity
  Eigen::SparseMatrix<Vector2> eye(nInteriorFace, nInteriorFace);
  eye.setIdentity();
  energyMatrix += 1e-4 * eye;

  // Store the solution here
  Eigen::VectorXcd solution;

  // If requested, align to principal curvatures
  if (alignCurvature) {

    Eigen::VectorXcd dirVec(nInteriorFace);
    for (Face f : mesh.faces()) {
      if (isInterior[f]) {

        // Compute something like the principal directions
        double weightSum = 0;
        Vector2 sum = 0;

        for (Halfedge he : f.adjacentHalfedges()) {

          double dihedralAngle = std::abs(geometry.dihedralAngles[he.edge()]);
          double weight = norm(geometry->vector(he));
          weightSum += weight;
          double angleCoord =
              angleInPlane(geometry->vector(f.halfedge()), geometry->vector(he), geometry.faceNormals[f]);
          Vector2 coord = std::exp(angleCoord * IM_I *
                                   (double)nSym); // nsym should be 2 or 4, checked in the funciton which calls this

          sum += coord * weight * dihedralAngle;
        }

        sum /= weightSum;

        // Normalize the curvature vectors. By doing so, we lose the property of adjusting the strength of the alignment
        // based on the strength of the curvature, but resolve any scaling issues between the magnitude of the normals
        // and the magnitude of the desired field.  Be careful when interpreting this as opposed to the usual direction
        // field optimization.
        dirVec[interiorFaceInd[f]] = unit(sum);
      }
    }


    double t = 0.1; // this is something of a magical constant, see "Globally
                    // Optimal Direction Fields", eqn 9
                    // NOTE: This value is different from the one used for vertex fields; seems to work better?

    std::cout << "Solving smoothest field dirichlet problem with curvature term..." << std::endl;
    Eigen::VectorXcd RHS = massMatrix * (t * dirVec + b);
    Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor> LHS = energyMatrix;
    solution = solveSquare(LHS, RHS);

  }
  // Otherwise find the general closest solution
  else {
    std::cout << "Solving smoothest field dirichlet problem..." << std::endl;
    Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor> LHS = energyMatrix;
    Eigen::VectorXcd RHS = massMatrix * b;
    solution = solveSquare(LHS, RHS);
  }


  // Copy the result to a FaceData object
  FaceData<Vector2> field(*mesh);
  for (Face f : mesh.faces()) {
    if (isInterior[f]) {
      field[f] = unit(solution[interiorFaceInd[f]]);
    } else {
      field[f] = unit(boundaryValues[f]);
    }
  }

  return field;
}

} // namespace

FaceData<Vector2> computeSmoothestFaceDirectionField(Geometry<Euclidean>* geometry, int nSym, bool alignCurvature) {

  std::cout << "Computing globally optimal direction field in faces" << std::endl;

  if (alignCurvature && !(nSym == 2 || nSym == 4)) {
    throw std::logic_error("ERROR: It only makes sense to align with curvature when nSym = 2 or "
                           "4");
  }

  // Dispatch to either the boundary of no boundary variant depending on the mesh type
  bool hasBoundary = false;
  for (Vertex v : geometry->getMesh()->vertices()) {
    hasBoundary |= v.isBoundary();
  }


  if (hasBoundary) {
    std::cout << "Mesh has boundary, computing dirichlet boundary condition solution" << std::endl;
    return computeSmoothestFaceDirectionField_boundary(geometry, nSym, alignCurvature);
  } else {
    std::cout << "Mesh has no boundary, computing unit-norm solution" << std::endl;
    return computeSmoothestFaceDirectionField_noBoundary(geometry, nSym, alignCurvature);
  }
}

FaceData<int> computeFaceIndex(IntrinsicGeometryInterface& geometry, const VertexData<Vector2>& directionField,
                               int nSym) {

  SurfaceMesh* mesh = geometry->getMesh();

  GeometryCache<Euclidean>& geometry = geometry->cache;
  geometry.requireFaceTransportCoefs();

  // Store the result here
  FaceData<int> indices(*mesh);

  // TODO haven't tested that this correctly reports the index when it is larger
  // than +-1

  for (Face f : mesh.faces()) {
    // Trace the direction field around the face and see how many times it
    // spins!
    double totalRot = 0;

    for (Halfedge he : f.adjacentHalfedges()) {
      // Compute the rotation along the halfedge implied by the field
      Vector2 x0 = directionField[he.vertex()];
      Vector2 x1 = directionField[he.twin().vertex()];
      Vector2 transport = std::pow(geometry.vertexTransportCoefs[he], nSym);

      // Find the difference in angle
      double theta0 = std::arg(transport * x0);
      double theta1 = std::arg(x1);
      double deltaTheta = regularizeAngle(theta1 - theta0 + PI) - PI; // regularize to [-PI,PI]

      totalRot += deltaTheta; // accumulate
    }

    // Compute the net rotation and corresponding index
    int index = static_cast<int>(std::round(totalRot / (2 * PI))); // should be very close to a multiple of 2PI
    indices[f] = index;
  }

  return indices;
}


VertexData<int> computeVertexIndex(IntrinsicGeometryInterface& geometry, const FaceData<Vector2>& directionField,
                                   int nSym) {

  SurfaceMesh* mesh = geometry->getMesh();
  GeometryCache<Euclidean>& geometry = geometry->cache;
  geometry.requireFaceTransportCoefs();

  // Store the result here
  VertexData<int> indices(*mesh);

  // TODO haven't tested that this correctly reports the index when it is larger
  // than +-1

  for (Vertex v : mesh.vertices()) {

    // Trace the direction field around the face and see how many times it
    // spins!
    double totalRot = 0;

    for (Halfedge he : v.incomingHalfedges()) {
      // Compute the rotation along the halfedge implied by the field
      Vector2 x0 = directionField[he.face()];
      Vector2 x1 = directionField[he.twin().face()];
      Vector2 transport = std::pow(geometry.faceTransportCoefs[he], nSym);

      // Find the difference in angle
      double theta0 = std::arg(transport * x0);
      double theta1 = std::arg(x1);
      double deltaTheta = std::arg(x1 / (transport * x0));

      totalRot += deltaTheta;
    }

    double angleDefect = geometry->angleDefect(v);
    totalRot += angleDefect * nSym;

    // Compute the net rotation and corresponding index
    int index = static_cast<int>(std::round(totalRot / (2 * PI))); // should be very close to a multiple of 2PI
    indices[v] = index;
  }

  return indices;
}
*/

} // namespace surface
} // namespace geometrycentral
