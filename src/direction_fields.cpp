#include <geometrycentral/direction_fields.h>

#include "geometrycentral/linear_solvers.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>

using std::cout;
using std::endl;

namespace geometrycentral {

// Anonymous namespace for helper functions
namespace {

VertexData<Complex> computeSmoothestVertexDirectionField_noBoundary(Geometry<Euclidean>* geometry, int nSym,
                                                                    bool alignCurvature) {
  HalfedgeMesh* mesh = geometry->getMesh();
  size_t N = mesh->nVertices();

  GeometryCache<Euclidean>& gc = geometry->cache;
  gc.requireVertexTransportCoefs();
  gc.requireEdgeCotanWeights();
  gc.requireVertexIndices();
  gc.requireVertexDualAreas();

  // Energy matrix
  Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor> energyMatrix(
      N, N); // have to use ColMajor because LU solver below demands it

  // Supposedly reserving space in the matrix makes construction real zippy
  // below
  Eigen::VectorXi nEntries(N);
  for (VertexPtr v : mesh->vertices()) {
    nEntries[gc.vertexIndices[v]] = v.degree() + 1;
  }
  energyMatrix.reserve(nEntries);

  // Mass matrix
  Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor> massMatrix(N, N);
  massMatrix.reserve(1);

  // === Build matrices

  // Build the mass matrix
  for (VertexPtr v : mesh->vertices()) {
    size_t i = gc.vertexIndices[v];
    massMatrix.insert(i, i) = gc.vertexDualAreas[v];
  }

  // Build the energy matrix
  for (VertexPtr v : mesh->vertices()) {
    size_t i = gc.vertexIndices[v];

    std::complex<double> weightISum = 0;
    for (HalfedgePtr he : v.incomingHalfedges()) {
      size_t j = gc.vertexIndices[he.vertex()];
      std::complex<double> rBar = std::pow(gc.vertexTransportCoefs[he], nSym);
      double weight = gc.edgeCotanWeights[he.edge()];
      energyMatrix.insert(i, j) = -weight * rBar;
      weightISum += weight;
    }

    energyMatrix.insert(i, i) = weightISum;
  }

  // Shift to avoid singularity
  Eigen::SparseMatrix<Complex> eye(N, N);
  eye.setIdentity();
  energyMatrix += 1e-4 * eye;

  // Store the solution here
  Eigen::VectorXcd solution;

  // If requested, align to principal curvatures
  if (alignCurvature) {

    gc.requirePrincipalDirections();

    Eigen::VectorXcd dirVec(N);
    if (nSym == 2) {
      for (VertexPtr v : mesh->vertices()) {
        dirVec[gc.vertexIndices[v]] = gc.principalDirections[v];
      }
    } else if (nSym == 4) {
      for (VertexPtr v : mesh->vertices()) {
        dirVec[gc.vertexIndices[v]] = std::pow(gc.principalDirections[v], 2);
      }
    }

    // Normalize the alignment field
    double scale = std::sqrt(std::abs((dirVec.adjoint() * massMatrix * dirVec)[0]));
    dirVec /= scale;

    double lambdaT = 0.0; // this is something of a magical constant, see
                          // "Globally Optimal Direction Fields", eqn 16

    // Eigen::VectorXcd RHS = massMatrix * dirVec;
    Eigen::VectorXcd RHS = dirVec;
    Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor> LHS = energyMatrix - lambdaT * massMatrix;
    solution = solveSquare(LHS, RHS);
  }
  // Otherwise find the smallest eigenvector
  else {
    std::cout << "Solving smoothest field eigenvalue problem..." << std::endl;
    solution = smallestEigenvectorSquare(energyMatrix, massMatrix);
  }

  // Copy the result to a VertexData vector
  VertexData<Complex> toReturn(mesh);
  for (VertexPtr v : mesh->vertices()) {
    toReturn[v] = solution[gc.vertexIndices[v]] / std::abs(solution[gc.vertexIndices[v]]);
  }

  return toReturn;
}

VertexData<Complex> computeSmoothestVertexDirectionField_boundary(Geometry<Euclidean>* geometry, int nSym,
                                                                  bool alignCurvature) {
  HalfedgeMesh* mesh = geometry->getMesh();
  size_t nInterior = mesh->nInteriorVertices();

  GeometryCache<Euclidean>& gc = geometry->cache;
  gc.requireVertexTransportCoefs();
  gc.requireEdgeCotanWeights();
  gc.requireVertexBases();
  gc.requireInteriorVertexIndices();
  gc.requireVertexDualAreas();

  // Compute the boundary values
  VertexData<std::complex<double>> boundaryValues(mesh);
  for (VertexPtr v : mesh->vertices()) {
    if (v.isBoundary()) {
      Vector3 b = geometry->boundaryNormal(v);
      Complex bC(dot(gc.vertexBases[v][0], b), dot(gc.vertexBases[v][1], b)); // TODO can do better
      bC = unit(bC);
      boundaryValues[v] = std::pow(bC, nSym);
    } else {
      boundaryValues[v] = 0;
    }
  }

  VertexData<size_t> vertInd = mesh->getInteriorVertexIndices();

  // Energy matrix
  Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor> energyMatrix(nInterior, nInterior);

  Eigen::VectorXi nEntries(nInterior);
  for (VertexPtr v : mesh->vertices()) {
    if (v.isBoundary()) {
      continue;
    }
    nEntries[gc.interiorVertexIndices[v]] = v.degree() + 1;
  }
  energyMatrix.reserve(nEntries);

  // Mass matrix
  Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor> massMatrix(nInterior, nInterior);
  massMatrix.reserve(1);

  // RHS
  Eigen::VectorXcd b(nInterior);

  // === Build matrices

  // Build the mass matrix and zero b
  for (VertexPtr v : mesh->vertices()) {
    if (v.isBoundary()) {
      continue;
    }
    size_t i = gc.interiorVertexIndices[v];
    b(i) = 0.0;
    massMatrix.insert(i, i) = gc.vertexDualAreas[v];
  }

  // Build the energy matrix
  for (VertexPtr v : mesh->vertices()) {
    if (v.isBoundary()) {
      continue;
    }
    size_t i = gc.interiorVertexIndices[v];

    std::complex<double> weightISum = 0;
    for (HalfedgePtr he : v.incomingHalfedges()) {
      std::complex<double> rBar = std::pow(gc.vertexTransportCoefs[he], nSym);
      double w = gc.edgeCotanWeights[he.edge()];

      // Interior-boundary term
      if (he.vertex().isBoundary()) {
        std::complex<double> bVal = boundaryValues[he.vertex()];
        b(i) += w * rBar * bVal;
      } else { // Interior-interior term
        size_t j = gc.interiorVertexIndices[he.vertex()];
        energyMatrix.insert(i, j) = -w * rBar;
      }
      weightISum += w;
    }

    energyMatrix.insert(i, i) = weightISum;
  }

  // Shift to avoid singularities
  Eigen::SparseMatrix<Complex> eye(nInterior, nInterior);
  eye.setIdentity();
  energyMatrix += 1e-4 * eye;

  // Compute the actual solution
  std::cout << "Solving linear problem..." << std::endl;

  // Store the solution here
  Eigen::VectorXcd solution;

  // If requested, align to principal curvatures
  if (alignCurvature) {

    gc.requirePrincipalDirections();

    Eigen::VectorXcd dirVec(nInterior);
    for (VertexPtr v : mesh->vertices()) {
      if (v.isBoundary()) {
        continue;
      }

      Complex directionVal = gc.principalDirections[v];
      if (nSym == 4) {
        directionVal = std::pow(directionVal, 2);
      }

      // Normalize the curvature vectors. By doing so, we lose the property of adjusting the strength of the alignment
      // based on the strength of the curvature, but resolve any scaling issues between the magnitude of the normals and
      // the magnitude of the desired field.  Be careful when interpreting this as opposed to the usual direction field
      // optimization.
      dirVec[gc.interiorVertexIndices[v]] = directionVal / std::abs(directionVal);
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
  VertexData<Complex> toReturn(mesh);
  for (VertexPtr v : mesh->vertices()) {
    if (v.isBoundary()) {
      toReturn[v] = boundaryValues[v];
    } else {
      toReturn[v] = unit(solution[gc.interiorVertexIndices[v]]);
    }
  }

  return toReturn;
}
}; // namespace

VertexData<Complex> computeSmoothestVertexDirectionField(Geometry<Euclidean>* geometry, int nSym, bool alignCurvature) {
  std::cout << "Computing globally optimal direction field" << std::endl;

  if (alignCurvature && !(nSym == 2 || nSym == 4)) {
    throw std::logic_error("ERROR: It only makes sense to align with curvature when nSym = 2 or "
                           "4");
  }

  // Dispatch to either the boundary of no boundary variant depending on the
  // mesh type
  bool hasBoundary = false;
  for (VertexPtr v : geometry->getMesh()->vertices()) {
    hasBoundary |= v.isBoundary();
  }

  if (hasBoundary) {
    std::cout << "Mesh has boundary, computing dirichlet boundary condition solution" << std::endl;
    return computeSmoothestVertexDirectionField_boundary(geometry, nSym, alignCurvature);
  } else {
    std::cout << "Mesh has no boundary, computing unit-norm solution" << std::endl;
    return computeSmoothestVertexDirectionField_noBoundary(geometry, nSym, alignCurvature);
  }
}

// Helpers for computing face-based direction fields
namespace {

FaceData<Complex> computeSmoothestFaceDirectionField_noBoundary(Geometry<Euclidean>* geometry, int nSym,
                                                                bool alignCurvature) {

  HalfedgeMesh* mesh = geometry->getMesh();
  unsigned int N = mesh->nFaces();

  GeometryCache<Euclidean>& gc = geometry->cache;
  gc.requireFaceTransportCoefs();
  gc.requireFaceNormals();
  gc.requireFaceAreas();
  gc.requireDihedralAngles();
  gc.requireFaceIndices();

  // === Allocate matrices
  // Energy matrix
  Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor> energyMatrix(N, N);
  energyMatrix.reserve(Eigen::VectorXi::Constant(N, 4));

  // Mass matrix
  Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor> massMatrix(N, N);
  massMatrix.reserve(Eigen::VectorXi::Constant(N, 1));


  // === Build matrices

  // Build the mass matrix
  for (FacePtr f : mesh->faces()) {
    size_t i = gc.faceIndices[f];
    massMatrix.insert(i, i) = gc.faceAreas[f];
  }

  // Build the energy matrix
  for (FacePtr f : mesh->faces()) {
    size_t i = gc.faceIndices[f];

    std::complex<double> weightISum = 0;
    for (HalfedgePtr he : f.adjacentHalfedges()) {

      if (!he.twin().isReal()) {
        continue;
      }

      FacePtr neighFace = he.twin().face();
      unsigned int j = gc.faceIndices[neighFace];

      // LC connection between the faces
      Complex rBar = std::pow(gc.faceTransportCoefs[he.twin()], nSym);

      double weight = 1; // FIXME TODO figure out weights
      energyMatrix.insert(i, j) = -weight * rBar;
      weightISum += weight;
    }

    energyMatrix.insert(i, i) = weightISum;
  }

  // Shift to avoid singularity
  Eigen::SparseMatrix<Complex> eye(N, N);
  eye.setIdentity();
  energyMatrix += 1e-4 * eye;

  // Store the solution here
  Eigen::VectorXcd solution;

  // If requested, align to principal curvatures
  if (alignCurvature) {

    Eigen::VectorXcd dirVec(N);
    for (FacePtr f : mesh->faces()) {

      // Compute something like the principal directions
      double weightSum = 0;
      Complex sum = 0;

      for (HalfedgePtr he : f.adjacentHalfedges()) {

        double dihedralAngle = std::abs(gc.dihedralAngles[he.edge()]);
        double weight = norm(geometry->vector(he));
        weightSum += weight;
        double angleCoord = angleInPlane(geometry->vector(f.halfedge()), geometry->vector(he), gc.faceNormals[f]);
        Complex coord = std::exp(angleCoord * IM_I *
                                 (double)nSym); // nsym should be 2 or 4, checked in the funciton which calls this

        sum += coord * weight * dihedralAngle;
      }

      sum /= weightSum;

      dirVec[gc.faceIndices[f]] = sum;
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
  FaceData<Complex> field(mesh);
  for (FacePtr f : mesh->faces()) {
    field[f] = solution[gc.faceIndices[f]] / std::abs(solution[gc.faceIndices[f]]);
  }

  return field;
}

FaceData<Complex> computeSmoothestFaceDirectionField_boundary(Geometry<Euclidean>* geometry, int nSym,
                                                              bool alignCurvature) {

  HalfedgeMesh* mesh = geometry->getMesh();

  GeometryCache<Euclidean>& gc = geometry->cache;
  gc.requireFaceTransportCoefs();
  gc.requireFaceNormals();
  gc.requireFaceAreas();
  gc.requireDihedralAngles();


  // Index interior faces
  size_t nInteriorFace = 0;
  FaceData<size_t> interiorFaceInd(mesh, -77);
  FaceData<char> isInterior(mesh);
  for (FacePtr f : mesh->faces()) {
    bool isBoundary = false;
    for (EdgePtr e : f.adjacentEdges()) {
      isBoundary |= e.isBoundary();
    }
    isInterior[f] = !isBoundary;
    if (!isBoundary) {
      interiorFaceInd[f] = nInteriorFace++;
    }
  }

  // Compute boundary values
  FaceData<Complex> boundaryValues(mesh);
  for (FacePtr f : mesh->faces()) {
    if (isInterior[f]) {
      boundaryValues[f] = 0;
    } else {
      Vector3 bVec = Vector3::zero();
      for (HalfedgePtr he : f.adjacentHalfedges()) {
        if (he.edge().isBoundary()) {
          bVec += geometry->vector(he).rotate_around(gc.faceNormals[f], -PI / 2.0);
        }
      }
      Complex bC(dot(gc.faceBases[f][0], bVec), dot(gc.faceBases[f][1], bVec));
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
  for (FacePtr f : mesh->faces()) {
    if (isInterior[f]) {
      size_t i = interiorFaceInd[f];
      massMatrix.insert(i, i) = gc.faceAreas[f];
    }
  }

  // Build the energy matrix
  for (FacePtr f : mesh->faces()) {
    if (isInterior[f]) {
      size_t i = interiorFaceInd[f];

      std::complex<double> weightISum = 0;
      for (HalfedgePtr he : f.adjacentHalfedges()) {

        FacePtr neighFace = he.twin().face();
        double weight = 1; // FIXME TODO figure out weights
        Complex rBar = std::pow(gc.faceTransportCoefs[he.twin()], nSym);

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
  Eigen::SparseMatrix<Complex> eye(nInteriorFace, nInteriorFace);
  eye.setIdentity();
  energyMatrix += 1e-4 * eye;

  // Store the solution here
  Eigen::VectorXcd solution;

  // If requested, align to principal curvatures
  if (alignCurvature) {

    Eigen::VectorXcd dirVec(nInteriorFace);
    for (FacePtr f : mesh->faces()) {
      if (isInterior[f]) {

        // Compute something like the principal directions
        double weightSum = 0;
        Complex sum = 0;

        for (HalfedgePtr he : f.adjacentHalfedges()) {

          double dihedralAngle = std::abs(gc.dihedralAngles[he.edge()]);
          double weight = norm(geometry->vector(he));
          weightSum += weight;
          double angleCoord = angleInPlane(geometry->vector(f.halfedge()), geometry->vector(he), gc.faceNormals[f]);
          Complex coord = std::exp(angleCoord * IM_I *
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


    double t = 0.1;  // this is something of a magical constant, see "Globally
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
  FaceData<Complex> field(mesh);
  for (FacePtr f : mesh->faces()) {
    if (isInterior[f]) {
      field[f] = unit(solution[interiorFaceInd[f]]);
    } else {
      field[f] = unit(boundaryValues[f]);
    }
  }

  return field;
}

} // namespace

FaceData<Complex> computeSmoothestFaceDirectionField(Geometry<Euclidean>* geometry, int nSym, bool alignCurvature) {

  std::cout << "Computing globally optimal direction field in faces" << std::endl;

  if (alignCurvature && !(nSym == 2 || nSym == 4)) {
    throw std::logic_error("ERROR: It only makes sense to align with curvature when nSym = 2 or "
                           "4");
  }

  // Dispatch to either the boundary of no boundary variant depending on the mesh type
  bool hasBoundary = false;
  for (VertexPtr v : geometry->getMesh()->vertices()) {
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


FaceData<int> computeFaceIndex(Geometry<Euclidean>* geometry, VertexData<Complex> directionField, int nSym) {
  HalfedgeMesh* mesh = geometry->getMesh();

  GeometryCache<Euclidean>& gc = geometry->cache;
  gc.requireFaceTransportCoefs();

  // Store the result here
  FaceData<int> indices(mesh);

  // TODO haven't tested that this correctly reports the index when it is larger
  // than +-1

  for (FacePtr f : mesh->faces()) {
    // Trace the direction field around the face and see how many times it
    // spins!
    double totalRot = 0;

    for (HalfedgePtr he : f.adjacentHalfedges()) {
      // Compute the rotation along the halfedge implied by the field
      Complex x0 = directionField[he.vertex()];
      Complex x1 = directionField[he.twin().vertex()];
      Complex transport = std::pow(gc.vertexTransportCoefs[he], nSym);

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


VertexData<int> computeVertexIndex(Geometry<Euclidean>* geometry, FaceData<Complex> directionField, int nSym) {

  HalfedgeMesh* mesh = geometry->getMesh();
  GeometryCache<Euclidean>& gc = geometry->cache;
  gc.requireFaceTransportCoefs();

  // Store the result here
  VertexData<int> indices(mesh);

  // TODO haven't tested that this correctly reports the index when it is larger
  // than +-1

  for (VertexPtr v : mesh->vertices()) {

    // Trace the direction field around the face and see how many times it
    // spins!
    double totalRot = 0;

    for (HalfedgePtr he : v.incomingHalfedges()) {
      // Compute the rotation along the halfedge implied by the field
      Complex x0 = directionField[he.face()];
      Complex x1 = directionField[he.twin().face()];
      Complex transport = std::pow(gc.faceTransportCoefs[he], nSym);

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

} // namespace geometrycentral
