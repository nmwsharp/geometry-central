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

VertexData<double> computeAngleDefects(Geometry<Euclidean>* geometry) {
  HalfedgeMesh* mesh = geometry->getMesh();

  VertexData<double> defect(mesh);

  for (VertexPtr v : mesh->vertices()) {
    // Sum up the wedge angle for adjacent triangles
    double angleSum = 0.0;
    HalfedgePtr prevHe{nullptr};
    for (HalfedgePtr he : v.incomingHalfedges()) {
      if (prevHe != nullptr) {
        Vector3 a = unit(-geometry->vector(he));
        Vector3 b = unit(-geometry->vector(prevHe));

        double angle = acos(clamp(dot(a, b), -1.0, 1.0)); // clamp to deal with numerical error

        // acos() can choose the wrong inverse for very sharp angles (only
        // happens with boundary vertices, assuming no triangles are degenerate)
        if (dot(cross(a, b), geometry->normal(v)) < 0) {
          angle = 2 * PI - angle;
        }

        angleSum += angle;
      }

      prevHe = he;
    }

    // Do the last one, which we missed above
    Vector3 a = unit(-geometry->vector(v.halfedge().twin()));
    Vector3 b = unit(-geometry->vector(prevHe));
    double angle = acos(clamp(dot(a, b), -1.0, 1.0)); // clamp to deal with numerical error
    if (dot(cross(a, b), geometry->normal(v)) < 0) {
      angle = 2 * PI - angle;
    }
    angleSum += angle;

    defect[v] = 2 * PI - angleSum;
  }

  return defect;
}

HalfedgeData<double> computeRescaledHalfedgeAngles(Geometry<Euclidean>* geometry) {
  HalfedgeMesh* mesh = geometry->getMesh();
  VertexData<double> defects = computeAngleDefects(geometry);
  return computeRescaledHalfedgeAngles(geometry, defects);
}

HalfedgeData<double> computeRescaledHalfedgeAngles(Geometry<Euclidean>* geometry,
                                                   const VertexData<double>& angleDefects) {
  HalfedgeMesh* mesh = geometry->getMesh();

  HalfedgeData<double> edgeAngles(mesh);

  for (VertexPtr v : mesh->vertices()) {
    double cumulativeAngleSum = 0.0;
    double totalAngleSum = 2 * PI - angleDefects[v];

    // Sum up the wedge angle for adjacent triangles
    double angleSum = 0.0;
    HalfedgePtr prevHe{nullptr};
    for (HalfedgePtr he : v.incomingHalfedges()) {
      if (prevHe == nullptr) {
        // Angle to the reference edge is 0 by definition
        edgeAngles[he] = 0.0;

      } else {
        Vector3 a = unit(-geometry->vector(he));
        Vector3 b = unit(-geometry->vector(prevHe));

        double angle = acos(clamp(dot(a, b), -1.0, 1.0)); // clamp to deal with numerical error

        // acos() can choose the wrong inverse for very sharp angles (only
        // happens with boundary vertices)
        if (dot(cross(a, b), geometry->normal(v)) < 0) {
          angle = 2 * PI - angle;
        }
        cumulativeAngleSum += angle;
        edgeAngles[he] = 2 * PI * (1.0 - cumulativeAngleSum / totalAngleSum); // subtract because we
                                                                              // measure angles in
                                                                              // the wrong direction
                                                                              // with this loop
      }

      prevHe = he;
    }
  }

  return edgeAngles;
}

HalfedgeData<double> computeTransportAngles(Geometry<Euclidean>* geometry) {
  HalfedgeMesh* mesh = geometry->getMesh();
  HalfedgeData<double> edgeAngles = computeRescaledHalfedgeAngles(geometry);
  return computeTransportAngles(geometry, edgeAngles);
}

HalfedgeData<double> computeTransportAngles(Geometry<Euclidean>* geometry,
                                            const HalfedgeData<double>& rescaledHalfedgeAngles) {
  HalfedgeMesh* mesh = geometry->getMesh();

  HalfedgeData<double> transportAngles(mesh);

  for (HalfedgePtr he : mesh->halfedges()) {
    transportAngles[he] = regularizeAngle(rescaledHalfedgeAngles[he] - rescaledHalfedgeAngles[he.twin()] + PI);
  }
  for (HalfedgePtr he : mesh->imaginaryHalfedges()) {
    transportAngles[he] = regularizeAngle(rescaledHalfedgeAngles[he] - rescaledHalfedgeAngles[he.twin()] + PI);
  }

  return transportAngles;
}

VertexData<Vector3> convertTangentAnglesToR3Vectors(Geometry<Euclidean>* geometry,
                                                    const VertexData<double>& tangentAngles) {
  HalfedgeMesh* mesh = geometry->getMesh();

  VertexData<Vector3> r3Vectors(mesh);

  for (VertexPtr v : mesh->vertices()) {
    // Project the reference halfedge in to the tangent plane
    Vector3 r = geometry->vector(v.halfedge());
    Vector3 N = geometry->normal(v);
    Vector3 referenceInTangentSpace = unit(r - dot(r, N) * N);

    // Rotate the reference vector to get the desired direction
    r3Vectors[v] = referenceInTangentSpace.rotate_around(N, tangentAngles[v]);
  }

  return r3Vectors;
}

VertexData<Vector3> convertComplexDirectionsToR3Vectors(Geometry<Euclidean>* geometry,
                                                        const VertexData<Complex>& directionField, int nSym) {
  HalfedgeMesh* mesh = geometry->getMesh();

  VertexData<Vector3> r3Vectors(mesh);
  for (VertexPtr v : mesh->vertices()) {
    // Compute the implied tangent angles
    std::complex<double> solVal = std::pow(directionField[v], 1.0 / nSym);
    double tangentAngle = std::arg(solVal);

    // Project the reference halfedge in to the tangent plane
    Vector3 r = geometry->vector(v.halfedge());
    Vector3 N = geometry->normal(v);
    Vector3 referenceInTangentSpace = unit(r - dot(r, N) * N);

    // Rotate the reference vector to get the desired direction
    r3Vectors[v] = referenceInTangentSpace.rotate_around(N, tangentAngle) * std::abs(solVal); // preserve scale
  }

  return r3Vectors;
}

// Anonymous namespace for helper functions
namespace {

VertexData<Complex> computeSmoothestVertexDirectionField_noBoundary(Geometry<Euclidean>* geometry, int nSym,
                                                                    bool alignCurvature) {
  HalfedgeMesh* mesh = geometry->getMesh();

  size_t N = mesh->nVertices();

  // Geometric values needed to build matrices
  HalfedgeData<double> transportAngles = computeTransportAngles(geometry);

  // === Allocate matrices
  VertexData<size_t> vertInd = mesh->getVertexIndices();

  // Energy matrix
  Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor> energyMatrix(
      N, N); // have to use ColMajor because LU solver below demands it
  // Supposedly reserving space in the matrix makes construction real zippy
  // below
  Eigen::VectorXi nEntries(N);
  for (VertexPtr v : mesh->vertices()) {
    nEntries[vertInd[v]] = v.degree() + 1;
  }
  energyMatrix.reserve(nEntries);

  // Mass matrix
  Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor> massMatrix(N, N);
  massMatrix.reserve(1);

  // === Build matrices

  // Build the mass matrix
  for (VertexPtr v : mesh->vertices()) {
    size_t i = vertInd[v];
    massMatrix.insert(i, i) = geometry->dualArea(v);
  }

  // Build the energy matrix
  for (VertexPtr v : mesh->vertices()) {
    size_t i = vertInd[v];

    std::complex<double> weightISum = 0;
    for (HalfedgePtr he : v.incomingHalfedges()) {
      size_t j = vertInd[he.vertex()];

      std::complex<double> rBar = std::exp(-nSym * transportAngles[he.twin()] * IM_I);

      double weight = geometry->cotanWeight(he.edge());
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
    VertexData<Complex> principalDirections;
    geometry->getPrincipalDirections(principalDirections);
    Eigen::VectorXcd dirVec(N);
    if (nSym == 2) {
      for (VertexPtr v : mesh->vertices()) {
        dirVec[vertInd[v]] = principalDirections[v];
      }
    } else if (nSym == 4) {
      for (VertexPtr v : mesh->vertices()) {
        dirVec[vertInd[v]] = std::pow(principalDirections[v], 2);
        // dirVec[vertInd[v]] = std::pow(principalDirections[v], 2) *
        // std::exp(PI * IM_I);
        // dirVec[vertInd[v]] /= std::abs(dirVec[vertInd[v]]);
        // cout << "princ dir = " << dirVec[vertInd[v]] << endl;
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
    solution = SquareSolver<Complex>::solve(LHS, RHS);
  }
  // Otherwise find the smallest eigenvector
  else {
    std::cout << "Solving smoothest field eigenvalue problem..." << std::endl;
    solution = smallestEigenvectorPositiveDefinite(energyMatrix, massMatrix);
  }

  // Copy the result to a VertexData vector
  VertexData<Complex> toReturn(mesh);
  for (VertexPtr v : mesh->vertices()) {
    // cout << "sol = " << solution[vertInd[v]] << endl;
    toReturn[v] = solution[vertInd[v]] / std::abs(solution[vertInd[v]]);
  }

  return toReturn;
}

VertexData<Complex> computeSmoothestVertexDirectionField_boundary(Geometry<Euclidean>* geometry, int nSym,
                                                                  bool alignCurvature) {
  HalfedgeMesh* mesh = geometry->getMesh();

  size_t nInterior = mesh->nInteriorVertices();

  // Compute the boundary values
  VertexData<std::complex<double>> boundaryValues(mesh);
  for (VertexPtr v : mesh->vertices()) {
    if (v.isBoundary()) {
      Vector3 b = geometry->boundaryNormal(v);
      boundaryValues[v] = std::pow(geometry->tangentVectorToComplexAngle(v, b), nSym);
    } else {
      boundaryValues[v] = 0;
    }
  }

  // Geometric values needed to build matrices
  HalfedgeData<double> transportAngles = computeTransportAngles(geometry);

  // === Allocate matrices
  VertexData<size_t> vertInd = mesh->getInteriorVertexIndices();

  // Energy matrix
  Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor> energyMatrix(
      nInterior,
      nInterior); // have to use ColMajor because LU solver below demands it
  // Supposedly reserving space in the matrix makes construction real zippy
  // below
  Eigen::VectorXi nEntries(nInterior);
  for (VertexPtr v : mesh->vertices()) {
    if (v.isBoundary()) {
      continue;
    }
    nEntries[vertInd[v]] = v.degree() + 1;
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
    size_t i = vertInd[v];
    b(i) = 0.0;
    massMatrix.insert(i, i) = geometry->dualArea(v);
  }

  // Build the energy matrix
  for (VertexPtr v : mesh->vertices()) {
    if (v.isBoundary()) {
      continue;
    }
    size_t i = vertInd[v];

    std::complex<double> weightISum = 0;
    for (HalfedgePtr he : v.incomingHalfedges()) {
      std::complex<double> rBar = std::exp(-nSym * transportAngles[he.twin()] * IM_I);
      double w = geometry->cotanWeight(he.edge());

      // Interior-boundary term
      if (he.vertex().isBoundary()) {
        std::complex<double> bVal = boundaryValues[he.vertex()];
        b(i) += w * rBar * bVal;
      } else { // Interior-interior term

        size_t j = vertInd[he.vertex()];
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
    VertexData<Complex> principalDirections;
    geometry->getPrincipalDirections(principalDirections);
    Eigen::VectorXcd dirVec(nInterior);
    for (VertexPtr v : mesh->vertices()) {
      if (v.isBoundary()) {
        continue;
      }

      Complex directionVal = principalDirections[v];
      if (nSym == 4) {
        directionVal = std::pow(directionVal, 2);
      }

      // Normalize the curvature vectors. By doing so, we lose the property of
      // adjusting the strength
      // of the alignment based on the strength of the curvature, but resolve
      // any scaling issues
      // between the magnitude of the normals and the magnitude of the desired
      // field.
      // Be careful when interpreting this as opposed to the usual direction
      // field optimization.
      dirVec[vertInd[v]] = directionVal / std::abs(directionVal);
    }

    double t = 0.01; // this is something of a magical constant, see "Globally
                     // Optimal Direction Fields", eqn 9

    Eigen::VectorXcd RHS = massMatrix * (t * dirVec + b);
    Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor> LHS = massMatrix * energyMatrix;
    solution = SquareSolver<Complex>::solve(LHS, RHS);
  }
  // Otherwise find the general closest solution
  else {
    std::cout << "Solving smoothest field dirichlet problem..." << std::endl;

    Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor> LHS = massMatrix * energyMatrix; // can be simplified
    Eigen::VectorXcd RHS = massMatrix * b;
    solution = SquareSolver<Complex>::solve(LHS, RHS);
  }

  // Copy the result to a VertexData vector for both the boudary and interior
  VertexData<Complex> toReturn(mesh);
  for (VertexPtr v : mesh->vertices()) {
    if (v.isBoundary()) {
      toReturn[v] = boundaryValues[v];
    } else {
      toReturn[v] = solution[vertInd[v]] / std::abs(solution[vertInd[v]]);
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
  gc.faceTransportCoefsQ.require();
  gc.faceNormalsQ.require();
  gc.dihedralAngleQ.require();

  // === Allocate matrices
  FaceData<size_t> faceInd = mesh->getFaceIndices();

  // Energy matrix
  Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor> energyMatrix(
      N, N); // use ColMajor because LU solver below demands it
  energyMatrix.reserve(
      Eigen::VectorXi::Constant(N, 4)); // Reserving space in the matrix makes construction real zippy below

  // Mass matrix
  Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor> massMatrix(N, N);
  massMatrix.reserve(Eigen::VectorXi::Constant(N, 1));


  // === Build matrices

  // Build the mass matrix
  for (FacePtr f : mesh->faces()) {
    size_t i = faceInd[f];
    massMatrix.insert(i, i) = geometry->area(f);
  }

  // Build the energy matrix
  for (FacePtr f : mesh->faces()) {
    size_t i = faceInd[f];

    std::complex<double> weightISum = 0;
    for (HalfedgePtr he : f.adjacentHalfedges()) {

      if (!he.twin().isReal()) {
        continue;
      }

      FacePtr neighFace = he.twin().face();
      unsigned int j = faceInd[neighFace];

      // LC connection between the faces
      Complex rBar = std::pow(gc.faceTransportCoefs[he.twin()], nSym);

      double weight = 1; // TODO figure out weights
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

        double dihedralAngle = std::abs(gc.dihedralAngle[he.edge()]);
        double weight = norm(geometry->vector(he));
        weightSum += weight;
        double angleCoord = angleInPlane(geometry->vector(f.halfedge()), geometry->vector(he), gc.faceNormals[f]);
        Complex coord = std::exp(angleCoord * IM_I *
                                 (double)nSym); // nsym should be 2 or 4, checked in the funciton which calls this

        sum += coord * weight * dihedralAngle;
      }

      sum /= weightSum;

      dirVec[faceInd[f]] = sum;
    }

    // Normalize the alignment field
    double scale = std::sqrt(std::abs((dirVec.adjoint() * massMatrix * dirVec)[0]));
    dirVec /= scale;

    double lambdaT = 0.0; // this is something of a magical constant, see "Globally Optimal Direction Fields", eqn 16

    // Eigen::VectorXcd RHS = massMatrix * dirVec;
    Eigen::VectorXcd RHS = dirVec;
    Eigen::SparseMatrix<std::complex<double>, Eigen::ColMajor> LHS = energyMatrix - lambdaT * massMatrix;
    solution = SquareSolver<Complex>::solve(LHS, RHS);

  }
  // Otherwise find the smallest eigenvector
  else {
    std::cout << "Solving smoothest field eigenvalue problem..." << std::endl;
    solution = smallestEigenvectorPositiveDefinite(energyMatrix, massMatrix);
  }


  // Copy the result to a FaceData object
  FaceData<Complex> field(mesh);
  for (FacePtr f : mesh->faces()) {
    field[f] = solution[faceInd[f]] / std::abs(solution[faceInd[f]]);
  }

  return field;
}

FaceData<Complex> computeSmoothestFaceDirectionField_boundary(Geometry<Euclidean>* geometry, int nSym,
                                                              bool alignCurvature) {

  // TODO implement
  HalfedgeMesh* mesh = geometry->getMesh();
  FaceData<Complex> field(mesh);

  throw std::runtime_error("Boundary fields not implemented");

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

  // Precompute transport angles
  // Note: These tend to get computed a bunch of times when working with
  // direction fields
  HalfedgeData<double> transportAngles = computeTransportAngles(geometry);

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
      Complex transport = std::exp(nSym * transportAngles[he] * IM_I);

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
  gc.faceTransportCoefsQ.require();

  // Store the result here
  VertexData<int> indices(mesh);

  // TODO haven't tested that this correctly reports the index when it is larger
  // than +-1

  for (VertexPtr v : mesh->vertices()) {

    // bool interesting = v.degree() == 3;
    // if (interesting) {
    //   cout << "!!! interesting !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    // }

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
      double deltaTheta = theta1 - theta0;

      totalRot += deltaTheta; // accumulate

      // if (interesting) {
      //   cout << "he: " << endl;
      //   cout << "  theta0 = " << (theta0 / (2 * PI)) << " theta1 = " << (theta1 / (2 * PI)) << endl;
      //   cout << "  deltaTheta/2PI = " << (deltaTheta / (2 * PI)) << endl;
      // }
    }

    // if (interesting) {
    //   cout << "totalRot (before reg) /2PI = " << (totalRot / (2 * PI)) << endl;
    // }
    totalRot = regularizeAngle(totalRot + PI) - PI; // regularize to [-PI,PI]

    double angleDefect = geometry->angleDefect(v);
    cout << endl;
    // cout << "totalRot = " << totalRot << endl;
    // cout << "angleDefect = " << angleDefect << endl;
    // cout << totalRot - angleDefect << endl;
    // cout << totalRot + angleDefect << endl;
    // cout << totalRot - angleDefect / nSym << endl;
    // cout << totalRot + angleDefect / nSym << endl;
    // cout << totalRot - angleDefect * nSym << endl;
    cout << totalRot + angleDefect * nSym << endl;
    totalRot += angleDefect * nSym;

    // Compute the net rotation and corresponding index
    // if (interesting) {
    cout << "totalRot/2PI = " << (totalRot / (2 * PI)) << endl;
    // }
    int index = static_cast<int>(std::round(totalRot / (2 * PI))); // should be very close to a multiple of 2PI
    indices[v] = index;
  }

  return indices;
}

} // namespace geometrycentral
