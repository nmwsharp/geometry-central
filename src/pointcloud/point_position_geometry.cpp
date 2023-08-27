#include "geometrycentral/pointcloud/point_position_geometry.h"

#include "geometrycentral/pointcloud/local_triangulation.h"
#include "geometrycentral/surface/intrinsic_mollification.h"
#include "geometrycentral/surface/simple_idt.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/tufted_laplacian.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

namespace geometrycentral {
namespace pointcloud {


// clang-format off
PointPositionGeometry::PointPositionGeometry(PointCloud& cloud_, const PointData<Vector3>& positions_)
    : cloud(cloud_), positions(positions_),
      
  // Construct the dependency graph of managed quantities and their callbacks

  pointIndicesQ             (&pointIndices,             std::bind(&PointPositionGeometry::computePointIndices, this),          quantities),
  neighborsQ                (&neighbors,                std::bind(&PointPositionGeometry::computeNeighbors, this),             quantities),
  normalsQ                  (&normals,                  std::bind(&PointPositionGeometry::computeNormals, this),               quantities),
  tangentBasisQ             (&tangentBasis,             std::bind(&PointPositionGeometry::computeTangentBasis, this),          quantities),
  tangentCoordinatesQ       (&tangentCoordinates,       std::bind(&PointPositionGeometry::computeTangentCoordinates, this),             quantities),
  tangentTransportQ         (&tangentTransport,         std::bind(&PointPositionGeometry::computeTangentTransport, this),             quantities),

  tuftedTriPair{&tuftedMesh, &tuftedGeom},
  tuftedTriangulationQ      (&tuftedTriPair,            std::bind(&PointPositionGeometry::computeTuftedTriangulation, this),   quantities),

  // operators
  laplacianQ                (&laplacian,                std::bind(&PointPositionGeometry::computeLaplacian, this),             quantities),
  connectionLaplacianQ      (&connectionLaplacian,      std::bind(&PointPositionGeometry::computeConnectionLaplacian, this),   quantities),
  gradientQ                 (&gradient,                 std::bind(&PointPositionGeometry::computeGradient, this),              quantities)

  {
  }
// clang-format on

PointPositionGeometry::PointPositionGeometry(PointCloud& cloud_)
    : PointPositionGeometry(cloud_, PointData<Vector3>(cloud_)) {}

PointPositionGeometry::~PointPositionGeometry() {}

void PointPositionGeometry::refreshQuantities() {
  for (DependentQuantity* q : quantities) {
    q->computed = false;
  }
  for (DependentQuantity* q : quantities) {
    q->ensureHaveIfRequired();
  }
}

void PointPositionGeometry::purgeQuantities() {
  for (DependentQuantity* q : quantities) {
    q->clearIfNotRequired();
  }
}


// Point indices
void PointPositionGeometry::computePointIndices() { pointIndices = cloud.getPointIndices(); }
void PointPositionGeometry::requirePointIndices() { pointIndicesQ.require(); }
void PointPositionGeometry::unrequirePointIndices() { pointIndicesQ.unrequire(); }

// Neighbors
void PointPositionGeometry::computeNeighbors() { neighbors.reset(new Neighborhoods(cloud, positions, kNeighborSize)); }
void PointPositionGeometry::requireNeighbors() { neighborsQ.require(); }
void PointPositionGeometry::unrequireNeighbors() { neighborsQ.unrequire(); }

// Normals
void PointPositionGeometry::computeNormals() {

  neighborsQ.ensureHave();

  normals = PointData<Vector3>(cloud);

  for (Point p : cloud.points()) {
    size_t nNeigh = neighbors->neighbors[p].size();
    Vector3 center = positions[p];
    Eigen::MatrixXd localMat(3, nNeigh);

    for (size_t iN = 0; iN < nNeigh; iN++) {
      Vector3 neighPos = positions[neighbors->neighbors[p][iN]] - center;
      localMat(0, iN) = neighPos.x;
      localMat(1, iN) = neighPos.y;
      localMat(2, iN) = neighPos.z;
    }

    // Smallest singular vector is best normal
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(localMat, Eigen::ComputeThinU);
    Eigen::Vector3d bestNormal = svd.matrixU().col(2);

    Vector3 N{bestNormal(0), bestNormal(1), bestNormal(2)};
    N = unit(N);
    normals[p] = N;
  }
}
void PointPositionGeometry::requireNormals() { normalsQ.require(); }
void PointPositionGeometry::unrequireNormals() { normalsQ.unrequire(); }

// Tangent basis
void PointPositionGeometry::computeTangentBasis() {
  normalsQ.ensureHave();

  tangentBasis = PointData<std::array<Vector3, 2>>(cloud);
  for (Point p : cloud.points()) {
    tangentBasis[p] = normals[p].buildTangentBasis();
  }
}
void PointPositionGeometry::requireTangentBasis() { tangentBasisQ.require(); }
void PointPositionGeometry::unrequireTangentBasis() { tangentBasisQ.unrequire(); }

// Tangent coordintes
void PointPositionGeometry::computeTangentCoordinates() {
  neighborsQ.ensureHave();
  tangentBasisQ.ensureHave();
  normalsQ.ensureHave();

  tangentCoordinates = PointData<std::vector<Vector2>>(cloud);
  for (Point p : cloud.points()) {
    size_t nNeigh = neighbors->neighbors[p].size();
    tangentCoordinates[p].resize(nNeigh);
    Vector3 center = positions[p];
    Vector3 normal = normals[p];
    Vector3 basisX = tangentBasis[p][0];
    Vector3 basisY = tangentBasis[p][1];

    for (size_t iN = 0; iN < nNeigh; iN++) {
      Vector3 vec = positions[neighbors->neighbors[p][iN]] - center;
      vec = vec.removeComponent(normal);

      Vector2 coord{dot(basisX, vec), dot(basisY, vec)};
      tangentCoordinates[p][iN] = coord;
    }
  }
}
void PointPositionGeometry::requireTangentCoordinates() { tangentCoordinatesQ.require(); }
void PointPositionGeometry::unrequireTangentCoordinates() { tangentCoordinatesQ.unrequire(); }

// Tangent transport
void PointPositionGeometry::computeTangentTransport() {
  neighborsQ.ensureHave();
  normalsQ.ensureHave();
  tangentBasisQ.ensureHave();

  tangentTransport = PointData<std::vector<Vector2>>(cloud);
  for (Point p : cloud.points()) {
    size_t nNeigh = neighbors->neighbors[p].size();
    tangentTransport[p].resize(nNeigh);
    for (size_t iN = 0; iN < nNeigh; iN++) {
      Point pN = neighbors->neighbors[p][iN];
      tangentTransport[p][iN] = transportBetween(p, pN);
    }
  }
}
void PointPositionGeometry::requireTangentTransport() { tangentTransportQ.require(); }
void PointPositionGeometry::unrequireTangentTransport() { tangentTransportQ.unrequire(); }

// Tufted triangulation
void PointPositionGeometry::computeTuftedTriangulation() {
  neighborsQ.ensureHave();
  tangentCoordinatesQ.ensureHave();

  using namespace surface;

  PointData<std::vector<std::array<Point, 3>>> localTriPoint = buildLocalTriangulations(cloud, *this, true);

  // == Make a mesh
  std::vector<std::vector<size_t>> allTris = handleToFlatInds(cloud, localTriPoint);
  std::vector<Vector3> posRaw(cloud.nPoints());
  for (size_t iP = 0; iP < posRaw.size(); iP++) {
    posRaw[iP] = positions[iP];
  }

  // Make a mesh, read off its
  std::unique_ptr<VertexPositionGeometry> posGeom;
  std::tie(tuftedMesh, posGeom) = makeSurfaceMeshAndGeometry(allTris, posRaw);
  posGeom->requireEdgeLengths();
  EdgeData<double> tuftedEdgeLengths = posGeom->edgeLengths;

  // Mollify
  mollifyIntrinsic(*tuftedMesh, tuftedEdgeLengths, 1e-5);

  // Build the cover
  buildIntrinsicTuftedCover(*tuftedMesh, tuftedEdgeLengths);

  flipToDelaunay(*tuftedMesh, tuftedEdgeLengths);

  // Create the geometry object
  tuftedGeom.reset(new EdgeLengthGeometry(*tuftedMesh, tuftedEdgeLengths));
}
void PointPositionGeometry::requireTuftedTriangulation() { tuftedTriangulationQ.require(); }
void PointPositionGeometry::unrequireTuftedTriangulation() { tuftedTriangulationQ.unrequire(); }

// === Operators / Matrices


// Laplacian
void PointPositionGeometry::computeLaplacian() {
  tuftedTriangulationQ.ensureHave();

  tuftedGeom->requireCotanLaplacian();
  laplacian = tuftedGeom->cotanLaplacian;

  tuftedGeom->unrequireCotanLaplacian();
  tuftedGeom->purgeQuantities(); // overkill?
}
void PointPositionGeometry::requireLaplacian() { laplacianQ.require(); }
void PointPositionGeometry::unrequireLaplacian() { laplacianQ.unrequire(); }


// Connection Laplacian
void PointPositionGeometry::computeConnectionLaplacian() {
  laplacianQ.ensureHave();
  normalsQ.ensureHave();

  // Build the conection Laplacian by hitting each entry of the Laplacian with a tangent rotation r_ij between the
  // corresponding vertices

  std::vector<Eigen::Triplet<double>> triplets;

  // Helper to expand complex values in to 2x2 real entries as we build the matrix
  auto addComplexCoef = [&](size_t i, size_t j, Vector2 val, bool conj) {
    double c = conj ? -1. : 1.;

    triplets.emplace_back(2 * i + 0, 2 * j + 0, val.x);
    triplets.emplace_back(2 * i + 0, 2 * j + 1, -val.y * c);
    triplets.emplace_back(2 * i + 1, 2 * j + 0, val.y);
    triplets.emplace_back(2 * i + 1, 2 * j + 1, val.x * c);
  };

  for (int k = 0; k < laplacian.outerSize(); ++k) {
    for (typename SparseMatrix<double>::InnerIterator it(laplacian, k); it; ++it) {

      double thisVal = it.value();
      size_t i = it.row();
      size_t j = it.col();

      if (i == j) continue;

      Point pI = cloud.point(i);
      Point pJ = cloud.point(j);

      Vector2 r_ij;
      bool inverted;
      std::tie(r_ij, inverted) = transportBetweenOriented(pJ, pI);

      // bool inverted = dot(normals[pI], normals[pJ]) < 0.;

      addComplexCoef(i, j, thisVal * r_ij, inverted);
      addComplexCoef(i, i, -thisVal * Vector2{1., 0}, false);
      // triplets.emplace_back(i, j, thisVal * r_ij);
      // triplets.emplace_back(i, i, -thisVal);
    }
  }

  // Build the matrix
  connectionLaplacian = SparseMatrix<double>(2 * cloud.nPoints(), 2 * cloud.nPoints());
  connectionLaplacian.setFromTriplets(triplets.begin(), triplets.end());
}
void PointPositionGeometry::requireConnectionLaplacian() { connectionLaplacianQ.require(); }
void PointPositionGeometry::unrequireConnectionLaplacian() { connectionLaplacianQ.unrequire(); }


// Gradient
void PointPositionGeometry::computeGradient() {
  pointIndicesQ.ensureHave();
  neighborsQ.ensureHave();
  tangentCoordinatesQ.ensureHave();

  using namespace Eigen;

  std::vector<Eigen::Triplet<std::complex<double>>> triplets;

  for (Point p : cloud.points()) {
    size_t iPt = pointIndices[p];

    std::vector<Point>& neigh = neighbors->neighbors[p];

    // Build a local linear system describing the gradient
    size_t nNeigh = neigh.size();
    MatrixXd lhsMat = MatrixXd::Zero(nNeigh, 2);
    MatrixXd rhsMat = MatrixXd::Zero(nNeigh, nNeigh + 1);

    for (size_t iN = 0; iN < nNeigh; iN++) {
      rhsMat(iN, 0) = -1.;
      rhsMat(iN, iN + 1) = 1.;

      lhsMat(iN, 0) = tangentCoordinates[p][iN].x;
      lhsMat(iN, 1) = tangentCoordinates[p][iN].y;
    }

    // Solve
    MatrixXd gradMat = lhsMat.colPivHouseholderQr().solve(rhsMat);

    // Copy coefficients in to global mat
    triplets.emplace_back(iPt, iPt, std::complex<double>{gradMat(0, 0), gradMat(1, 0)});
    for (size_t iN = 0; iN < nNeigh; iN++) {
      Point n = neigh[iN];
      size_t nInd = pointIndices[n];
      triplets.emplace_back(iPt, nInd, std::complex<double>{gradMat(0, iN + 1), gradMat(1, iN + 1)});
    }
  }

  // Build the matrix
  gradient = SparseMatrix<std::complex<double>>(cloud.nPoints(), cloud.nPoints());
  gradient.setFromTriplets(triplets.begin(), triplets.end());
}
void PointPositionGeometry::requireGradient() { gradientQ.require(); }
void PointPositionGeometry::unrequireGradient() { gradientQ.unrequire(); }


// === Helpers

Vector2 PointPositionGeometry::transportBetween(Point pSource, Point pTarget) {

  Vector3 sourceN = normals[pSource];
  Vector3 sourceBasisX = tangentBasis[pSource][0];
  // Vector3 sourceBasisY = tangentBasis[pSource][1];
  Vector3 targetN = normals[pTarget];
  Vector3 targetBasisX = tangentBasis[pTarget][0];
  Vector3 targetBasisY = tangentBasis[pTarget][1];

  // Find rotation that aligns normals
  Vector3 axis = cross(targetN, sourceN);
  if (norm(axis) > 1e-6) {
    axis = unit(axis);
  } else {
    axis = sourceBasisX;
  }
  double angle = angleInPlane(sourceN, targetN, axis);

  // Rotate the axes in to the plane of this vertex
  Vector3 sourceXInTarget3 = sourceBasisX.rotateAround(axis, angle);
  Vector2 sourceXInTarget{dot(sourceXInTarget3, targetBasisX), dot(sourceXInTarget3, targetBasisY)};

  return sourceXInTarget;
}

std::tuple<Vector2, bool> PointPositionGeometry::transportBetweenOriented(Point pSource, Point pTarget) {

  Vector3 sourceN = normals[pSource];
  Vector3 sourceBasisX = tangentBasis[pSource][0];
  // Vector3 sourceBasisY = tangentBasis[pSource][1];
  Vector3 targetN = normals[pTarget];
  Vector3 targetBasisX = tangentBasis[pTarget][0];
  Vector3 targetBasisY = tangentBasis[pTarget][1];

  // Flip orientation
  bool inverted = false;
  if (dot(sourceN, targetN) < 0.) {
    targetN *= -1;
    targetBasisY *= -1;
    inverted = true;
  }

  // Find rotation that aligns normals
  Vector3 axis = cross(targetN, sourceN);
  if (norm(axis) > 1e-6) {
    axis = unit(axis);
  } else {
    axis = sourceBasisX;
  }
  double angle = angleInPlane(sourceN, targetN, axis);

  // Rotate the axes in to the plane of this vertex
  Vector3 sourceXInTarget3 = sourceBasisX.rotateAround(axis, angle);
  Vector2 sourceXInTarget{dot(sourceXInTarget3, targetBasisX), dot(sourceXInTarget3, targetBasisY)};

  return std::make_tuple(sourceXInTarget, inverted);
}

} // namespace pointcloud
} // namespace geometrycentral
