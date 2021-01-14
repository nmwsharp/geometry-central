#include "geometrycentral/pointcloud/geometry3D.h"

namespace geometrycentral {
namespace pointcloud {


// clang-format off
Geometry3D::Geometry3D(PointCloud& cloud_)
    : cloud(cloud_), positions(cloud),
      
  // Construct the dependency graph of managed quantities and their callbacks

  pointIndicesQ           (&pointIndices,         std::bind(&Geometry3D::computePointIndices, this),             quantities),
  neighborsQ              (&neighbors,            std::bind(&Geometry3D::computeNeighbors, this),             quantities),
  normalsQ                (&normals,              std::bind(&Geometry3D::computeNormals, this),             quantities),
  tangentBasisQ           (&tangentBasis,         std::bind(&Geometry3D::computeTangentBasis, this),             quantities),
  tangentCoordinatesQ     (&tangentCoordinates,   std::bind(&Geometry3D::computeTangentCoordinates, this),             quantities)

  {
  }
// clang-format on

Geometry3D::~Geometry3D() {}

void Geometry3D::refreshQuantities() {
  for (DependentQuantity* q : quantities) {
    q->computed = false;
  }
  for (DependentQuantity* q : quantities) {
    q->ensureHaveIfRequired();
  }
}

void Geometry3D::purgeQuantities() {
  for (DependentQuantity* q : quantities) {
    q->clearIfNotRequired();
  }
}


// Point indices
void Geometry3D::computePointIndices() { pointIndices = cloud.getPointIndices(); }
void Geometry3D::requirePointIndices() { pointIndicesQ.require(); }
void Geometry3D::unrequirePointIndices() { pointIndicesQ.unrequire(); }

// Neighbors
void Geometry3D::computeNeighbors() { neighbors.reset(new Neighborhoods(cloud, positions, kNeighborSize)); }
void Geometry3D::requireNeighbors() { neighborsQ.require(); }
void Geometry3D::unrequireNeighbors() { neighborsQ.unrequire(); }

// Normals
void Geometry3D::computeNormals() {
  using namespace Eigen;

  neighborsQ.ensureHave();

  normals = PointData<Vector3>(cloud);

  for (Point p : cloud.points()) {
    size_t nNeigh = neighbors->neighbors[p].size();
    Vector3 center = positions[p];
    MatrixXd localMat(3, nNeigh);

    for (size_t iN = 0; iN < nNeigh; iN++) {
      Vector3 neighPos = positions[neighbors->neighbors[p][iN]] - center;
      localMat(0, iN) = neighPos.x;
      localMat(1, iN) = neighPos.y;
      localMat(2, iN) = neighPos.z;
    }

    // Smallest singular vector is best normal
    JacobiSVD<MatrixXd> svd(localMat, ComputeThinU);
    Vector3d bestNormal = svd.matrixU().col(2);

    Vector3 N{bestNormal(0), bestNormal(1), bestNormal(2)};
    N = unit(N);
    normals[p] = N;
  }
}
void Geometry3D::requireNormals() { normalsQ.require(); }
void Geometry3D::unrequireNormals() { normalsQ.unrequire(); }

// Tangent basis
void Geometry3D::computeTangentBasis() {
  normalsQ.ensureHave();

  tangentBasis = PointData<std::array<Vector3, 2>>(cloud);
  for (Point p : cloud.points()) {
    tangentBasis[p] = normals[p].buildTangentBasis();
  }
}
void Geometry3D::requireTangentBasis() { tangentBasisQ.require(); }
void Geometry3D::unrequireTangentBasis() { tangentBasisQ.unrequire(); }

// Tangent coordintes
void Geometry3D::computeTangentCoordinates() {
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
void Geometry3D::requireTangentCoordinates() { tangentCoordinatesQ.require(); }
void Geometry3D::unrequireTangentCoordinates() { tangentCoordinatesQ.unrequire(); }


} // namespace pointcloud
} // namespace geometrycentral
