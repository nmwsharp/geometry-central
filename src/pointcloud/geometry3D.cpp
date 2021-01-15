#include "geometrycentral/pointcloud/geometry3D.h"

#include "geometrycentral/pointcloud/local_triangulation.h"
#include "geometrycentral/surface/intrinsic_mollification.h"
#include "geometrycentral/surface/simple_idt.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/tufted_laplacian.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

namespace geometrycentral {
namespace pointcloud {


// clang-format off
Geometry3D::Geometry3D(PointCloud& cloud_)
    : cloud(cloud_), positions(cloud),
      
  // Construct the dependency graph of managed quantities and their callbacks

  pointIndicesQ             (&pointIndices,             std::bind(&Geometry3D::computePointIndices, this),          quantities),
  neighborsQ                (&neighbors,                std::bind(&Geometry3D::computeNeighbors, this),             quantities),
  normalsQ                  (&normals,                  std::bind(&Geometry3D::computeNormals, this),               quantities),
  tangentBasisQ             (&tangentBasis,             std::bind(&Geometry3D::computeTangentBasis, this),          quantities),
  tangentCoordinatesQ       (&tangentCoordinates,       std::bind(&Geometry3D::computeTangentCoordinates, this),             quantities),
  tangentTransportQ         (&tangentTransport,         std::bind(&Geometry3D::computeTangentTransport, this),             quantities),

  tuftedTriPair{&tuftedMesh, &tuftedGeom},
  tuftedTriangulationQ      (&tuftedTriPair,            std::bind(&Geometry3D::computeTuftedTriangulation, this),   quantities),

  // operators
  laplacianQ                (&laplacian,                std::bind(&Geometry3D::computeLaplacian, this),             quantities),
  connectionLaplacianQ      (&connectionLaplacian,      std::bind(&Geometry3D::computeConnectionLaplacian, this),   quantities),
  gradientQ                 (&gradient,                 std::bind(&Geometry3D::computeGradient, this),              quantities)

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

// Tangent transport
void Geometry3D::computeTangentTransport() {
  neighborsQ.ensureHave();
  tangentBasisQ.ensureHave();

  // TODO
}
void Geometry3D::requireTangentTransport() { tangentTransportQ.require(); }
void Geometry3D::unrequireTangentTransport() { tangentTransportQ.unrequire(); }

// Tufted triangulation
void Geometry3D::computeTuftedTriangulation() {
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

  // Flip to delaunay
  flipToDelaunay(*tuftedMesh, tuftedEdgeLengths);

  // Create the geometry object
  tuftedGeom.reset(new EdgeLengthGeometry(*tuftedMesh, tuftedEdgeLengths));
}
void Geometry3D::requireTuftedTriangulation() { tuftedTriangulationQ.require(); }
void Geometry3D::unrequireTuftedTriangulation() { tuftedTriangulationQ.unrequire(); }

// === Operators / Matrices


// Laplacian
void Geometry3D::computeLaplacian() {
  tuftedTriangulationQ.ensureHave();

  tuftedGeom->requireCotanLaplacian();
  laplacian = tuftedGeom->cotanLaplacian;

  tuftedGeom->unrequireCotanLaplacian();
  tuftedGeom->purgeQuantities(); // overkill?
}
void Geometry3D::requireLaplacian() { laplacianQ.require(); }
void Geometry3D::unrequireLaplacian() { laplacianQ.unrequire(); }


// Connection Laplacian
void Geometry3D::computeConnectionLaplacian() {
  laplacianQ.ensureHave();

}
void Geometry3D::requireConnectionLaplacian() { connectionLaplacianQ.require(); }
void Geometry3D::unrequireConnectionLaplacian() { connectionLaplacianQ.unrequire(); }


// Gradient
void Geometry3D::computeGradient() {}
void Geometry3D::requireGradient() { gradientQ.require(); }
void Geometry3D::unrequireGradient() { gradientQ.unrequire(); }

} // namespace pointcloud
} // namespace geometrycentral
