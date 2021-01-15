#pragma once

#include "geometrycentral/pointcloud/neighborhoods.h"
#include "geometrycentral/pointcloud/point_cloud.h"
#include "geometrycentral/surface/edge_length_geometry.h"
#include "geometrycentral/utilities/dependent_quantity.h"
#include "geometrycentral/utilities/vector2.h"
#include "geometrycentral/utilities/vector3.h"

namespace geometrycentral {
namespace pointcloud {

class Geometry3D {
public:
  // Data
  Geometry3D(PointCloud& mesh);
  virtual ~Geometry3D();

  // == Members
  PointCloud& cloud;


  // == Utility methods

  // Recompute all require'd quantities from input data. Call this after e.g. repositioning a vertex or mutating the
  // mesh
  void refreshQuantities();

  // Clear out any cached quantities which were previously computed but are not currently required.
  void purgeQuantities();

  // Construct a geometry object on another point cloud identical to this one
  std::unique_ptr<Geometry3D> reinterpretTo(PointCloud& targetCloud);

  // Hide copy and move constructors; users are more likely to use them accidentally than intentionally.
  // See the explicit copy() function in derived classes.
  Geometry3D(const Geometry3D& other) = delete;
  Geometry3D& operator=(const Geometry3D& other) = delete;
  Geometry3D(Geometry3D&& other) = delete;
  Geometry3D& operator=(Geometry3D&& other) = delete;


  // Essential data
  PointData<Vector3> positions;
  unsigned int kNeighborSize = 30; // do NOT change except bfore calling refreshQuantities();


  // === Quantities

  // TODO switch all the neighborhood stuff to a flat representation w/ containers rather than ested vectors

  // Point indices
  PointData<size_t> pointIndices;
  void requirePointIndices();
  void unrequirePointIndices();

  // Neighbors (fixed size for now, according to kNeighborSize)
  std::unique_ptr<Neighborhoods> neighbors;
  void requireNeighbors();
  void unrequireNeighbors();

  // Normals
  PointData<Vector3> normals;
  void requireNormals();
  void unrequireNormals();

  // Tangent basis
  PointData<std::array<Vector3, 2>> tangentBasis;
  void requireTangentBasis();
  void unrequireTangentBasis();

  // Neighborhood tangent coordinates
  PointData<std::vector<Vector2>> tangentCoordinates;
  void requireTangentCoordinates();
  void unrequireTangentCoordinates();
  
  // Rotations to align tangent space
  // tangentTransport[i][j] holds the rotation which maps a vector in the tangent space of i to that of j.
  PointData<std::vector<Vector2>> tangentTransport;
  void requireTangentTransport();
  void unrequireTangentTransport();

  // Tufted IDT triangulation / geometry. The vertices are in correspondence with the (compressed) point cloud.
  std::unique_ptr<surface::SurfaceMesh> tuftedMesh;
  std::unique_ptr<surface::EdgeLengthGeometry> tuftedGeom;
  void requireTuftedTriangulation();
  void unrequireTuftedTriangulation();


  // === Operators / Matrices

  // Laplacian
  // (uses Tufted Intrinsic Laplacian: Sharp & Crane 2020 @ SGP)
  Eigen::SparseMatrix<double> laplacian;
  void requireLaplacian();
  void unrequireLaplacian();

  // Connection Laplacian
  // uses Tufted Intrinsic Laplacian as above
  Eigen::SparseMatrix<std::complex<double>> connectionLaplacian;
  void requireConnectionLaplacian();
  void unrequireConnectionLaplacian();

  // Gradient
  Eigen::SparseMatrix<std::complex<double>> gradient;
  void requireGradient();
  void unrequireGradient();

protected:
  // All of the quantities available (subclasses will also add quantities to this list)
  // Note that this is a vector of non-owning pointers; the quantities are generally value members in the class, so
  // there is no need to delete these.
  std::vector<DependentQuantity*> quantities;

  // === Implementation details for quantities

  DependentQuantityD<PointData<size_t>> pointIndicesQ;
  virtual void computePointIndices();

  DependentQuantityD<std::unique_ptr<Neighborhoods>> neighborsQ;
  virtual void computeNeighbors();

  DependentQuantityD<PointData<Vector3>> normalsQ;
  virtual void computeNormals();

  DependentQuantityD<PointData<std::array<Vector3, 2>>> tangentBasisQ;
  virtual void computeTangentBasis();

  DependentQuantityD<PointData<std::vector<Vector2>>> tangentCoordinatesQ;
  virtual void computeTangentCoordinates();
  
  DependentQuantityD<PointData<std::vector<Vector2>>> tangentTransportQ;
  virtual void computeTangentTransport();

  std::pair<std::unique_ptr<surface::SurfaceMesh>*, std::unique_ptr<surface::EdgeLengthGeometry>*> tuftedTriPair;
  DependentQuantityD<std::pair<std::unique_ptr<surface::SurfaceMesh>*, std::unique_ptr<surface::EdgeLengthGeometry>*>>
      tuftedTriangulationQ;
  virtual void computeTuftedTriangulation();

  // === Operators / Matrices

  // Laplacian
  DependentQuantityD<Eigen::SparseMatrix<double>> laplacianQ;
  virtual void computeLaplacian();

  // Connection Laplacian
  DependentQuantityD<Eigen::SparseMatrix<std::complex<double>>> connectionLaplacianQ;
  virtual void computeConnectionLaplacian();

  // Gradient
  DependentQuantityD<Eigen::SparseMatrix<std::complex<double>>> gradientQ;
  virtual void computeGradient();
};


} // namespace pointcloud
} // namespace geometrycentral
