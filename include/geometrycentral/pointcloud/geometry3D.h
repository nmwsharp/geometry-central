#pragma once


#include "geometrycentral/pointcloud/point_cloud.h"
#include "geometrycentral/utilities/dependent_quantity.h"
#include "geometrycentral/utilities/vector2.h"
#include "geometrycentral/utilities/vector3.h"

namespace geometrycentral {
namespace pointcloud {

class Geometry3D {
public:
  // Data
  Geometry3D(PointCloud& mesh);

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

  // Point indices
  PointData<size_t> pointIndices;
  void requirePointIndices();
  void unrequirePointIndices();


protected:
  // All of the quantities available (subclasses will also add quantities to this list)
  // Note that this is a vector of non-owning pointers; the quantities are generally value members in the class, so
  // there is no need to delete these.
  std::vector<DependentQuantity*> quantities;

  // === Implementation details for quantities

  DependentQuantityD<PointData<size_t>> pointIndicesQ;
  virtual void computePointIndices();
};


} // namespace pointcloud
} // namespace geometrycentral
