#include "geometrycentral/pointcloud/geometry3D.h"

namespace geometrycentral {
namespace pointcloud {


// clang-format off
Geometry3D::Geometry3D(PointCloud& cloud_)
    : cloud(cloud_), positions(cloud),
      
  // Construct the dependency graph of managed quantities and their callbacks

  pointIndicesQ           (&pointIndices,         std::bind(&Geometry3D::computePointIndices, this),          quantities)

  {
  }
// clang-format on

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


// Vertex indices
void Geometry3D::computePointIndices() { pointIndices = cloud.getPointIndices(); }
void Geometry3D::requirePointIndices() { pointIndicesQ.require(); }
void Geometry3D::unrequirePointIndices() { pointIndicesQ.unrequire(); }

} // namespace pointcloud
} // namespace geometrycentral
