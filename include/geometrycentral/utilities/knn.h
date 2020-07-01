#pragma once

#include "geometrycentral/utilities/vector3.h"

#include <cstddef>
#include <memory>
#include <vector>

namespace geometrycentral {


class NearestNeighborFinder {
public:
  NearestNeighborFinder(const std::vector<Vector3>& points);
  ~NearestNeighborFinder();

  // Return the indices of points in the input set
  std::vector<size_t> kNearest(Vector3 query, size_t k);

  // Source index refers to a point in the input set, which will not appear in the output
  std::vector<size_t> kNearestNeighbors(size_t sourceInd, size_t k);

  // Return all neighbors within ball. `rad` should be the actual distance, not squared distance; we square internally
  // to pass to nanoflann
  std::vector<size_t> radiusSearch(Vector3 query, double rad);

private:
  // "PImpl" idiom
  class KNNImpl;
  std::unique_ptr<KNNImpl> impl;
};
  

} // namespace geometrycentral
