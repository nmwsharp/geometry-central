#include "geometrycentral/pointcloud/point_cloud.h"

#include "geometrycentral/utilities/combining_hash_functions.h"


namespace geometrycentral {
namespace pointcloud {

// ==========================================================
// ================      Construction      ==================
// ==========================================================

PointCloud::PointCloud(size_t nPts) {

  // Initialize counts
  nPointsCount = nPts;
  nPointsCapacityCount = nPts;
  nPointsFillCount = nPts;

  // All points are initially valid
  pointValid = std::vector<char>(nPts, true);

  isCompressedFlag = true;
}

PointCloud::~PointCloud() {}


// ==========================================================
// ================       Utilities        ==================
// ==========================================================

void PointCloud::printStatistics() const { std::cout << "Point cloud with # points =  " << nPoints() << std::endl; }


PointData<size_t> PointCloud::getPointIndices() {
  PointData<size_t> indices(*this);
  size_t i = 0;
  for (Point p : points()) {
    indices[p] = i;
    i++;
  }
  return indices;
}


std::unique_ptr<PointCloud> PointCloud::copy() const {
  PointCloud* newMesh = new PointCloud(0);
  copyInternalFields(*newMesh);
  return std::unique_ptr<PointCloud>(newMesh);
}


void PointCloud::copyInternalFields(PointCloud& target) const {
  // == Copy _all_ the fields!

  // Raw data buffers (underlying std::vectors duplicate storage automatically)
  target.pointValid = pointValid;

  // counts and flags
  target.nPointsCount = nPointsCount;
  target.nPointsCapacityCount = nPointsCapacityCount;
  target.nPointsFillCount = nPointsFillCount;

  target.isCompressedFlag = isCompressedFlag;

  // Note: _don't_ copy callbacks lists! New mesh has new callbacks
}


void PointCloud::validateConnectivity() {

  // Sanity check sizes and counts
  if (nPointsCount > nPointsFillCount) throw std::logic_error("point count > point fill");
  if (nPointsFillCount > nPointsCapacityCount) throw std::logic_error("point fill > point capacity");

  // Check for overflow / other unreasonable values
  if (nPointsCount > std::numeric_limits<uint64_t>::max() / 2) throw std::logic_error("point count overflow");
  if (nPointsFillCount > std::numeric_limits<uint64_t>::max() / 2) throw std::logic_error("point fill count overflow");
  if (nPointsCapacityCount > std::numeric_limits<uint64_t>::max() / 2)
    throw std::logic_error("point capacity count overflow");


  { // Count valid points
    size_t pCount = 0;
    for (Point p : points()) {
      pCount++;
    }
    if (pCount != nPoints()) throw std::logic_error("number of points does not match recount");
  }
}


Point PointCloud::getNewPoint() {

  // The boring case, when no resize is needed
  if (nPointsFillCount < nPointsCapacityCount) {
    // No work needed
  }
  // The intesting case, where vectors resize
  else {
    size_t newCapacity = nPointsCapacityCount * 2;

    // Resize internal arrays
    pointValid.resize(newCapacity);

    // Mark new points as dead
    for (size_t iP = nPointsCapacityCount; iP < pointValid.size(); iP++) {
      pointValid[iP] = false;
    }

    nPointsCapacityCount = newCapacity;

    // Invoke relevant callback functions
    for (auto& f : pointExpandCallbackList) {
      f(newCapacity);
    }
  }

  pointValid[nPointsFillCount] = true;
  nPointsFillCount++;
  nPointsCount++;

  modificationTick++;
  isCompressedFlag = false;
  return Point(this, nPointsFillCount - 1);
}

void PointCloud::deleteElement(Point v) {
  size_t iV = v.getIndex();

  pointValid[iV] = false;
  nPointsCount--;

  modificationTick++;
  isCompressedFlag = false;
}

void PointCloud::compressPoints() {
  // Build the compressing shift
  std::vector<size_t> newIndMap;                                // maps new ind -> old ind
  std::vector<size_t> oldIndMap(nPointsFillCount, INVALID_IND); // maps old ind -> new ind
  for (size_t i = 0; i < nPointsFillCount; i++) {
    if (!pointIsDead(i)) {
      oldIndMap[i] = newIndMap.size();
      newIndMap.push_back(i);
    }
  }


  // Permute & resize all per-point arrays
  pointValid = applyPermutation(pointValid, newIndMap);

  // Update counts
  nPointsFillCount = nPointsCount;
  nPointsCapacityCount = nPointsCount;

  // Invoke callbacks
  for (auto& f : pointPermuteCallbackList) {
    f(newIndMap);
  }
}


void PointCloud::compress() {
  if (isCompressed()) {
    return;
  }

  compressPoints();
  isCompressedFlag = true;
}


} // namespace pointcloud
} // namespace geometrycentral
