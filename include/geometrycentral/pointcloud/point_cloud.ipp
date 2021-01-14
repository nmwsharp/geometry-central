#pragma once

#include "geometrycentral/numerical/linear_algebra_utilities.h"

namespace geometrycentral {
namespace pointcloud {


// Methods for getting number of mesh elements
inline size_t PointCloud::nPoints() const { return nPointsCount; }

// Capacities
inline size_t PointCloud::nPointsCapacity() const { return nPointsCapacityCount; }

// Detect dead elements
inline bool PointCloud::pointIsDead(size_t iP) const { return !pointValid[iP]; }

// Methods for iterating over mesh elements w/ range-based for loops ===========
inline PointSet PointCloud::points() { return PointSet(this, 0, nPointsFillCount); }

// Methods for accessing elements by index =====================================
// Note that these are only valid when the mesh is compressed.
inline Point PointCloud::point(size_t index) { return Point(this, index); }

// Misc utility methods =====================================
inline bool PointCloud::isCompressed() const { return isCompressedFlag; }


} // namespace pointcloud
} // namespace geometrycentral
