#pragma once

#include "geometrycentral/numerical/linear_algebra_types.h"
#include "geometrycentral/pointcloud/point_cloud_element_types.h"
#include "geometrycentral/utilities/mesh_data.h"
#include "geometrycentral/utilities/utilities.h"

#include <list>
#include <memory>
#include <vector>

// NOTE: ipp includes at bottom of file

namespace geometrycentral {
namespace pointcloud {

// Typedefs and forward declarations

template <typename T>
using PointData = MeshData<Point, T>;


// ==========================================================
// ===================    Point Cloud    ====================
// ==========================================================

class PointCloud {

public:
  // Initialize a point cloud with N points
  PointCloud(size_t nPts);
  virtual ~PointCloud();

  // Number of points
  size_t nPoints() const;

  // Methods for range-based for loops
  // Example: for(Point p : cloud.points()) { ... }
  PointSet points();

  // Methods for accessing elements by index
  // only valid when the cloud is compressed
  Point point(size_t index);

  // Methods for obtaining canonical indices for cloud elements
  // When the cloud is compressed, these will be equivalent to `vertex.getIndex()`, etc. However, even when the cloud is
  // not compressed they will still provide a dense enumeration. Of course in some situations, custom indices might
  // instead be needed, this is just a default dense enumeration.
  PointData<size_t> getPointIndices();

  // == Utility functions
  void printStatistics() const; // print info about element counts to std::cout

  std::unique_ptr<PointCloud> copy() const;

  // Compress the cloud
  bool isCompressed() const;
  void compress();

  // == Mutation routines


  // == Callbacks that will be invoked on mutation to keep containers/iterators/etc valid.

  // Expansion callbacks
  // Argument is the new size of the element list. Elements up to this index may now be used (but _might_ not be
  // in use immediately).
  std::list<std::function<void(size_t)>> pointExpandCallbackList;

  // Compression callbacks
  // Argument is a permutation to a apply, such that d_new[i] = d_old[p[i]]. THe length of the permutation is hte size
  // of the new index space. Any elements with p[i] == INVALID_IND are unused in the new index space.
  std::list<std::function<void(const std::vector<size_t>&)>> pointPermuteCallbackList;

  // Delete callbacks
  // (this unfortunately seems to be necessary; objects which have registered their callbacks above
  // need to know not to try to de-register them if the cloud has been deleted)
  std::list<std::function<void()>> meshDeleteCallbackList;

  // Check capacity. Needed when implementing expandable containers for mutable meshes to ensure the contain can
  // hold a sufficient number of elements before the next resize event.
  size_t nPointsCapacity() const;

  // == Debugging, etc

  // Performs a sanity checks on the structure; throws on fail
  void validateConnectivity();

protected:
  // Note: these data structures / conventions are essentially a stripped-down version of surface_mesh.h for the sake of
  // consistency.

  // = Core arrays which hold the connectivity (here, just point validity)
  std::vector<char> pointValid; // true if not dead

  // Auxilliary arrays which cache other useful information

  // Track element counts. These are the actual number of valid elements,
  // not the size of the buffer that holds them.
  size_t nPointsCount = 0;

  // == Track the capacity and fill size of our buffers.
  // These give the capacity of the currently allocated buffer.
  // Note that this is _not_ defined to be std::vector::capacity(), it's the largest size such that arr[i] is legal (aka
  // arr.size()).
  size_t nPointsCapacityCount = 0;

  // These give the number of filled elements in the currently allocated buffer. This will also be the maximal index of
  // any element (except the weirdness of boundary loop faces). As elements get marked dead, nPointsCount decreases
  // but nPointsFillCount does not, so it denotes the end of the region in the buffer where elements have been
  // stored.
  size_t nPointsFillCount = 0;

  // The cloud is _compressed_ if all of the index spaces are dense. E.g. if thare are |N| points, then the points
  // are densely indexed from 0 ... |N|-1. The cloud can become not-compressed as deletions mark elements with
  // tombstones--this is how we support constant time deletion. Call compress() to re-index and return to usual dense
  // indexing.
  bool isCompressedFlag = true;

  uint64_t modificationTick = 1; // Increment every time the cloud is mutated in any way. Used to track staleness.

  // Hide copy and move constructors, we don't wanna mess with that
  PointCloud(const PointCloud& other) = delete;
  PointCloud& operator=(const PointCloud& other) = delete;
  PointCloud(PointCloud&& other) = delete;
  PointCloud& operator=(PointCloud&& other) = delete;

  // TODO although the skeleton is here, inserting/deleting points has not actually been fleshed out and tested.

  // Used to resize the cloud. Expands and shifts vectors as necessary.
  Point getNewPoint();

  // Detect dead elements
  bool pointIsDead(size_t ind) const;
  void deleteElement(Point p);

  // Compression and other helpers
  void compressPoints();
  void copyInternalFields(PointCloud& target) const;

  // Elements need direct access in to members to traverse
  friend class Point;
  friend struct PointRangeF;
};

} // namespace pointcloud
} // namespace geometrycentral

// clang-format off
// preserve ordering
#include "geometrycentral/pointcloud/point_cloud_logic_templates.ipp"
#include "geometrycentral/pointcloud/point_cloud_element_types.ipp"
#include "geometrycentral/pointcloud/point_cloud.ipp"
// clang-format on
