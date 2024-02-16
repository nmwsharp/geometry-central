#include "geometrycentral/utilities/knn.h"

#include "nanoflann.hpp"

using std::vector;

namespace geometrycentral {

namespace { // anonymous helpers

// Helper class to implement interface for nanoflann

struct Vector3Adaptor {
  std::vector<Vector3> rawPoints;
  inline size_t kdtree_get_point_count() const { return rawPoints.size(); }
  inline double kdtree_get_pt(const size_t idx, int dim) const { return rawPoints[idx][dim]; }
  template <class BBOX>
  bool kdtree_get_bbox(BBOX& bb) const {
    return false;
  }
};
} // namespace


class NearestNeighborFinder::KNNImpl {
public:
  // == Types
  typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, Vector3Adaptor>, Vector3Adaptor, 3>
      KD_Tree_T;

  // == Constructors
  KNNImpl(const std::vector<Vector3>& points) : data{points}, tree(3, data) { tree.buildIndex(); }

  // == Members

  // Adapted buffer holding data
  Vector3Adaptor data;

  // construct a kd-tree index:
  KD_Tree_T tree;

  // == Methods

  std::vector<size_t> kNearest(Vector3 query, size_t k) {
    if (k > data.rawPoints.size()) throw std::runtime_error("k is greater than number of points");
    std::vector<size_t> outInds(k);
    std::vector<double> outDistSq(k);
    tree.knnSearch(&query[0], k, &outInds[0], &outDistSq[0]);
    return outInds;
  }

  std::vector<size_t> kNearestNeighbors(size_t sourceInd, size_t k) {
    if ((k + 1) > data.rawPoints.size()) throw std::runtime_error("k+1 is greater than number of points");

    std::vector<size_t> outInds(k + 1);
    std::vector<double> outDistSq(k + 1);
    tree.knnSearch(&data.rawPoints[sourceInd][0], k + 1, &outInds[0], &outDistSq[0]);

    // remove source from list
    bool found = false;
    for (size_t i = 0; i < outInds.size(); i++) {
      if (outInds[i] == sourceInd) {
        outInds.erase(outInds.begin() + i);
        // outDistSq.erase(outDistSq.begin() + i);
        found = true;
        break;
      }
    }

    // if the source didn't appear, just remove the last point
    if (!found) {
      outInds.pop_back();
      outDistSq.pop_back();
    }

    return outInds;
  }

  std::vector<size_t> radiusSearch(Vector3 query, double rad) {
    // nanoflann wants a SQUARED raidus
    double radSq = rad * rad;

    std::vector<std::pair<size_t, double>> outPairs;
    tree.radiusSearch(&query[0], radSq, outPairs, nanoflann::SearchParams());

    // copy in to an array off indices
    std::vector<size_t> outInds(outPairs.size());
    for (size_t i = 0; i < outInds.size(); i++) {
      outInds[i] = outPairs[i].first;
    }

    return outInds;
  }
};


NearestNeighborFinder::NearestNeighborFinder(const std::vector<Vector3>& points) { impl.reset(new KNNImpl(points)); }
NearestNeighborFinder::~NearestNeighborFinder() = default;

std::vector<size_t> NearestNeighborFinder::kNearest(Vector3 query, size_t k) { return impl->kNearest(query, k); }
std::vector<size_t> NearestNeighborFinder::kNearestNeighbors(size_t sourceInd, size_t k) {
  return impl->kNearestNeighbors(sourceInd, k);
}

std::vector<size_t> NearestNeighborFinder::radiusSearch(Vector3 query, double rad) {
  return impl->radiusSearch(query, rad);
}

} // namespace geometrycentral
