
// === Helpers which will allow us abstract over types
// The corresponding template declarations are given in halfedge_element_types.h

namespace geometrycentral {


template <>
inline size_t nElements<pointcloud::Point>(pointcloud::PointCloud* mesh) {
  return mesh->nPoints();
}

template <>
inline size_t elementCapacity<pointcloud::Point>(pointcloud::PointCloud* mesh) {
  return mesh->nPointsCapacity();
}

template <>
inline size_t dataIndexOfElement<pointcloud::Point>(pointcloud::PointCloud* mesh, pointcloud::Point e) {
  return e.getIndex();
}

template <>
inline pointcloud::PointSet iterateElements<pointcloud::Point>(pointcloud::PointCloud* mesh) {
  return mesh->points();
}

template <>
inline std::list<std::function<void(size_t)>>& getExpandCallbackList<pointcloud::Point>(pointcloud::PointCloud* mesh) {
  return mesh->pointExpandCallbackList;
}

template <>
inline std::list<std::function<void(const std::vector<size_t>&)>>&
getPermuteCallbackList<pointcloud::Point>(pointcloud::PointCloud* mesh) {
  return mesh->pointPermuteCallbackList;
}

template <>
inline std::string typeShortName<pointcloud::Point>() {
  return "p";
}

} // namespace geometrycentral
