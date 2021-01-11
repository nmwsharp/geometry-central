#pragma once

// Implementations for halfedge_mesh_types.ipp

// Make the element types hashable (this _should_ be doable for just the parent class, but I couldn't sort out how)
namespace std {
template <>
struct hash<geometrycentral::pointcloud::Point> {
  std::size_t operator()(const geometrycentral::pointcloud::Point& e) const { return std::hash<size_t>{}(e.getIndex()); }
};
} // namespace std

namespace geometrycentral {
namespace pointcloud {

// ==========================================================
// ================        Point           ==================
// ==========================================================

// Constructors
// (see note in header, these should be inherited but aren't due to compiler issues)
inline Point::Point() {}
inline Point::Point(PointCloud* mesh_, size_t ind_) : Element(mesh_, ind_) {}
// inline Point::Point(const DynamicElement<Point>& e) : Element(e.getMesh(), e.getIndex()) {}

// Navigators
//inline Halfedge Point::halfedge() const { return Halfedge(mesh, mesh->vHalfedge(ind)); }

// Properties
//inline bool Point::isBoundary() const { return !halfedge().twin().isInterior(); }

// Navigation iterators
//inline NavigationSetBase<VertexIncomingHalfedgeNavigator> Point::incomingHalfedges() const {
  //return NavigationSetBase<VertexIncomingHalfedgeNavigator>(halfedge().prevOrbitFace());
//}

// == Range iterators
inline bool PointRangeF::elementOkay(const PointCloud& mesh, size_t ind) { return !mesh.pointIsDead(ind); }

// == Navigation iterators

} // namespace pointcloud
} // namespace geometrycentral
