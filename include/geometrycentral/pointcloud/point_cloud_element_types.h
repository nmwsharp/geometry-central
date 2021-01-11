#pragma once

#include "geometrycentral/utilities/element.h"
#include "geometrycentral/utilities/element_iterators.h"
#include "geometrycentral/utilities/mesh_data.h"
#include "geometrycentral/utilities/utilities.h"

#include <array>
#include <cstddef>
#include <iostream>
#include <list>
#include <typeindex>
#include <unordered_set>

namespace geometrycentral {
namespace pointcloud {

// === Types and inline methods for the halfedge mesh pointer and datatypes
class PointCloud;

class Point;

// struct VertexAdjacentVertexNavigator; forward-eclare navigators

// ==========================================================
// ================        Vertex          ==================
// ==========================================================

class Point: public Element<Point, PointCloud> {
public:
  // Constructors
  // inheriting constructor would work here, and replace the constructors below, but gcc-5 erroneously rejects the combo
  // with CRTP :( perhaps resurrect here and in other elements below once gcc-5 is sufficiently old
  // using Element<Vertex>::Element;
  Point();                             // construct an empty (null) element
  Point(PointCloud* mesh, size_t ind); // construct pointing to the i'th element of that type on a mesh.
  // Vertex(const DynamicElement<Vertex>& e); // construct from a dynamic element of matching type

  // Navigators
  // bool isDead() const;

  // Properties
  // bool isBoundary() const;

  // Iterators
  // NavigationSetBase<VertexAdjacentVertexNavigator> adjacentVertices() const;
};

// using DynamicVertex = DynamicElement<Vertex>;

// == Range iterators

// All vertices
struct PointRangeF {
  static bool elementOkay(const PointCloud& mesh, size_t ind);
  typedef Point Etype;
  typedef PointCloud ParentMeshT;
};
typedef RangeSetBase<PointRangeF> PointSet;


// ==========================================================
// ===============   Navigation Iterators   =================
// ==========================================================

// == Point


} // namespace pointcloud

// Declare specializations of the logic templates. This is important, because these need to be declared before any of
// the templates using them are instantiated.


template <>
inline size_t nElements<pointcloud::Point>(pointcloud::PointCloud* mesh);

template <>
inline size_t elementCapacity<pointcloud::Point>(pointcloud::PointCloud* mesh);

template <>
inline size_t dataIndexOfElement<pointcloud::Point>(pointcloud::PointCloud* mesh, pointcloud::Point e);

template <>
struct ElementSetType<pointcloud::Point> {
  typedef pointcloud::PointSet type;
};

template <>
inline pointcloud::PointSet iterateElements<pointcloud::Point>(pointcloud::PointCloud* mesh);

template <>
inline std::list<std::function<void(size_t)>>& getExpandCallbackList<pointcloud::Point>(pointcloud::PointCloud* mesh);

template <>
inline std::list<std::function<void(const std::vector<size_t>&)>>&
getPermuteCallbackList<pointcloud::Point>(pointcloud::PointCloud* mesh);

template <>
inline std::string typeShortName<pointcloud::Point>();


} // namespace geometrycentral
