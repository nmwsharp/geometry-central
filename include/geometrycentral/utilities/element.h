#pragma once

#include "geometrycentral/utilities/utilities.h"


#include <functional>
#include <iostream>
#include <list>
#include <tuple>

namespace geometrycentral {

// === Templated helper functions
// Each element type must provide specializations for all of these functions

// clang-format off

// Current count of this element in the mesh
template <typename E> size_t nElements(typename E::ParentMeshT* mesh) { return INVALID_IND; }

// Capacity of element type in mesh (containers should be at least this big before next resize)
template <typename E> size_t elementCapacity(typename E::ParentMeshT* mesh) { return INVALID_IND; }

// Canonical index for this element (not always an enumeration)
template <typename E> size_t dataIndexOfElement(typename E::ParentMeshT* mesh, E e) { return INVALID_IND; }

// The set type used to iterate over all elements of this type
template <typename E> struct ElementSetType { typedef std::tuple<> type; }; // nonsense default value

template <typename E> typename ElementSetType<E>::type iterateElements(typename E::ParentMeshT* mesh) { return INVALID_IND; }

template <typename E> std::list<std::function<void(size_t)>>& getExpandCallbackList(typename E::ParentMeshT* mesh);

template <typename E> std::list<std::function<void(const std::vector<size_t>&)>>& getPermuteCallbackList(typename E::ParentMeshT* mesh);

template <typename E> std::string typeShortName() { return "X"; }
// clang-format on

// Forward declare
template <typename N>
class NavigationSetBase;

// ==========================================================
// ================      Base Element      ==================
// ==========================================================

// Forward-declare dynamic equivalent so we can declare a conversion constructor
// template <typename S>
// class DynamicElement;

// == Base type for shared logic between elements.
//
// This class uses the "curiously recurring template pattern" (CRTP) because it allows to implement more shared
// functionality in this template. Instantiations of the template will be templted on its child types, like `Vertex :
// public Element<Vertex>`. Because the parent class will know its child type at compile time, it can customize
// functionality based on that child type using the helpers above. This is essentially "compile time polymorphism",
// which allows us to share common functionality without paying a virtual function runtime cost.
template <typename T, typename M>
class Element {

public:
  using ParentMeshT = M;

  Element();                              // construct an empty (null) element
  Element(ParentMeshT* mesh, size_t ind); // construct pointing to the i'th element of that type on a mesh.
  // Element(const DynamicElement<T>& e);          // construct from a dynamic element of matching type

  inline bool operator==(const Element<T, M>& other) const;
  inline bool operator!=(const Element<T, M>& other) const;
  inline bool operator>(const Element<T, M>& other) const;
  inline bool operator>=(const Element<T, M>& other) const;
  inline bool operator<(const Element<T, M>& other) const;
  inline bool operator<=(const Element<T, M>& other) const;

  // Get the "index" associated with the element.
  // Note that these are not always a dense enumeration, and probably should not be accessed by "users" unless you are
  // monkeying around the halfedge mesh datastructure in some deep way. Generally prefer
  // `SurfaceMesh::getVertexIndices()` (etc) if you are looking for a set of indices for a linear algebra problem or
  // something.
  size_t getIndex() const;

  // Get the parent mesh on which the element is defined.
  ParentMeshT* getMesh() const;

  bool isDead() const;

protected:
  ParentMeshT* mesh = nullptr;
  size_t ind = INVALID_IND;

  // Friends
  friend struct std::hash<Element<T, M>>;
};

} // namespace geometrycentral

namespace std {
template <typename T, typename M>
std::ostream& operator<<(std::ostream& output, const geometrycentral::Element<T, M>& e);
template <typename T, typename M>
std::string to_string(const geometrycentral::Element<T, M>& e);
} // namespace std


namespace geometrycentral {

/*
// The equivalent dynamic pointers. These should be rarely used, but are guaranteed to be preserved through _all_ mesh
// operations, including compress().
template <typename S>
class DynamicElement : public S {
public:
  using ParentMeshT = typename S::ParentMeshT;

  DynamicElement();                                // construct an empty (null) element
  DynamicElement(ParentMeshT* mesh_, size_t ind_); // construct from an index as usual
  DynamicElement(const S& e);                      // construct from a non-dynamic element


  DynamicElement(const DynamicElement& other);
  DynamicElement(DynamicElement&& other);
  DynamicElement& operator=(const DynamicElement<S>& other);
  DynamicElement& operator=(DynamicElement<S>&& other) noexcept;

  ~DynamicElement();

  // returns a new plain old static element. Useful for chaining.
  S decay() const;

private:
  // References to the callbacks which keep the element valid. Keep these around to de-register on destruction.
  std::list<std::function<void(const std::vector<size_t>&)>>::iterator permuteCallbackIt;
  std::list<std::function<void()>>::iterator deleteCallbackIt;

  void registerWithMesh();
  void deregisterWithMesh();
};
*/


} // namespace geometrycentral

#include "geometrycentral/utilities/element.ipp"
