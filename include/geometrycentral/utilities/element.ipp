#pragma once

// Implementations for halfedge_mesh_types.ipp

namespace geometrycentral {

// clang-format off

// ==========================================================
// ================      Base Element      ==================
// ==========================================================

// Constructors
template<typename T, typename M> 
Element<T, M>::Element() {}
template<typename T, typename M> 
Element<T, M>::Element(ParentMeshT* mesh_, size_t ind_) : mesh(mesh_), ind(ind_) {}
//template<typename T, typename M> 
//Element<T, M>::Element(const DynamicElement<T>& e) : mesh(e.getMesh()), ind(e.getIndex()) {}

// Comparators
template<typename T, typename M> 
inline bool Element<T, M>::operator==(const Element<T, M>& other) const { return ind == other.ind; }
template<typename T, typename M> 
inline bool Element<T, M>::operator!=(const Element<T, M>& other) const { return !(*this == other); }
template<typename T, typename M> 
inline bool Element<T, M>::operator>(const Element<T, M>& other) const { return ind > other.ind; }
template<typename T, typename M> 
inline bool Element<T, M>::operator>=(const Element<T, M>& other) const { return ind >= other.ind; }
template<typename T, typename M> 
inline bool Element<T, M>::operator<(const Element<T, M>& other) const { return ind < other.ind; }
template<typename T, typename M> 
inline bool Element<T, M>::operator<=(const Element<T, M>& other) const { return ind <= other.ind; }

template <typename T, typename M>
size_t Element<T, M>::getIndex() const { return ind; }

template <typename T, typename M>
M* Element<T, M>::getMesh() const { return mesh; }
}

namespace std {
template <typename T, typename M>
inline ostream& operator<<(ostream& output, const geometrycentral::Element<T, M>& e) {
  output << geometrycentral::typeShortName<T>() << "_" << e.getIndex();
  return output;
}

template <typename T, typename M>
inline std::string to_string(const geometrycentral::Element<T, M>& e) {
  ostringstream output;
  output << e;
  return output.str();
}
}

namespace geometrycentral {

// Dynamic element
/*
template<typename S> 
DynamicElement<S>::DynamicElement() {}

template<typename S> 
DynamicElement<S>::DynamicElement(ParentMeshT* mesh_, size_t ind_) : S(mesh_, ind_) {
  registerWithMesh();
}

template<typename S> 
DynamicElement<S>::DynamicElement(const S& e) : S(e) {
  registerWithMesh();
}
  
template <typename S>
DynamicElement<S>::DynamicElement(const DynamicElement& other) : S(other.mesh, other.ind) {
  registerWithMesh();
}

template <typename S>
DynamicElement<S>::DynamicElement(DynamicElement&& other) : S(other.mesh, other.ind) {
  registerWithMesh();
}

template <typename S>
DynamicElement<S>& DynamicElement<S>::operator=(const DynamicElement<S>& other) {
  deregisterWithMesh();
  this->mesh = other.mesh;
  this->ind = other.ind;
  registerWithMesh();
  return *this;
}

template <typename S>
DynamicElement<S>& DynamicElement<S>::operator=(DynamicElement<S>&& other) noexcept {
  deregisterWithMesh();
  this->mesh = other.mesh;
  this->ind = other.ind;
  registerWithMesh();
  return *this;
}

template<typename S> 
DynamicElement<S>::~DynamicElement() {
  deregisterWithMesh();
}

template<typename S> 
void DynamicElement<S>::registerWithMesh() {
  
  // Callback function on permutation
  std::function<void(const std::vector<size_t>&)> permuteFunc = [this](const std::vector<size_t>& perm) {
      throw std::runtime_error("not implemented");
  };

  // Callback function on mesh delete
  std::function<void()> deleteFunc = [this]() {
    // Ensures that we don't try to remove with iterators on deconstruct of this object
    this->mesh = nullptr;
  };

  permuteCallbackIt = getPermuteCallbackList<S>(this->mesh).insert(getPermuteCallbackList<S>(this->mesh).end(), permuteFunc);
  deleteCallbackIt = this->mesh->meshDeleteCallbackList.insert(this->mesh->meshDeleteCallbackList.end(), deleteFunc);
}

template<typename S> 
void DynamicElement<S>::deregisterWithMesh() {
  if (this->mesh == nullptr) return;
  getPermuteCallbackList<S>(this->mesh).erase(permuteCallbackIt);
  this->mesh->meshDeleteCallbackList.erase(deleteCallbackIt);
}
*/

// clang-format on

} // namespace geometrycentral
