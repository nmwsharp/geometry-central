#pragma once

// === Implementations for datatypes which hold data stored on the mesh

namespace geometrycentral {

// === Actual function implementations

template <typename E, typename T>
MeshData<E, T>::MeshData() {}

template <typename E, typename T>
MeshData<E, T>::MeshData(ParentMeshT& parentMesh) : mesh(&parentMesh) {
  data.resize(elementCapacity<E>(mesh));
  data.setConstant(defaultValue);

  registerWithMesh();
}

template <typename E, typename T>
MeshData<E, T>::MeshData(ParentMeshT& parentMesh, T initVal) : mesh(&parentMesh), defaultValue(initVal) {
  data.resize(elementCapacity<E>(mesh));
  data.setConstant(defaultValue);

  registerWithMesh();
}

template <typename E, typename T>
MeshData<E, T>::MeshData(ParentMeshT& parentMesh, const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector)
    : MeshData(parentMesh) {
  fromVector(vector);
}

template <typename E, typename T>
MeshData<E, T>::MeshData(ParentMeshT& parentMesh, const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector,
                         const MeshData<E, size_t>& indexer)
    : MeshData(parentMesh) {
  fromVector(vector, indexer);
}

template <typename E, typename T>
MeshData<E, T>::MeshData(const MeshData<E, T>& other)
    : mesh(other.mesh), defaultValue(other.defaultValue), data(other.data) {
  registerWithMesh();
}

template <typename E, typename T>
MeshData<E, T>::MeshData(MeshData<E, T>&& other) noexcept
    : mesh(other.mesh), defaultValue(other.defaultValue), data(std::move(other.data)) {
  registerWithMesh();
}


template <typename E, typename T>
MeshData<E, T>& MeshData<E, T>::operator=(const MeshData<E, T>& other) {
  deregisterWithMesh();
  mesh = other.mesh;
  defaultValue = other.defaultValue;
  data = other.data;
  registerWithMesh();

  return *this;
}

template <typename E, typename T>
MeshData<E, T>& MeshData<E, T>::operator=(MeshData<E, T>&& other) noexcept {
  deregisterWithMesh();
  mesh = other.mesh;
  defaultValue = other.defaultValue;
  data = std::move(other.data);
  registerWithMesh();

  return *this;
}

template <typename E, typename T>
MeshData<E, T>::~MeshData() {
  deregisterWithMesh();
}


template <typename E, typename T>
void MeshData<E, T>::registerWithMesh() {

  // Used during default initialization
  if (mesh == nullptr) return;

  // Callback function on expansion
  std::function<void(size_t)> expandFunc = [&](size_t newSize) {
    size_t oldSize = data.size();
    Eigen::Matrix<T, Eigen::Dynamic, 1> newData(newSize);
    for (size_t i = 0; i < oldSize; i++) {
      newData[i] = data[i];
    }
    for (size_t i = oldSize; i < newSize; i++) {
      newData[i] = defaultValue;
    }
    data = newData;
  };


  // Callback function on compression
  std::function<void(const std::vector<size_t>&)> permuteFunc = [this](const std::vector<size_t>& perm) {
    // inline applyPermutation() for Eigen vectors
    Eigen::Matrix<T, Eigen::Dynamic, 1> newData(perm.size());
    for (size_t i = 0; i < perm.size(); i++) {
      newData[i] = data[perm[i]];
    }
    data = newData;
  };


  // Callback function on mesh delete
  std::function<void()> deleteFunc = [this]() {
    // Ensures that we don't try to remove with iterators on deconstruct of this object
    mesh = nullptr;
  };

  expandCallbackIt = getExpandCallbackList<E>(mesh).insert(getExpandCallbackList<E>(mesh).begin(), expandFunc);
  permuteCallbackIt = getPermuteCallbackList<E>(mesh).insert(getPermuteCallbackList<E>(mesh).end(), permuteFunc);
  deleteCallbackIt = mesh->meshDeleteCallbackList.insert(mesh->meshDeleteCallbackList.end(), deleteFunc);
}

template <typename E, typename T>
void MeshData<E, T>::deregisterWithMesh() {

  // Used during destruction of default-initializated object, for instance
  if (mesh == nullptr) return;

  getExpandCallbackList<E>(mesh).erase(expandCallbackIt);
  getPermuteCallbackList<E>(mesh).erase(permuteCallbackIt);
  mesh->meshDeleteCallbackList.erase(deleteCallbackIt);
}

template <typename E, typename T>
void MeshData<E, T>::fill(T val) {
  for (E e : iterateElements<E>(mesh)) {
    (*this)[e] = val;
  }
}

template <typename E, typename T>
inline void MeshData<E, T>::clear() {
  deregisterWithMesh();
  mesh = nullptr;
  defaultValue = T();
  data = Eigen::Matrix<T, Eigen::Dynamic, 1>();
}

template <typename E, typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> MeshData<E, T>::toVector() const {
  Eigen::Matrix<T, Eigen::Dynamic, 1> result(nElements<E>(mesh));
  size_t i = 0;
  for (E e : iterateElements<E>(mesh)) {
    result(i) = (*this)[e];
    i++;
  }
  return result;
}

template <typename E, typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> MeshData<E, T>::toVector(const MeshData<E, size_t>& indexer) const {
  size_t outSize = 0;
  for (E e : iterateElements<E>(mesh)) {
    if (indexer[e] != std::numeric_limits<size_t>::max()) outSize++;
  }
  Eigen::Matrix<T, Eigen::Dynamic, 1> result(outSize);
  for (E e : iterateElements<E>(mesh)) {
    if (indexer[e] != std::numeric_limits<size_t>::max()) {
      result(indexer[e]) = (*this)[e];
    }
  }
  return result;
}

template <typename E, typename T>
void MeshData<E, T>::fromVector(const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector) {
  if ((size_t)vector.rows() != nElements<E>(mesh)) throw std::runtime_error("Vector size does not match mesh size.");
  size_t i = 0;
  for (E e : iterateElements<E>(mesh)) {
    (*this)[e] = vector(i);
    i++;
  }
}

template <typename E, typename T>
void MeshData<E, T>::fromVector(const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector, const MeshData<E, size_t>& indexer) {
  for (E e : iterateElements<E>(mesh)) {
    if (indexer[e] != std::numeric_limits<size_t>::max()) {
      (*this)[e] = vector(indexer[e]);
    }
  }
}

template <typename E, typename T>
inline T& MeshData<E, T>::operator[](E e) {
#ifndef NDEBUG
  // These checks are a bit much to do on every access, so disable in release mode.
  assert(mesh != nullptr && "MeshData is uninitialized.");
  assert(e.getMesh() == mesh && "Attempted to access MeshData with member from wrong mesh");
#endif
  size_t i = dataIndexOfElement(mesh, e);
  return data[i];
}

template <typename E, typename T>
inline const T& MeshData<E, T>::operator[](E e) const {
#ifndef NDEBUG
  // These checks are a bit much to do on every access, so disable in release mode.
  assert(mesh != nullptr && "MeshData is uninitialized.");
  assert(e.getMesh() == mesh && "Attempted to access MeshData with member from wrong mesh");
#endif
  size_t i = dataIndexOfElement(mesh, e);
  return data[i];
}

template <typename E, typename T>
inline T& MeshData<E, T>::operator[](size_t i) {
#ifndef NDEBUG
  assert(i < (size_t)data.size() && "Attempted to access MeshData with out of bounds index");
#endif
  return data[i];
}

template <typename E, typename T>
inline const T& MeshData<E, T>::operator[](size_t i) const {
#ifndef NDEBUG
  assert(i < (size_t)data.size() && "Attempted to access MeshData with out of bounds index");
#endif
  return data[i];
}

template <typename E, typename T>
inline size_t MeshData<E, T>::size() const {
  if (mesh == nullptr) return 0;
  return elementCapacity<E>(mesh);
}

template <typename E, typename T>
typename MeshData<E,T>::DATA_T& MeshData<E, T>::raw() {
  return data;
}

template <typename E, typename T>
const typename MeshData<E,T>::DATA_T& MeshData<E, T>::raw() const {
  return data;
}

template <typename E, typename T>
typename MeshData<E, T>::ParentMeshT* MeshData<E, T>::getMesh() const {
  return mesh;
}

template <typename E, typename T>
inline MeshData<E, T> MeshData<E, T>::reinterpretTo(ParentMeshT& targetMesh) const {
  GC_SAFETY_ASSERT(nElements<E>(mesh) == nElements<E>(&targetMesh),
                   "meshes must have same number of elements to reinterpret");
  MeshData<E, T> newData(targetMesh, defaultValue);
  newData.data = data;
  return newData;
}

template <typename E, typename T>
inline void MeshData<E, T>::setDefault(T newDefault) {

  // set the new default
  defaultValue = newDefault;

  // need to ensure that any allocated values beyond the end of the currently used array get populated with this new
  // default (they are currently filled with the old default)
  // TODO this is pretty sloppy, and probably indicates a flaw in the design of this class... The default should
  // probably just be const instead.
  // TODO could avoid this (and some things in fill()) with an element_fill<E>() template
  size_t maxInd = 0;
  for (E e : iterateElements<E>(mesh)) {
    maxInd = std::max(maxInd, e.getIndex());
  }
  for (size_t i = maxInd + 1; i < (size_t)data.size(); i++) {
    data[i] = defaultValue;
  }
}

template <typename E, typename T>
inline T MeshData<E, T>::getDefault() const {
  return defaultValue;
}

// === Arithmetic overloads ===

template <typename E, typename T, typename U>
void checkMeshCompatible(const MeshData<E, T>& lhs, const MeshData<E, U>& rhs) {
  GC_SAFETY_ASSERT(lhs.getMesh() != nullptr && rhs.getMesh() != nullptr, "arguments must both be initialized");
  GC_SAFETY_ASSERT(lhs.getMesh() == rhs.getMesh(), "arguments be defined on same mesh");
}

template <typename E, typename T>
void checkMeshValid(const MeshData<E, T>& val) {
  GC_SAFETY_ASSERT(val.getMesh() != nullptr, "argument must be initialized");
}

} // namespace geometrycentral
