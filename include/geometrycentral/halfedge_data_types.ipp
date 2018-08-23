#pragma once

// === Implementations for datatypes which hold data stored on the mesh

// Note: these implementations assuming that the default indexing convention for elements will be the same as the
// iteration order.

namespace geometrycentral {

// === Helpers which will allow us abstract over element access
// clang-format off

// Get number of elemenents
template <typename E> size_t nElements(HalfedgeMesh* mesh)            { return std::numeric_limits<size_t>::max(); }
template<> inline size_t nElements<VertexPtr  >(HalfedgeMesh* mesh)   { return mesh->nVertices();   }
template<> inline size_t nElements<FacePtr    >(HalfedgeMesh* mesh)   { return mesh->nFaces() + mesh->nBoundaryLoops(); }
template<> inline size_t nElements<EdgePtr    >(HalfedgeMesh* mesh)   { return mesh->nEdges();      }
template<> inline size_t nElements<HalfedgePtr>(HalfedgeMesh* mesh)   { return mesh->nHalfedges() + mesh->nImaginaryHalfedges();  }
template<> inline size_t nElements<CornerPtr  >(HalfedgeMesh* mesh)   { return mesh->nCorners();    }

// Index from element
template <typename E> size_t dataIndexOfElement(HalfedgeMesh* mesh, E e)            { return std::numeric_limits<size_t>::max(); }
template<> inline size_t dataIndexOfElement<VertexPtr  >(HalfedgeMesh* mesh, VertexPtr e)    { return e - mesh->vertex(0); }
template<> inline size_t dataIndexOfElement<FacePtr    >(HalfedgeMesh* mesh, FacePtr e)      { 
  if (e.isReal()) {
    return e - mesh->face(0);
  } else {
    return e - mesh->boundaryLoop(0);
  }
}
template<> inline size_t dataIndexOfElement<EdgePtr    >(HalfedgeMesh* mesh, EdgePtr e)      { return e - mesh->edge(0);      }
template<> inline size_t dataIndexOfElement<HalfedgePtr>(HalfedgeMesh* mesh, HalfedgePtr e)  { return e - mesh->halfedge(0);  }
template<> inline size_t dataIndexOfElement<CornerPtr  >(HalfedgeMesh* mesh, CornerPtr e)    { return e - mesh->corner(0);    }

// Set iterator type
template <class E> struct ElementSetType            { typedef E                 type; };
template <> struct ElementSetType<VertexPtr     >   { typedef VertexPtrSet      type; };
template <> struct ElementSetType<FacePtr       >   { typedef FacePtrSet        type; };
template <> struct ElementSetType<EdgePtr       >   { typedef EdgePtrSet        type; };
template <> struct ElementSetType<HalfedgePtr   >   { typedef HalfedgePtrSet    type; };
template <> struct ElementSetType<CornerPtr     >   { typedef CornerPtrSet      type; };

// Iterate through elements
template <typename E> typename ElementSetType<E>::type iterateElements(HalfedgeMesh* mesh)  { return std::numeric_limits<size_t>::max(); }
template<> inline VertexPtrSet      iterateElements<VertexPtr  >(HalfedgeMesh* mesh)   { return mesh->vertices();      }
template<> inline FacePtrSet        iterateElements<FacePtr    >(HalfedgeMesh* mesh)   { return mesh->faces();         }
template<> inline EdgePtrSet        iterateElements<EdgePtr    >(HalfedgeMesh* mesh)   { return mesh->edges();         }
template<> inline HalfedgePtrSet    iterateElements<HalfedgePtr>(HalfedgeMesh* mesh)   { return mesh->allHalfedges();  }
template<> inline CornerPtrSet      iterateElements<CornerPtr  >(HalfedgeMesh* mesh)   { return mesh->corners();       }


// clang-format on


// === Actual function implementations

template <typename E, typename T>
MeshData<E, T>::MeshData(HalfedgeMesh* parentMesh) : mesh(parentMesh) {
  if (parentMesh != nullptr) {
    data.resize(nElements<E>(parentMesh));
    fill(defaultValue);
  }
}

template <typename E, typename T>
MeshData<E, T>::MeshData(HalfedgeMesh* parentMesh, T initVal) : mesh(parentMesh), defaultValue(initVal) {
  if (parentMesh != nullptr) {
    data.resize(nElements<E>(parentMesh));
    fill(defaultValue);
  }
}

template <typename E, typename T>
MeshData<E, T>::MeshData(HalfedgeMesh* parentMesh, const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector)
    : MeshData(parentMesh) {
  fromVector(vector);
}

template <typename E, typename T>
MeshData<E, T>::MeshData(HalfedgeMesh* parentMesh, const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector,
                         const MeshData<E, size_t>& indexer)
    : MeshData(parentMesh) {
  fromVector(vector, indexer);
}

template <typename E, typename T>
void MeshData<E, T>::fill(T val) {
  std::fill(data.begin(), data.end(), val);
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
  assert(mesh != nullptr && "MeshData is uninitialized.");
  assert(e->parentMesh == mesh && "Attempted to access MeshData with member from wrong mesh");
#endif
  size_t i = dataIndexOfElement(mesh, e);
  return data[i];
}

template <typename E, typename T>
inline const T& MeshData<E, T>::operator[](E e) const {
#ifndef NDEBUG
  assert(mesh != nullptr && "MeshData is uninitialized.");
  assert(e->parentMesh == mesh && "Attempted to access MeshData with member from wrong mesh");
#endif
  size_t i = dataIndexOfElement(mesh, e);
  return data[i];
}

template <typename E, typename T>
inline T& MeshData<E, T>::operator[](size_t i) {
#ifndef NDEBUG
  assert(i < size() && "Attempted to access MeshData with out of bounds index");
#endif
  return data[i];
}

template <typename E, typename T>
inline const T& MeshData<E, T>::operator[](size_t i) const {
#ifndef NDEBUG
  assert(i < size() && "Attempted to access MeshData with out of bounds index");
#endif
  return data[i];
}

template <typename E, typename T>
inline size_t MeshData<E, T>::size() const {
  return nElements<E>(mesh);
}


} // namespace geometrycentral
