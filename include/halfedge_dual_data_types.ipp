#pragma once

#include "halfedge_data_macros.h"

namespace geometrycentral {

// Data on vertices
template <typename T>
DualVertexData<T>::DualVertexData(HalfedgeDual* parentMesh)
    : dualMesh(parentMesh) {
  if (parentMesh != nullptr) {
    data.resize(parentMesh->nVertices());
  }
}

template <typename T>
DualVertexData<T>::DualVertexData(HalfedgeMesh* parentMesh, T initVal)
    : DualVertexData(parentMesh) {
  fill(initVal);
}

template <typename T>
void DualVertexData<T>::fill(T val) {
  std::fill(data.begin(), data.end(), val);
}

template <typename T>
inline T& DualVertexData<T>::operator[](DualVertexPtr v) {
#ifndef NDEBUG
  assert(v->parentMesh == dualMesh &&
         "Attempted to access data with member from wrong dualMesh");
#endif
  unsigned int i = v - dualMesh->vertex(0);
  return data[i];
}

template <typename T>
inline const T& DualVertexData<T>::operator[](DualVertexPtr v) const {
#ifndef NDEBUG
  assert(v->parentMesh == dualMesh &&
         "Attempted access data with member from wrong dualMesh");
#endif
  unsigned int i = v - dualMesh->vertex(0);
  return data[i];
}

GC_INTERNAL_GENERATE_DATATYPE_OPERATOR_DEFINITIONS(DualVertexData, dualMesh)

// Data on edges
template <typename T>
DualEdgeData<T>::DualEdgeData(HalfedgeDual* parentMesh) : dualMesh(parentMesh) {
  if (parentMesh != nullptr) {
    data.resize(parentMesh->nEdges());
  }
}

template <typename T>
DualEdgeData<T>::DualEdgeData(HalfedgeMesh* parentMesh, T initVal)
    : DualEdgeData(parentMesh) {
  fill(initVal);
}

template <typename T>
void DualEdgeData<T>::fill(T val) {
  std::fill(data.begin(), data.end(), val);
}

template <typename T>
inline T& DualEdgeData<T>::operator[](DualEdgePtr e) {
#ifndef NDEBUG
  assert(e->parentMesh == dualMesh &&
         "Attempted access data with member from wrong dualMesh");
#endif
  unsigned int i = e - dualMesh->edge(0);
  return data[i];
}

template <typename T>
inline const T& DualEdgeData<T>::operator[](DualEdgePtr e) const {
#ifndef NDEBUG
  assert(e->parentMesh == dualMesh &&
         "Attempted access data with member from wrong dualMesh");
#endif
  unsigned int i = e - dualMesh->edge(0);
  return data[i];
}

GC_INTERNAL_GENERATE_DATATYPE_OPERATOR_DEFINITIONS(DualEdgeData, dualMesh)

// Data on (real) faces
template <typename T>
DualFaceData<T>::DualFaceData(HalfedgeDual* parentMesh) : dualMesh(parentMesh) {
  if (parentMesh != nullptr) {
    data.resize(parentMesh->nFaces());
  }
}

template <typename T>
DualFaceData<T>::DualFaceData(HalfedgeMesh* parentMesh, T initVal)
    : DualFaceData(parentMesh) {
  fill(initVal);
}

template <typename T>
void DualFaceData<T>::fill(T val) {
  std::fill(data.begin(), data.end(), val);
}

template <typename T>
inline T& DualFaceData<T>::operator[](DualFacePtr f) {
#ifndef NDEBUG
  assert(f->parentMesh == dualMesh->mesh &&
         "Attempted access data with member from wrong dualMesh");
#endif
  unsigned int i = f - dualMesh->face(0);
  return data[i];
}

template <typename T>
inline const T& DualFaceData<T>::operator[](DualFacePtr f) const {
#ifndef NDEBUG
  assert(f->parentMesh == dualMesh->mesh &&
         "Attempted access data with member from wrong dualMesh");
#endif
  unsigned int i = f - dualMesh->face(0);
  return data[i];
}

GC_INTERNAL_GENERATE_DATATYPE_OPERATOR_DEFINITIONS(DualFaceData, dualMesh)

// Data on (real and imaginary) dual halfedges
template <typename T>
DualHalfedgeData<T>::DualHalfedgeData(HalfedgeDual* parentMesh)
    : dualMesh(parentMesh) {
  if (parentMesh != nullptr) {
    data.resize(parentMesh->nHalfedges());
  }
}

template <typename T>
DualHalfedgeData<T>::DualHalfedgeData(HalfedgeMesh* parentMesh, T initVal)
    : DualHalfedgeData(parentMesh) {
  fill(initVal);
}

template <typename T>
void DualHalfedgeData<T>::fill(T val) {
  std::fill(data.begin(), data.end(), val);
}

template <typename T>
inline T& DualHalfedgeData<T>::operator[](DualHalfedgePtr he) {
#ifndef NDEBUG
  assert(he->parentMesh == dualMesh &&
         "Attempted access data with member from wrong dualMesh");
#endif
  unsigned int i = he - dualMesh->halfedge(0);
  return data[i];
}

template <typename T>
inline const T& DualHalfedgeData<T>::operator[](DualHalfedgePtr he) const {
#ifndef NDEBUG
  assert(he->parentMesh == dualMesh &&
         "Attempted access data with member from wrong dualMesh");
#endif
  unsigned int i = he - dualMesh->halfedge(0);
  return data[i];
}

GC_INTERNAL_GENERATE_DATATYPE_OPERATOR_DEFINITIONS(DualHalfedgeData, dualMesh)

}  // namespace geometrycentral