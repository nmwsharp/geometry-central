#pragma once

#include <cassert>

#include "Eigen/Core"

// === Datatypes which hold data stored on the mesh


namespace geometrycentral {

// Geneneric datatype, specialized as VertexData (etc) below
// E is the element pointer type (eg VertexPtr)
// T is the data type that it holds (eg double)
template <typename E, typename T>
class MeshData {
private:
  HalfedgeMesh* mesh = nullptr;
  std::vector<T> data;
  T defaultValue;

public:
  MeshData() {}
  MeshData(HalfedgeMesh* parentMesh);
  MeshData(HalfedgeMesh* parentMesh, T initVal);
  MeshData(HalfedgeMesh* parentMesh, const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector);
  MeshData(HalfedgeMesh* parentMesh, const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector,
           const MeshData<E, size_t>& indexer);

  // Acess with an element pointer
  T& operator[](E e);
  const T& operator[](E e) const;

  // Access with an index
  T& operator[](size_t e);
  const T& operator[](size_t e) const;

  // Get the size of the container
  // (note: logical size, like nVertices(), not actual size of vector buffer)
  size_t size() const;

  // Fill with some value
  void fill(T val);

  // Convert to and from vector types
  Eigen::Matrix<T, Eigen::Dynamic, 1> toVector() const;
  Eigen::Matrix<T, Eigen::Dynamic, 1> toVector(const MeshData<E, size_t>& indexer) const;
  void fromVector(const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector);
  void fromVector(const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector, const MeshData<E, size_t>& indexer);

};

// === Typdefs for the usual VertexData<> etc
template <typename T>
using VertexData = MeshData<VertexPtr, T>;

template <typename T>
using FaceData = MeshData<FacePtr, T>;

template <typename T>
using EdgeData = MeshData<EdgePtr, T>;

template <typename T>
using HalfedgeData = MeshData<HalfedgePtr, T>;

template <typename T>
using CornerData = MeshData<CornerPtr, T>;

} // namespace geometrycentral
