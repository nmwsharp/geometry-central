#pragma once

#include <cassert>

#include "halfedge_data_macros.h"

#include "Eigen/Core"

// === Datatypes which hold data stored on the mesh


namespace geometrycentral {

// Data on vertices
template <typename T> class VertexData {
private:
  HalfedgeMesh* mesh = nullptr;
  std::vector<T> data;

public:
  VertexData() {}
  VertexData(HalfedgeMesh* parentMesh);
  VertexData(HalfedgeMesh* parentMesh, T initVal);
  VertexData(HalfedgeMesh* parentMesh, const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector);
  VertexData(HalfedgeMesh* parentMesh, const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector,
             const VertexData<size_t>& indexer);

  T& operator[](VertexPtr v);
  const T& operator[](VertexPtr v) const;

  void fill(T val);
  Eigen::Matrix<T, Eigen::Dynamic, 1> toVector() const;
  Eigen::Matrix<T, Eigen::Dynamic, 1> toVector(const VertexData<size_t>& indexer) const;
  void fromVector(const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector);
  void fromVector(const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector, const VertexData<size_t>& indexer);

  GC_INTERNAL_GENERATE_DATATYPE_OPERATOR_DECLARATIONS(VertexData)
};

// Data on (real) faces
template <typename T> class FaceData {
private:
  HalfedgeMesh* mesh;
  std::vector<T> data;
  size_t realSize;

public:
  FaceData() {}
  FaceData(HalfedgeMesh* parentMesh);
  FaceData(HalfedgeMesh* parentMesh, T initVal);
  FaceData(HalfedgeMesh* parentMesh, const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector);
  FaceData(HalfedgeMesh* parentMesh, const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector, const FaceData<size_t>& indexer);

  T& operator[](FacePtr f);
  const T& operator[](FacePtr f) const;

  void fill(T val);
  Eigen::Matrix<T, Eigen::Dynamic, 1> toVector() const;
  Eigen::Matrix<T, Eigen::Dynamic, 1> toVector(const FaceData<size_t>& indexer) const;
  void fromVector(const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector);
  void fromVector(const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector, const FaceData<size_t>& indexer);

  GC_INTERNAL_GENERATE_DATATYPE_OPERATOR_DECLARATIONS(FaceData)
};

// Data on edges
template <typename T> class EdgeData {
private:
  HalfedgeMesh* mesh;
  std::vector<T> data;

public:
  EdgeData() {}
  EdgeData(HalfedgeMesh* parentMesh);
  EdgeData(HalfedgeMesh* parentMesh, T initVal);
  EdgeData(HalfedgeMesh* parentMesh, const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector);
  EdgeData(HalfedgeMesh* parentMesh, const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector, const EdgeData<size_t>& indexer);

  T& operator[](EdgePtr e);
  const T& operator[](EdgePtr e) const;

  void fill(T val);
  Eigen::Matrix<T, Eigen::Dynamic, 1> toVector() const;
  Eigen::Matrix<T, Eigen::Dynamic, 1> toVector(const EdgeData<size_t>& indexer) const;
  void fromVector(const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector);
  void fromVector(const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector, const EdgeData<size_t>& indexer);

  GC_INTERNAL_GENERATE_DATATYPE_OPERATOR_DECLARATIONS(EdgeData)
};

// Data on (real and imaginary) halfedges
template <typename T> class HalfedgeData {
private:
  HalfedgeMesh* mesh;
  std::vector<T> data;
  size_t realSize;

public:
  HalfedgeData() {}
  HalfedgeData(HalfedgeMesh* parentMesh);
  HalfedgeData(HalfedgeMesh* parentMesh, T initVal);
  HalfedgeData(HalfedgeMesh* parentMesh, const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector);
  HalfedgeData(HalfedgeMesh* parentMesh, const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector,
               const HalfedgeData<size_t>& indexer);

  T& operator[](HalfedgePtr he);
  const T& operator[](HalfedgePtr he) const;

  void fill(T val);
  Eigen::Matrix<T, Eigen::Dynamic, 1> toVector() const;
  Eigen::Matrix<T, Eigen::Dynamic, 1> toVector(const HalfedgeData<size_t>& indexer) const;
  void fromVector(const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector);
  void fromVector(const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector, const HalfedgeData<size_t>& indexer);

  GC_INTERNAL_GENERATE_DATATYPE_OPERATOR_DECLARATIONS(HalfedgeData)
};

// Data on corners
template <typename T> class CornerData {
private:
  HalfedgeMesh* mesh;
  std::vector<T> data;
  unsigned int realSize;

public:
  CornerData() {}
  CornerData(HalfedgeMesh* parentMesh);
  CornerData(HalfedgeMesh* parentMesh, T initVal);
  CornerData(HalfedgeMesh* parentMesh, const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector);
  CornerData(HalfedgeMesh* parentMesh, const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector,
             const CornerData<size_t>& indexer);

  T& operator[](CornerPtr c);
  const T& operator[](CornerPtr c) const;

  void fill(T val);
  Eigen::Matrix<T, Eigen::Dynamic, 1> toVector() const;
  Eigen::Matrix<T, Eigen::Dynamic, 1> toVector(const CornerData<size_t>& indexer) const;
  void fromVector(const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector);
  void fromVector(const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector, const CornerData<size_t>& indexer);

  GC_INTERNAL_GENERATE_DATATYPE_OPERATOR_DECLARATIONS(CornerData)
};

} // namespace geometrycentral
