#pragma once

#include <cassert>

#include "halfedge_data_macros.h"

// Data on dual vertices
template <typename T>
class DualVertexData {
 private:
  HalfedgeDual* dualMesh = nullptr;
  std::vector<T> data;

 public:
  DualVertexData() {}
  DualVertexData(HalfedgeDual* parentMesh);
  DualVertexData(HalfedgeMesh* parentMesh, T initVal);

  T& operator[](DualVertexPtr v);
  const T& operator[](DualVertexPtr v) const;

  void fill(T val);

  GC_INTERNAL_GENERATE_DATATYPE_OPERATOR_DECLARATIONS(DualVertexData)
};

// Data on (real) faces
template <typename T>
class DualFaceData {
 private:
  HalfedgeDual* dualMesh;
  std::vector<T> data;

 public:
  DualFaceData() {}
  DualFaceData(HalfedgeDual* parentMesh);
  DualFaceData(HalfedgeMesh* parentMesh, T initVal);

  T& operator[](DualFacePtr f);
  const T& operator[](DualFacePtr f) const;

  void fill(T val);

  GC_INTERNAL_GENERATE_DATATYPE_OPERATOR_DECLARATIONS(DualFaceData)
};

// Data on edges
template <typename T>
class DualEdgeData {
 private:
  HalfedgeDual* dualMesh;
  std::vector<T> data;

 public:
  DualEdgeData() {}
  DualEdgeData(HalfedgeDual* parentMesh);
  DualEdgeData(HalfedgeMesh* parentMesh, T initVal);

  T& operator[](DualEdgePtr e);
  const T& operator[](DualEdgePtr e) const;

  void fill(T val);

  GC_INTERNAL_GENERATE_DATATYPE_OPERATOR_DECLARATIONS(DualEdgeData)
};

// Data on (real and imaginary) halfedges
template <typename T>
class DualHalfedgeData {
 private:
  HalfedgeDual* dualMesh;
  std::vector<T> data;

 public:
  DualHalfedgeData() {}
  DualHalfedgeData(HalfedgeDual* parentMesh);
  DualHalfedgeData(HalfedgeMesh* parentMesh, T initVal);

  T& operator[](DualHalfedgePtr he);
  const T& operator[](DualHalfedgePtr he) const;

  void fill(T val);

  GC_INTERNAL_GENERATE_DATATYPE_OPERATOR_DECLARATIONS(DualHalfedgeData)
};
