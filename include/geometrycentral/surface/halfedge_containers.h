#pragma once

#include "geometrycentral/utilities/dependent_quantity.h"

#include "geometrycentral/surface/halfedge_element_types.h"

#include "Eigen/Core"

#include <cassert>

// === Datatypes which hold data stored on the mesh

namespace geometrycentral {
namespace surface {


// Geneneric datatype, specialized as VertexData (etc) below
// E is the element pointer type (eg Vertex)
// T is the data type that it holds (eg double)
template <typename E, typename T>
class MeshData {
private:
  // The mesh that this data is defined on
  HalfedgeMesh* mesh = nullptr;

  // A default value used to initialize all entries (both on creation, and if the container is expanded due to mesh
  // modification).
  T defaultValue;

  // The raw buffer which holds the data.
  // As a mesh is being modified, data.size() might be larger than the number of elements. Don't attempt any direct
  // access to this buffer.
  std::vector<T> data;

  // Mutability behavior:
  // From the user's point of view, this container can always be accessed with a valid element pointer, no matter what
  // resizing is done. This is implemented by mirroring the callbacks used in the HalfedgeMesh class on resizing; the
  // data<> vector here is always the same size as the corresponding vector of elements in the mesh class.
  // Accessing with an integer index is only meaningful when the backing mesh is compressed.

  // Manage a callback on the mesh object used to keep the container valid on resize events. Should be called once on
  // construction and once on destruction, respectively.
  std::list<std::function<void(size_t)>>::iterator expandCallbackIt;
  std::list<std::function<void(const std::vector<size_t>&)>>::iterator permuteCallbackIt;
  std::list<std::function<void()>>::iterator deleteCallbackIt;
  void registerWithMesh();
  void deregisterWithMesh();

public:
  MeshData();
  MeshData(HalfedgeMesh& parentMesh);
  MeshData(HalfedgeMesh& parentMesh, T initVal);
  MeshData(HalfedgeMesh& parentMesh, const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector);
  MeshData(HalfedgeMesh& parentMesh, const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector,
           const MeshData<E, size_t>& indexer);

  // Rule of 5
  MeshData(const MeshData<E, T>& other);                      // copy constructor
  MeshData(MeshData<E, T>&& other) noexcept;                  // move constructor
  MeshData<E, T>& operator=(const MeshData<E, T>& other);     // copy assignment
  MeshData<E, T>& operator=(MeshData<E, T>&& other) noexcept; // move assignment
  ~MeshData();                                                // destructor

  // Access with an element pointer
  T& operator[](E e);
  const T& operator[](E e) const;

  // Access with an index. Returns the item associated with the i'th element pointer (in the natural ordering).
  // The underlying mesh must be compressed.
  T& operator[](size_t e);
  const T& operator[](size_t e) const;

  // Get the size of the container
  // (note: logical size, like nVertices(), not actual size of vector buffer)
  size_t size() const;

  // Fill with some value
  void fill(T val);

  // Clear out storage.
  // Essentially resets to MeshData<>(), can no longer be used to hold data.
  void clear();

  // Convert to and from vector types
  Eigen::Matrix<T, Eigen::Dynamic, 1> toVector() const;
  Eigen::Matrix<T, Eigen::Dynamic, 1> toVector(const MeshData<E, size_t>& indexer) const;
  void fromVector(const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector);
  void fromVector(const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector, const MeshData<E, size_t>& indexer);

  // Naively reinterpret the data as residing on another mesh, constructing a new container
  MeshData<E, T> reinterpretTo(HalfedgeMesh& targetMesh);
};

// === Typdefs for the usual VertexData<> etc
template <typename T>
using VertexData = MeshData<Vertex, T>;

template <typename T>
using FaceData = MeshData<Face, T>;

template <typename T>
using EdgeData = MeshData<Edge, T>;

template <typename T>
using HalfedgeData = MeshData<Halfedge, T>;

template <typename T>
using CornerData = MeshData<Corner, T>;

template <typename T>
using BoundaryLoopData = MeshData<BoundaryLoop, T>;

} // namespace surface
} // namespace geometrycentral
