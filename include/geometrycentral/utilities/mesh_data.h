#pragma once

#include "geometrycentral/utilities/element.h"

// Need to include all of the relevant mesh types, so that we get forward declarations of template specializations which
// are used below.
#include "geometrycentral/surface/halfedge_element_types.h"

#include <Eigen/Core>
#include <cassert>

// === Datatypes which hold data stored on the mesh

namespace geometrycentral {


// Geneneric datatype, specialized as VertexData (etc) below
// E is the element pointer type (eg Vertex)
// T is the data type that it holds (eg double)
template <typename E, typename T>
class MeshData {
public:
  // The type of the raw buffer which holds the data.
  // Note that this is Eigen's vector type. In particular, it should (???) implement alignment policies which may
  // improve vectorization and performance.
  using DATA_T = Eigen::Matrix<T, Eigen::Dynamic, 1>;

protected:
  // The mesh that this data is defined on
  using ParentMeshT = typename E::ParentMeshT;
  ParentMeshT* mesh = nullptr;

  // A default value used to initialize all entries (both on creation, and if the container is expanded due to mesh
  // modification).
  T defaultValue = T();

  // The raw buffer which holds the data.
  // As a mesh is being modified, data.size() might be larger than the number of elements. Don't attempt any direct
  // access to this buffer.
  DATA_T data;

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
  MeshData(ParentMeshT& parentMesh);
  MeshData(ParentMeshT& parentMesh, T initVal);
  // here `vector` should be a length-nElements dense array (which may be different from the raw internal storage of
  // this class when the mesh is compressed)
  MeshData(ParentMeshT& parentMesh, const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector);
  MeshData(ParentMeshT& parentMesh, const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector,
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

  // Raw access to the underlying buffer
  DATA_T& raw();
  const DATA_T& raw() const;

  // Access to the underlying mesh object
  ParentMeshT* getMesh() const;

  // Fill with some value
  void fill(T val);

  // Clear out storage.
  // Essentially resets to MeshData<>(), can no longer be used to hold data.
  void clear();

  // Convert to and from (Eigen) vector types
  Eigen::Matrix<T, Eigen::Dynamic, 1> toVector() const;
  Eigen::Matrix<T, Eigen::Dynamic, 1> toVector(const MeshData<E, size_t>& indexer) const;
  void fromVector(const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector);
  void fromVector(const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector, const MeshData<E, size_t>& indexer);

  // Naively reinterpret the data as residing on another mesh, constructing a new container
  MeshData<E, T> reinterpretTo(ParentMeshT& targetMesh) const;

  void setDefault(T newDefault);
  T getDefault() const;
};

// === Arithmetic overloads ===

// Allow things like adding meshdata containers, multiplying by a scalar, etc.

// helper used in many functions below
template <typename E, typename T, typename U>
void checkMeshCompatible(const MeshData<E, T>& lhs, const MeshData<E, U>& rhs);
template <typename E, typename T>
void checkMeshValid(const MeshData<E, T>& val);

// Yes, macros. The alternative seems to be _a lot_ of duplicated code, and I don't know of anyway to template on
// the operator in C++. The macros below define all the necessary functions for each unary/binary operator. At the
// bottom we instantiate all the common arithmetic and logical operators. We don't bother separating declaration from
// definition, since it's not human-readable anyway.

#define GC_INTERNAL_MESHDATA_UNARY(OP)                                                                                 \
  template <typename E, typename T>                                                                                    \
  MeshData<E, decltype(OP std::declval<T>())> operator OP(const MeshData<E, T>& lhs) {                                 \
    checkMeshValid(lhs);                                                                                               \
    MeshData<E, decltype(OP std::declval<T>())> result(*lhs.getMesh());                                                \
    for (Eigen::Index i = 0; i < lhs.raw().rows(); i++) result.raw()[i] = OP lhs.raw()[i];                             \
    return result;                                                                                                     \
  }


#define GC_INTERNAL_MESHDATA_BINARY_BASIC(OP)                                                                          \
  template <typename E, typename T, typename U>                                                                        \
  MeshData<E, decltype(std::declval<T>() OP std::declval<U>())> operator OP(const MeshData<E, T>& lhs,                 \
                                                                            const MeshData<E, U>& rhs) {               \
    checkMeshCompatible(lhs, rhs);                                                                                     \
    MeshData<E, decltype(std::declval<T>() OP std::declval<U>())> result(*lhs.getMesh());                              \
    for (Eigen::Index i = 0; i < lhs.raw().rows(); i++) result.raw()[i] = lhs.raw()[i] OP rhs.raw()[i];                \
    return result;                                                                                                     \
  }                                                                                                                    \
  template <typename E, typename T, typename U>                                                                        \
  MeshData<E, decltype(std::declval<T>() OP std::declval<U>())> operator OP(const MeshData<E, T>& lhs, const U& rhs) { \
    checkMeshValid(lhs);                                                                                               \
    MeshData<E, decltype(std::declval<T>() OP std::declval<U>())> result(*lhs.getMesh());                              \
    for (Eigen::Index i = 0; i < lhs.raw().rows(); i++) result.raw()[i] = lhs.raw()[i] OP rhs;                         \
    return result;                                                                                                     \
  }                                                                                                                    \
  template <typename E, typename T, typename U>                                                                        \
  MeshData<E, decltype(std::declval<T>() OP std::declval<U>())> operator OP(const T& lhs, const MeshData<E, U>& rhs) { \
    checkMeshValid(rhs);                                                                                               \
    MeshData<E, decltype(std::declval<T>() OP std::declval<U>())> result(*rhs.getMesh());                              \
    for (Eigen::Index i = 0; i < rhs.raw().rows(); i++) result.raw()[i] = lhs OP rhs.raw()[i];                         \
    return result;                                                                                                     \
  }

// for operators where we can also call e.g. += on each scalar (for e.g. &&= this doesn't exist)
#define GC_INTERNAL_MESHDATA_BINARY_EQUALS_ASSIGN(OP)                                                                  \
  GC_INTERNAL_MESHDATA_BINARY_BASIC(OP)                                                                                \
  template <typename E, typename T, typename U>                                                                        \
  MeshData<E, decltype(std::declval<T>() OP std::declval<U>())>& operator OP##=(MeshData<E, T>& lhs,                   \
                                                                                const MeshData<E, U>& rhs) {           \
    checkMeshCompatible(lhs, rhs);                                                                                     \
    for (Eigen::Index i = 0; i < lhs.raw().rows(); i++) lhs.raw()[i] OP## = rhs.raw()[i];                              \
    return lhs;                                                                                                        \
  }                                                                                                                    \
  template <typename E, typename T, typename U>                                                                        \
  MeshData<E, decltype(std::declval<T>() OP std::declval<U>())>& operator OP##=(MeshData<E, T>& lhs, const U& rhs) {   \
    checkMeshValid(lhs);                                                                                               \
    for (Eigen::Index i = 0; i < lhs.raw().rows(); i++) lhs.raw()[i] OP## = rhs;                                       \
    return lhs;                                                                                                        \
  }


GC_INTERNAL_MESHDATA_UNARY(+)
GC_INTERNAL_MESHDATA_UNARY(-)
GC_INTERNAL_MESHDATA_UNARY(!)
GC_INTERNAL_MESHDATA_UNARY(~)

GC_INTERNAL_MESHDATA_BINARY_EQUALS_ASSIGN(+)
GC_INTERNAL_MESHDATA_BINARY_EQUALS_ASSIGN(-)
GC_INTERNAL_MESHDATA_BINARY_EQUALS_ASSIGN(*)
GC_INTERNAL_MESHDATA_BINARY_EQUALS_ASSIGN(/)
GC_INTERNAL_MESHDATA_BINARY_EQUALS_ASSIGN(%)
GC_INTERNAL_MESHDATA_BINARY_EQUALS_ASSIGN(&)
GC_INTERNAL_MESHDATA_BINARY_EQUALS_ASSIGN(|)
GC_INTERNAL_MESHDATA_BINARY_EQUALS_ASSIGN(^)
GC_INTERNAL_MESHDATA_BINARY_EQUALS_ASSIGN(<<)
GC_INTERNAL_MESHDATA_BINARY_EQUALS_ASSIGN(>>)
GC_INTERNAL_MESHDATA_BINARY_BASIC(&&)
GC_INTERNAL_MESHDATA_BINARY_BASIC(||)

#undef GC_INTERNAL_MESHDATA_UNARY
#undef GC_INTERNAL_MESHDATA_BINARY_BASIC
#undef GC_INTERNAL_MESHDATA_BINARY_EQUALS_ASSIGN

} // namespace geometrycentral

#include "geometrycentral/utilities/mesh_data.ipp"
