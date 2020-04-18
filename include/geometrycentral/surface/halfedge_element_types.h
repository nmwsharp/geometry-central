#pragma once

#include "geometrycentral/utilities/utilities.h"

#include <cstddef>
#include <functional>
#include <iostream>
#include <list>

namespace geometrycentral {
namespace surface {

// === Types and inline methods for the halfedge mesh pointer and datatypes
class HalfedgeMesh;


// === Types for mesh elements
class Vertex;
class Halfedge;
class Corner;
class Edge;
class Face;
class BoundaryLoop;


// === Forward declare iterator set types (pointers may need to return these to
// support range-based for loops, and we cannot foward declare typedef.)

template <typename N>
class NavigationSetBase;
struct VertexAdjacentVertexNavigator;
struct VertexIncomingHalfedgeNavigator;
struct VertexOutgoingHalfedgeNavigator;
struct VertexAdjacentFaceNavigator;
struct VertexAdjacentEdgeNavigator;
struct VertexAdjacentCornerNavigator;
struct FaceAdjacentHalfedgeNavigator;
struct FaceAdjacentVertexNavigator;
struct FaceAdjacentEdgeNavigator;
struct FaceAdjacentFaceNavigator;
struct FaceAdjacentCornerNavigator;
struct BoundaryLoopAdjacentHalfedgeNavigator;
struct BoundaryLoopAdjacentVertexNavigator;
struct BoundaryLoopAdjacentEdgeNavigator;
struct BoundaryLoopAdjacentFaceNavigator;
struct BoundaryLoopAdjacentCornerNavigator;


// === Templated helper functions
// The corresponding specializations are given in halfedge_logic_templates.ipp
// clang-format off

// Enums, which we use to template the base case of all actual element types
enum class ElementType { Vertex = 0, Halfedge, Corner, Edge, Face, BoundaryLoop };

// Current count of this element in the mesh
template <typename E> size_t nElements(HalfedgeMesh* mesh) { return INVALID_IND; }

// Capacity of element type in mesh (containers should be at least this big before next resize)
template <typename E> size_t elementCapacity(HalfedgeMesh* mesh) { return INVALID_IND; }

// Canonical index for this element (not always an enumeration)
template <typename E> size_t dataIndexOfElement(HalfedgeMesh* mesh, E e) { return INVALID_IND; }

// The set type used to iterate over all elements of this type
template <typename E> struct ElementSetType { typedef std::tuple<> type; }; // nonsense default value

template <typename E> typename ElementSetType<E>::type iterateElements(HalfedgeMesh* mesh) { return INVALID_IND; }

template <typename E> std::list<std::function<void(size_t)>>& getExpandCallbackList(HalfedgeMesh* mesh);

template <typename E> std::list<std::function<void(const std::vector<size_t>&)>>& getPermuteCallbackList(HalfedgeMesh* mesh);

template <typename E> std::string typeShortName() { return "X"; }
// clang-format on

// ==========================================================
// ================      Base Element      ==================
// ==========================================================


// Need to pre-declare template friends to avoid generating a non-template version
// (see https://isocpp.org/wiki/faq/templates#template-friends)
template <typename T>
class Element;
} // namespace surface
} // namespace geometrycentral

namespace std {
template <typename T>
std::ostream& operator<<(std::ostream& o, const geometrycentral::surface::Element<T>& x);
template <typename T>
std::string to_string(const geometrycentral::surface::Element<T>& e);
} // namespace std

namespace geometrycentral {
namespace surface {

// Forward-declare dynamic equivalent so we can declare a conversion constructor
template <typename S>
class DynamicElement;

// == Base type for shared logic between elements.
//
// This class uses the "curiously recurring template pattern" (CRTP) because it allows to implement more shared
// functionality in this template. Instantiations of the template will be templted on its child types, like `Vertex :
// public Element<Vertex>`. Because the parent class will know its child type at compile time, it can customize
// functionality based on that child type using the helpers above. This is essentially "compile time polymorphism",
// which allows us to share common functionality without paying a virtual function runtime cost.
template <typename T>
class Element {

public:
  Element();                               // construct an empty (null) element
  Element(HalfedgeMesh* mesh, size_t ind); // construct pointing to the i'th element of that type on a mesh.
  Element(const DynamicElement<T>& e);     // construct from a dynamic element of matching type

  inline bool operator==(const Element<T>& other) const;
  inline bool operator!=(const Element<T>& other) const;
  inline bool operator>(const Element<T>& other) const;
  inline bool operator>=(const Element<T>& other) const;
  inline bool operator<(const Element<T>& other) const;
  inline bool operator<=(const Element<T>& other) const;

  // Get the "index" associated with the element.
  // Note that these are not always a dense enumeration, and generally should not be accessed by "users" unless you are
  // monkeying around the HalfedgeMesh datastructure in some deep and scary way. Generally prefer
  // `HalfedgeMesh::getVertexIndices()` (etc) if you are looking for a set of indices for a linear algebra problem or
  // something.
  size_t getIndex() const;

  // Get the halfedge mesh on which the element is defined.
  HalfedgeMesh* getMesh() const;

  bool isDead() const;

protected:
  HalfedgeMesh* mesh = nullptr;
  size_t ind = INVALID_IND;

  // Friends
  friend std::ostream& std::operator<<<>(std::ostream& output, const Element<T>& e);
  friend struct std::hash<Element<T>>;
  friend std::string std::to_string<>(const Element<T>& e);
};
} // namespace surface
} // namespace geometrycentral

namespace std {
template <typename T>
std::ostream& operator<<(std::ostream& output, const geometrycentral::surface::Element<T>& e);
template <typename T>
std::string to_string(const geometrycentral::surface::Element<T>& e);
} // namespace std


namespace geometrycentral {
namespace surface {

// The equivalent dynamic pointers. These should be rarely used, but are guaranteed to be preserved through _all_ mesh
// operations, including compress().
template <typename S>
class DynamicElement : public S {
public:
  DynamicElement();                                 // construct an empty (null) element
  DynamicElement(HalfedgeMesh* mesh_, size_t ind_); // construct from an index as usual
  DynamicElement(const S& e);                       // construct from a non-dynamic element


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


// == Base range iterator
// All range iterators have the form "advance through indices, skipping invalid elements". The two classes below
// encapsulate that functionality, allowing us to just specify the element type and "valid" function for each.
template <typename F>
class RangeIteratorBase {

public:
  RangeIteratorBase(HalfedgeMesh* mesh_, size_t iStart_, size_t iEnd_);
  const RangeIteratorBase& operator++();
  bool operator==(const RangeIteratorBase& other) const;
  bool operator!=(const RangeIteratorBase& other) const;
  typename F::Etype operator*() const;

private:
  HalfedgeMesh* mesh;
  size_t iCurr;
  size_t iEnd;
};

template <typename F>
class RangeSetBase {
public:
  RangeSetBase(HalfedgeMesh* mesh_, size_t iStart_, size_t iEnd_);
  RangeIteratorBase<F> begin() const;
  RangeIteratorBase<F> end() const;

private:
  HalfedgeMesh* mesh;
  size_t iStart, iEnd;
};


// ==========================================================
// ================        Vertex          ==================
// ==========================================================

class Vertex : public Element<Vertex> {
public:
  // Constructors
  // inheriting constructor would work here, and replace the constructors below, but gcc-5 erroneously rejects the combo
  // with CRTP :( perhaps resurrect here and in other elements below once gcc-5 is sufficiently old
  // using Element<Vertex>::Element;
  Vertex();                                // construct an empty (null) element
  Vertex(HalfedgeMesh* mesh, size_t ind);  // construct pointing to the i'th element of that type on a mesh.
  Vertex(const DynamicElement<Vertex>& e); // construct from a dynamic element of matching type

  // Navigators
  Halfedge halfedge() const;
  Corner corner() const;
  bool isDead() const;

  // Properties
  bool isBoundary() const;
  size_t degree() const;
  size_t faceDegree() const;

  // Iterators
  NavigationSetBase<VertexAdjacentVertexNavigator> adjacentVertices() const;
  NavigationSetBase<VertexIncomingHalfedgeNavigator> incomingHalfedges() const;
  NavigationSetBase<VertexOutgoingHalfedgeNavigator> outgoingHalfedges() const;
  NavigationSetBase<VertexAdjacentCornerNavigator> adjacentCorners() const;
  NavigationSetBase<VertexAdjacentEdgeNavigator> adjacentEdges() const;
  NavigationSetBase<VertexAdjacentFaceNavigator> adjacentFaces() const;
};

using DynamicVertex = DynamicElement<Vertex>;

// == Range iterators

// All vertices
struct VertexRangeF {
  static bool elementOkay(const HalfedgeMesh& mesh, size_t ind);
  typedef Vertex Etype;
};
typedef RangeSetBase<VertexRangeF> VertexSet;

// ==========================================================
// ================        Halfedge        ==================
// ==========================================================

class Halfedge : public Element<Halfedge> {
public:
  // Constructors
  Halfedge();                                  // construct an empty (null) element
  Halfedge(HalfedgeMesh* mesh, size_t ind);    // construct pointing to the i'th element of that type on a mesh.
  Halfedge(const DynamicElement<Halfedge>& e); // construct from a dynamic element of matching type

  // Navigators
  Halfedge twin() const;
  Halfedge next() const;
  Corner corner() const;
  Vertex vertex() const;
  Edge edge() const;
  Face face() const;
  bool isDead() const;

  // Super-navigators
  Halfedge prevOrbitFace() const;
  Halfedge prevOrbitVertex() const;

  // Properties
  bool isInterior() const;
};

using DynamicHalfedge = DynamicElement<Halfedge>;

// == Range iterators

// All halfedges
struct HalfedgeRangeF {
  static bool elementOkay(const HalfedgeMesh& mesh, size_t ind);
  typedef Halfedge Etype;
};
typedef RangeSetBase<HalfedgeRangeF> HalfedgeSet;

// Interior halfedges
struct HalfedgeInteriorRangeF {
  static bool elementOkay(const HalfedgeMesh& mesh, size_t ind);
  typedef Halfedge Etype;
};
typedef RangeSetBase<HalfedgeInteriorRangeF> HalfedgeInteriorSet;

// Exterior halfedges
struct HalfedgeExteriorRangeF {
  static bool elementOkay(const HalfedgeMesh& mesh, size_t ind);
  typedef Halfedge Etype;
};
typedef RangeSetBase<HalfedgeExteriorRangeF> HalfedgeExteriorSet;

// ==========================================================
// ================        Corner          ==================
// ==========================================================

// Implmentation note: The `ind` parameter for a corner will be the index of a halfedge, which should always be real.

class Corner : public Element<Corner> {
public:
  // Constructors
  Corner();                                // construct an empty (null) element
  Corner(HalfedgeMesh* mesh, size_t ind);  // construct pointing to the i'th element of that type on a mesh.
  Corner(const DynamicElement<Corner>& e); // construct from a dynamic element of matching type

  // Navigators
  Halfedge halfedge() const;
  Vertex vertex() const;
  Face face() const;
  bool isDead() const;
};

using DynamicCorner = DynamicElement<Corner>;

// == Range iterators

// All corners
struct CornerRangeF {
  static bool elementOkay(const HalfedgeMesh& mesh, size_t ind);
  typedef Corner Etype;
};
typedef RangeSetBase<CornerRangeF> CornerSet;

// ==========================================================
// ================          Edge          ==================
// ==========================================================

class Edge : public Element<Edge> {
public:
  // Constructors
  Edge();                               // construct an empty (null) element
  Edge(HalfedgeMesh* mesh, size_t ind); // construct pointing to the i'th element of that type on a mesh.
  Edge(const DynamicElement<Edge>& e);  // construct from a dynamic element of matching type

  // Navigators
  Halfedge halfedge() const;
  bool isDead() const;

  // Properties
  bool isBoundary() const;
};

using DynamicEdge = DynamicElement<Edge>;

// == Range iterators

// All edges
struct EdgeRangeF {
  static bool elementOkay(const HalfedgeMesh& mesh, size_t ind);
  typedef Edge Etype;
};
typedef RangeSetBase<EdgeRangeF> EdgeSet;

// ==========================================================
// ================          Face          ==================
// ==========================================================

// Implmentation note: The `ind` parameter for a face might correspond to a boundary loop. The boundary loops have face
// IDs which are at the very end of the face buffer, but can still index in to face-valued arrays/functions in
// HalfedgeMesh (they _cannot_ index in to FaceData<> containers).

class Face : public Element<Face> {
public:
  // Constructors
  Face();                               // construct an empty (null) element
  Face(HalfedgeMesh* mesh, size_t ind); // construct pointing to the i'th element of that type on a mesh.
  Face(const DynamicElement<Face>& e);  // construct from a dynamic element of matching type

  // Navigators
  Halfedge halfedge() const;
  BoundaryLoop asBoundaryLoop() const;
  bool isDead() const;

  // Properties
  bool isBoundaryLoop() const;
  bool isTriangle() const;
  size_t degree() const;

  // Iterators
  NavigationSetBase<FaceAdjacentVertexNavigator> adjacentVertices() const;
  NavigationSetBase<FaceAdjacentHalfedgeNavigator> adjacentHalfedges() const;
  NavigationSetBase<FaceAdjacentCornerNavigator> adjacentCorners() const;
  NavigationSetBase<FaceAdjacentEdgeNavigator> adjacentEdges() const;
  NavigationSetBase<FaceAdjacentFaceNavigator> adjacentFaces() const;
};

using DynamicFace = DynamicElement<Face>;

// == Range iterators

// All faces
struct FaceRangeF {
  static bool elementOkay(const HalfedgeMesh& mesh, size_t ind);
  typedef Face Etype;
};
typedef RangeSetBase<FaceRangeF> FaceSet;

// ==========================================================
// ================     Boundary Loop      ==================
// ==========================================================

// Implementation note: the `ind` parameter for a boundary loop is index from the back of the face index space, from [0,
// nBoundaryLoopFillCount).

class BoundaryLoop : public Element<BoundaryLoop> {
public:
  // Constructors
  BoundaryLoop();                                      // construct an empty (null) element
  BoundaryLoop(HalfedgeMesh* mesh, size_t ind);        // construct pointing to the i'th element of that type on a mesh.
  BoundaryLoop(const DynamicElement<BoundaryLoop>& e); // construct from a dynamic element of matching type

  Halfedge halfedge() const;
  Face asFace() const;
  bool isDead() const;

  // Properties
  size_t degree() const;

  // Iterators
  NavigationSetBase<BoundaryLoopAdjacentVertexNavigator> adjacentVertices() const;
  NavigationSetBase<BoundaryLoopAdjacentHalfedgeNavigator> adjacentHalfedges() const;
  NavigationSetBase<BoundaryLoopAdjacentEdgeNavigator> adjacentEdges() const;
};

using DynamicBoundaryLoop = DynamicElement<BoundaryLoop>;

// == Range iterators

// All boundary loops
struct BoundaryLoopRangeF {
  static bool elementOkay(const HalfedgeMesh& mesh, size_t ind);
  typedef BoundaryLoop Etype;
};
typedef RangeSetBase<BoundaryLoopRangeF> BoundaryLoopSet;


} // namespace surface
} // namespace geometrycentral
