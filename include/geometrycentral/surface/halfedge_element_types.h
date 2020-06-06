#pragma once

#include "geometrycentral/utilities/element.h"
#include "geometrycentral/utilities/element_iterators.h"
#include "geometrycentral/utilities/mesh_data.h"
#include "geometrycentral/utilities/utilities.h"

#include <cstddef>
#include <functional>
#include <iostream>
#include <list>

namespace geometrycentral {
namespace surface {

// === Types and inline methods for the halfedge mesh pointer and datatypes
class HalfedgeMesh;

class Vertex;
class Halfedge;
class Corner;
class Edge;
class Face;
class BoundaryLoop;

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

// ==========================================================
// ================        Vertex          ==================
// ==========================================================

class Vertex : public Element<Vertex, HalfedgeMesh> {
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
  typedef HalfedgeMesh ParentMeshT;
};
typedef RangeSetBase<VertexRangeF> VertexSet;


// ==========================================================
// ================        Halfedge        ==================
// ==========================================================

class Halfedge : public Element<Halfedge, HalfedgeMesh> {
public:
  // Constructors
  Halfedge();                                  // construct an empty (null) element
  Halfedge(HalfedgeMesh* mesh, size_t ind);    // construct pointing to the i'th element of that type on a mesh.
  Halfedge(const DynamicElement<Halfedge>& e); // construct from a dynamic element of matching type

  // Navigators
  Halfedge twin() const;
  Halfedge sibling() const;
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
  typedef HalfedgeMesh ParentMeshT;
};
typedef RangeSetBase<HalfedgeRangeF> HalfedgeSet;

// Interior halfedges
struct HalfedgeInteriorRangeF {
  static bool elementOkay(const HalfedgeMesh& mesh, size_t ind);
  typedef Halfedge Etype;
  typedef HalfedgeMesh ParentMeshT;
};
typedef RangeSetBase<HalfedgeInteriorRangeF> HalfedgeInteriorSet;

// Exterior halfedges
struct HalfedgeExteriorRangeF {
  static bool elementOkay(const HalfedgeMesh& mesh, size_t ind);
  typedef Halfedge Etype;
  typedef HalfedgeMesh ParentMeshT;
};
typedef RangeSetBase<HalfedgeExteriorRangeF> HalfedgeExteriorSet;


// ==========================================================
// ================        Corner          ==================
// ==========================================================

// Implmentation note: The `ind` parameter for a corner will be the index of a halfedge, which should always be real.

class Corner : public Element<Corner, HalfedgeMesh> {
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
  typedef HalfedgeMesh ParentMeshT;
};
typedef RangeSetBase<CornerRangeF> CornerSet;


// ==========================================================
// ================          Edge          ==================
// ==========================================================

class Edge : public Element<Edge, HalfedgeMesh> {
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
  bool isManifold() const;
  //size_t degree() const;
};

using DynamicEdge = DynamicElement<Edge>;

// == Range iterators

// All edges
struct EdgeRangeF {
  static bool elementOkay(const HalfedgeMesh& mesh, size_t ind);
  typedef Edge Etype;
  typedef HalfedgeMesh ParentMeshT;
};
typedef RangeSetBase<EdgeRangeF> EdgeSet;


// ==========================================================
// ================          Face          ==================
// ==========================================================

// Implmentation note: The `ind` parameter for a face might correspond to a boundary loop. The boundary loops have face
// IDs which are at the very end of the face buffer, but can still index in to face-valued arrays/functions in
// HalfedgeMesh (they _cannot_ index in to FaceData<> containers).

class Face : public Element<Face, HalfedgeMesh> {
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
  typedef HalfedgeMesh ParentMeshT;
};
typedef RangeSetBase<FaceRangeF> FaceSet;


// ==========================================================
// ================     Boundary Loop      ==================
// ==========================================================

// Implementation note: the `ind` parameter for a boundary loop is index from the back of the face index space, from [0,
// nBoundaryLoopFillCount).

class BoundaryLoop : public Element<BoundaryLoop, HalfedgeMesh> {
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
  typedef HalfedgeMesh ParentMeshT;
};
typedef RangeSetBase<BoundaryLoopRangeF> BoundaryLoopSet;


// ==========================================================
// ===============   Navigation Iterators   =================
// ==========================================================

// == Vertex

// Helper class to store some special extra state for the vertex iterators
struct VertexNeighborIteratorState {

  VertexNeighborIteratorState(Halfedge currHe);
  VertexNeighborIteratorState(HalfedgeMesh* mesh, Vertex v);

  bool useArray = false;

  // if useArray == false, this is populated
  Halfedge currHe = Halfedge();

  // if useArray == true, these are is populated
  HalfedgeMesh* mesh = nullptr;
  size_t degree = INVALID_IND;
  size_t indStart = INVALID_IND;
  size_t currInd = INVALID_IND;

  Halfedge getCurrHalfedge() const;
  void advance();
  bool operator==(const VertexNeighborIteratorState& rhs) const;
};

// Adjacent vertices
struct VertexAdjacentVertexNavigator {
  void advance();
  bool isValid() const;
  typedef VertexNeighborIteratorState Etype;
  Etype currE;
  typedef Vertex Rtype;
  Rtype getCurrent() const;
};

// Adjacent incoming halfedges
struct VertexIncomingHalfedgeNavigator {
  void advance();
  bool isValid() const;
  typedef VertexNeighborIteratorState Etype;
  Etype currE;
  typedef Halfedge Rtype;
  Rtype getCurrent() const;
};

// Adjacent outgoing halfedges
struct VertexOutgoingHalfedgeNavigator {
  void advance();
  bool isValid() const;
  typedef VertexNeighborIteratorState Etype;
  Etype currE;
  typedef Halfedge Rtype;
  Rtype getCurrent() const;
};

// Adjacent corners
struct VertexAdjacentCornerNavigator {
  void advance();
  bool isValid() const;
  typedef VertexNeighborIteratorState Etype;
  Etype currE;
  typedef Corner Rtype;
  Rtype getCurrent() const;
};


// Adjacent edges
struct VertexAdjacentEdgeNavigator {
  void advance();
  bool isValid() const;
  typedef VertexNeighborIteratorState Etype;
  Etype currE;
  typedef Edge Rtype;
  Rtype getCurrent() const;
};

// Adjacent faces
struct VertexAdjacentFaceNavigator {
  void advance();
  bool isValid() const;
  typedef VertexNeighborIteratorState Etype;
  Etype currE;
  typedef Face Rtype;
  Rtype getCurrent() const;
};


// == Face

// Adjacent vertices
struct FaceAdjacentVertexNavigator {
  void advance();
  bool isValid() const;
  typedef Halfedge Etype;
  Etype currE;
  typedef Vertex Rtype;
  Rtype getCurrent() const;
};

// Adjacent halfedge
struct FaceAdjacentHalfedgeNavigator {
  void advance();
  bool isValid() const;
  typedef Halfedge Etype;
  Etype currE;
  typedef Halfedge Rtype;
  Rtype getCurrent() const;
};

// Adjacent corner
struct FaceAdjacentCornerNavigator {
  void advance();
  bool isValid() const;
  typedef Halfedge Etype;
  Etype currE;
  typedef Corner Rtype;
  Rtype getCurrent() const;
};

// Adjacent edge
struct FaceAdjacentEdgeNavigator {
  void advance();
  bool isValid() const;
  typedef Halfedge Etype;
  Etype currE;
  typedef Edge Rtype;
  Rtype getCurrent() const;
};

// Adjacent face
struct FaceAdjacentFaceNavigator {
  void advance();
  bool isValid() const;
  typedef Halfedge Etype;
  Etype currE;
  typedef Face Rtype;
  Rtype getCurrent() const;
};


// == Faces

// Adjacent vertex
struct BoundaryLoopAdjacentVertexNavigator {
  void advance();
  bool isValid() const;
  typedef Halfedge Etype;
  Etype currE;
  typedef Vertex Rtype;
  Rtype getCurrent() const;
};

// Adjacent halfedge
struct BoundaryLoopAdjacentHalfedgeNavigator {
  void advance();
  bool isValid() const;
  typedef Halfedge Etype;
  Etype currE;
  typedef Halfedge Rtype;
  Rtype getCurrent() const;
};

// Adjacent edge
struct BoundaryLoopAdjacentEdgeNavigator {
  void advance();
  bool isValid() const;
  typedef Halfedge Etype;
  Etype currE;
  typedef Edge Rtype;
  Rtype getCurrent() const;
};


} // namespace surface

// Declare specializations of the logic templates. This is important, because these need to be declared before any of
// the templates using them are instantiated.

// clang-format off

template<> inline size_t nElements<surface::Vertex       >(surface::HalfedgeMesh* mesh); 
template<> inline size_t nElements<surface::Face         >(surface::HalfedgeMesh* mesh); 
template<> inline size_t nElements<surface::Edge         >(surface::HalfedgeMesh* mesh); 
template<> inline size_t nElements<surface::Halfedge     >(surface::HalfedgeMesh* mesh); 
template<> inline size_t nElements<surface::Corner       >(surface::HalfedgeMesh* mesh); 
template<> inline size_t nElements<surface::BoundaryLoop >(surface::HalfedgeMesh* mesh); 

template<> inline size_t elementCapacity<surface::Vertex      >(surface::HalfedgeMesh* mesh);
template<> inline size_t elementCapacity<surface::Face        >(surface::HalfedgeMesh* mesh);
template<> inline size_t elementCapacity<surface::Edge        >(surface::HalfedgeMesh* mesh);
template<> inline size_t elementCapacity<surface::Halfedge    >(surface::HalfedgeMesh* mesh);
template<> inline size_t elementCapacity<surface::Corner      >(surface::HalfedgeMesh* mesh);
template<> inline size_t elementCapacity<surface::BoundaryLoop>(surface::HalfedgeMesh* mesh);

template<> inline size_t dataIndexOfElement<surface::Vertex          >(surface::HalfedgeMesh* mesh, surface::Vertex e           );
template<> inline size_t dataIndexOfElement<surface::Face            >(surface::HalfedgeMesh* mesh, surface::Face e             );
template<> inline size_t dataIndexOfElement<surface::Edge            >(surface::HalfedgeMesh* mesh, surface::Edge e             );
template<> inline size_t dataIndexOfElement<surface::Halfedge        >(surface::HalfedgeMesh* mesh, surface::Halfedge e         );
template<> inline size_t dataIndexOfElement<surface::Corner          >(surface::HalfedgeMesh* mesh, surface::Corner e           );
template<> inline size_t dataIndexOfElement<surface::BoundaryLoop    >(surface::HalfedgeMesh* mesh, surface::BoundaryLoop e     );

template<> struct ElementSetType<surface::Vertex        >   { typedef surface::VertexSet       type; };
template<> struct ElementSetType<surface::Face          >   { typedef surface::FaceSet         type; };
template<> struct ElementSetType<surface::Edge          >   { typedef surface::EdgeSet         type; };
template<> struct ElementSetType<surface::Halfedge      >   { typedef surface::HalfedgeSet     type; };
template<> struct ElementSetType<surface::Corner        >   { typedef surface::CornerSet       type; };
template<> struct ElementSetType<surface::BoundaryLoop  >   { typedef surface::BoundaryLoopSet type; };

template<> inline surface::VertexSet         iterateElements<surface::Vertex      >(surface::HalfedgeMesh* mesh);
template<> inline surface::HalfedgeSet       iterateElements<surface::Halfedge    >(surface::HalfedgeMesh* mesh);
template<> inline surface::CornerSet         iterateElements<surface::Corner      >(surface::HalfedgeMesh* mesh);
template<> inline surface::EdgeSet           iterateElements<surface::Edge        >(surface::HalfedgeMesh* mesh);
template<> inline surface::FaceSet           iterateElements<surface::Face        >(surface::HalfedgeMesh* mesh);
template<> inline surface::BoundaryLoopSet   iterateElements<surface::BoundaryLoop>(surface::HalfedgeMesh* mesh);

template<> inline std::list<std::function<void(size_t)>>& getExpandCallbackList<surface::Vertex      >(surface::HalfedgeMesh* mesh);
template<> inline std::list<std::function<void(size_t)>>& getExpandCallbackList<surface::Halfedge    >(surface::HalfedgeMesh* mesh);
template<> inline std::list<std::function<void(size_t)>>& getExpandCallbackList<surface::Corner      >(surface::HalfedgeMesh* mesh);
template<> inline std::list<std::function<void(size_t)>>& getExpandCallbackList<surface::Edge        >(surface::HalfedgeMesh* mesh);
template<> inline std::list<std::function<void(size_t)>>& getExpandCallbackList<surface::Face        >(surface::HalfedgeMesh* mesh);
template<> inline std::list<std::function<void(size_t)>>& getExpandCallbackList<surface::BoundaryLoop>(surface::HalfedgeMesh* mesh);

template<> inline std::list<std::function<void(const std::vector<size_t>&)>>& getPermuteCallbackList<surface::Vertex       >(surface::HalfedgeMesh* mesh);
template<> inline std::list<std::function<void(const std::vector<size_t>&)>>& getPermuteCallbackList<surface::Halfedge     >(surface::HalfedgeMesh* mesh);
template<> inline std::list<std::function<void(const std::vector<size_t>&)>>& getPermuteCallbackList<surface::Corner       >(surface::HalfedgeMesh* mesh);
template<> inline std::list<std::function<void(const std::vector<size_t>&)>>& getPermuteCallbackList<surface::Edge         >(surface::HalfedgeMesh* mesh);
template<> inline std::list<std::function<void(const std::vector<size_t>&)>>& getPermuteCallbackList<surface::Face         >(surface::HalfedgeMesh* mesh);
template<> inline std::list<std::function<void(const std::vector<size_t>&)>>& getPermuteCallbackList<surface::BoundaryLoop >(surface::HalfedgeMesh* mesh);

template<> inline std::string typeShortName<surface::Vertex       >();
template<> inline std::string typeShortName<surface::Halfedge     >();
template<> inline std::string typeShortName<surface::Corner       >();
template<> inline std::string typeShortName<surface::Edge         >();
template<> inline std::string typeShortName<surface::Face         >();
template<> inline std::string typeShortName<surface::BoundaryLoop >();

// clang-format on

} // namespace geometrycentral
