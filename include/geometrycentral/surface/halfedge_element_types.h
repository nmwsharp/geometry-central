#pragma once

#include "geometrycentral/utilities/element.h"
#include "geometrycentral/utilities/element_iterators.h"
#include "geometrycentral/utilities/mesh_data.h"
#include "geometrycentral/utilities/utilities.h"

#include <cstddef>
#include <iostream>
#include <list>
#include <typeindex>
#include <array>
#include <unordered_set>

namespace geometrycentral {
namespace surface {

// === Types and inline methods for the halfedge mesh pointer and datatypes
class SurfaceMesh;

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
struct EdgeAdjacentHalfedgeNavigator;
struct EdgeAdjacentInteriorHalfedgeNavigator;
struct EdgeAdjacentFaceNavigator;
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

class Vertex : public Element<Vertex, SurfaceMesh> {
public:
  // Constructors
  // inheriting constructor would work here, and replace the constructors below, but gcc-5 erroneously rejects the combo
  // with CRTP :( perhaps resurrect here and in other elements below once gcc-5 is sufficiently old
  // using Element<Vertex>::Element;
  Vertex();                              // construct an empty (null) element
  Vertex(SurfaceMesh* mesh, size_t ind); // construct pointing to the i'th element of that type on a mesh.
  // Vertex(const DynamicElement<Vertex>& e); // construct from a dynamic element of matching type

  // Navigators
  Halfedge halfedge() const;
  Corner corner() const;
  bool isDead() const;

  // Properties
  bool isBoundary() const;
  bool isManifold() const;
  bool isManifoldAndOriented() const; 
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

// using DynamicVertex = DynamicElement<Vertex>;

// == Range iterators

// All vertices
struct VertexRangeF {
  static bool elementOkay(const SurfaceMesh& mesh, size_t ind);
  typedef Vertex Etype;
  typedef SurfaceMesh ParentMeshT;
};
typedef RangeSetBase<VertexRangeF> VertexSet;


// ==========================================================
// ================        Halfedge        ==================
// ==========================================================

class Halfedge : public Element<Halfedge, SurfaceMesh> {
public:
  // Constructors
  Halfedge();                              // construct an empty (null) element
  Halfedge(SurfaceMesh* mesh, size_t ind); // construct pointing to the i'th element of that type on a mesh.
  // Halfedge(const DynamicElement<Halfedge>& e); // construct from a dynamic element of matching type

  // Navigators
  Halfedge twin() const;
  Halfedge sibling() const;
  Halfedge nextOutgoingNeighbor() const; // next halfedge which has the same tail vertex as this, form a cycle
  Halfedge nextIncomingNeighbor() const; // next halfedge which has the same tip vertex as this, form a cycle
  Halfedge next() const;
  Corner corner() const;
  Vertex vertex() const;
  Vertex tipVertex() const;
  Vertex tailVertex() const;
  Edge edge() const;
  Face face() const;
  bool isDead() const;

  // Super-navigators
  Halfedge prevOrbitFace() const;
  Halfedge prevOrbitVertex() const; // only meaningful if manifold

  // Properties
  bool isInterior() const;
  bool orientation() const; // True if the halfedge has the same orientation as its edge
                            // Remember that halfedge orientation means "points in direction". If two faces have the
                            // same orientation, their halfedges along a shared edge will have opposite orientations.
};

// using DynamicHalfedge = DynamicElement<Halfedge>;

// == Range iterators

// All halfedges
struct HalfedgeRangeF {
  static bool elementOkay(const SurfaceMesh& mesh, size_t ind);
  typedef Halfedge Etype;
  typedef SurfaceMesh ParentMeshT;
};
typedef RangeSetBase<HalfedgeRangeF> HalfedgeSet;

// Interior halfedges
struct HalfedgeInteriorRangeF {
  static bool elementOkay(const SurfaceMesh& mesh, size_t ind);
  typedef Halfedge Etype;
  typedef SurfaceMesh ParentMeshT;
};
typedef RangeSetBase<HalfedgeInteriorRangeF> HalfedgeInteriorSet;

// Exterior halfedges
struct HalfedgeExteriorRangeF {
  static bool elementOkay(const SurfaceMesh& mesh, size_t ind);
  typedef Halfedge Etype;
  typedef SurfaceMesh ParentMeshT;
};
typedef RangeSetBase<HalfedgeExteriorRangeF> HalfedgeExteriorSet;


// ==========================================================
// ================        Corner          ==================
// ==========================================================

// Implmentation note: The `ind` parameter for a corner will be the index of a halfedge, which should always be
// interior.

class Corner : public Element<Corner, SurfaceMesh> {
public:
  // Constructors
  Corner();                              // construct an empty (null) element
  Corner(SurfaceMesh* mesh, size_t ind); // construct pointing to the i'th element of that type on a mesh.
  // Corner(const DynamicElement<Corner>& e); // construct from a dynamic element of matching type

  // Navigators
  Halfedge halfedge() const;
  Vertex vertex() const;
  Face face() const;
  bool isDead() const;
};

// using DynamicCorner = DynamicElement<Corner>;

// == Range iterators

// All corners
struct CornerRangeF {
  static bool elementOkay(const SurfaceMesh& mesh, size_t ind);
  typedef Corner Etype;
  typedef SurfaceMesh ParentMeshT;
};
typedef RangeSetBase<CornerRangeF> CornerSet;


// ==========================================================
// ================          Edge          ==================
// ==========================================================

class Edge : public Element<Edge, SurfaceMesh> {
public:
  // Constructors
  Edge();                              // construct an empty (null) element
  Edge(SurfaceMesh* mesh, size_t ind); // construct pointing to the i'th element of that type on a mesh.
  // Edge(const DynamicElement<Edge>& e); // construct from a dynamic element of matching type

  // Navigators
  Halfedge halfedge() const;
  Vertex otherVertex(Vertex v) const;
  Vertex firstVertex() const;
  Vertex secondVertex() const;
  std::array<Halfedge,4> diamondBoundary() const;
  bool isDead() const;

  // Properties
  bool isBoundary() const;
  bool isManifold() const;
  bool isOriented() const;
  size_t degree() const;

  // Iterators
  NavigationSetBase<EdgeAdjacentHalfedgeNavigator> adjacentHalfedges() const;
  NavigationSetBase<EdgeAdjacentInteriorHalfedgeNavigator> adjacentInteriorHalfedges() const;
  NavigationSetBase<EdgeAdjacentFaceNavigator> adjacentFaces() const;
  std::array<Vertex, 2> adjacentVertices() const;
};

// using DynamicEdge = DynamicElement<Edge>;

// == Range iterators

// All edges
struct EdgeRangeF {
  static bool elementOkay(const SurfaceMesh& mesh, size_t ind);
  typedef Edge Etype;
  typedef SurfaceMesh ParentMeshT;
};
typedef RangeSetBase<EdgeRangeF> EdgeSet;


// ==========================================================
// ================          Face          ==================
// ==========================================================

// Implmentation note: The `ind` parameter for a face might correspond to a boundary loop. The boundary loops have face
// IDs which are at the very end of the face buffer, but can still index in to face-valued arrays/functions in
// SurfaceMesh (they _cannot_ index in to FaceData<> containers).

class Face : public Element<Face, SurfaceMesh> {
public:
  // Constructors
  Face();                              // construct an empty (null) element
  Face(SurfaceMesh* mesh, size_t ind); // construct pointing to the i'th element of that type on a mesh.
  // Face(const DynamicElement<Face>& e); // construct from a dynamic element of matching type

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

// using DynamicFace = DynamicElement<Face>;

// == Range iterators

// All faces
struct FaceRangeF {
  static bool elementOkay(const SurfaceMesh& mesh, size_t ind);
  typedef Face Etype;
  typedef SurfaceMesh ParentMeshT;
};
typedef RangeSetBase<FaceRangeF> FaceSet;


// ==========================================================
// ================     Boundary Loop      ==================
// ==========================================================

// Implementation note: the `ind` parameter for a boundary loop is index from the back of the face index space, from [0,
// nBoundaryLoopFillCount).

class BoundaryLoop : public Element<BoundaryLoop, SurfaceMesh> {
public:
  // Constructors
  BoundaryLoop();                              // construct an empty (null) element
  BoundaryLoop(SurfaceMesh* mesh, size_t ind); // construct pointing to the i'th element of that type on a mesh.
  // BoundaryLoop(const DynamicElement<BoundaryLoop>& e); // construct from a dynamic element of matching type

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

// using DynamicBoundaryLoop = DynamicElement<BoundaryLoop>;

// == Range iterators

// All boundary loops
struct BoundaryLoopRangeF {
  static bool elementOkay(const SurfaceMesh& mesh, size_t ind);
  typedef BoundaryLoop Etype;
  typedef SurfaceMesh ParentMeshT;
};
typedef RangeSetBase<BoundaryLoopRangeF> BoundaryLoopSet;


// ==========================================================
// ===============   Navigation Iterators   =================
// ==========================================================

// == Vertex

// Helper class to store some special extra state for the vertex iterators
// For implicit-twin meshes, it just does the usual twin.next() and currHe is always an outgoing halfedge. For general
// meshes, it iterates first through the outgoing, then the incoming halfedges (always stored in currHe).
struct VertexNeighborIteratorState {
  VertexNeighborIteratorState(Halfedge currHeOutgoing, bool useImplicitTwin);

  const bool useImplicitTwin;
  Halfedge currHe = Halfedge();

  // if useImplicitTwin == false, this is populated
  bool processingIncoming = false;
  Halfedge firstHe = Halfedge();

  void advance();
  bool isHalfedgeCanonical() const; // this currently pointing at the one canonical halfedge along an edge
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
  typedef Halfedge Etype;
  Etype currE;
  typedef Halfedge Rtype;
  Rtype getCurrent() const;
};

// Adjacent outgoing halfedges
struct VertexOutgoingHalfedgeNavigator {
  void advance();
  bool isValid() const;
  typedef Halfedge Etype;
  Etype currE;
  typedef Halfedge Rtype;
  Rtype getCurrent() const;
};

// Adjacent corners
struct VertexAdjacentCornerNavigator {
  void advance();
  bool isValid() const;
  typedef Halfedge Etype;
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
  typedef Halfedge Etype;
  Etype currE;
  typedef Face Rtype;
  Rtype getCurrent() const;
};

// == Edge


// Adjacent halfedge
struct EdgeAdjacentHalfedgeNavigator {
  void advance();
  bool isValid() const;
  typedef Halfedge Etype;
  Etype currE;
  typedef Halfedge Rtype;
  Rtype getCurrent() const;
};

// Adjacent interior halfedge
struct EdgeAdjacentInteriorHalfedgeNavigator {
  void advance();
  bool isValid() const;
  typedef Halfedge Etype;
  Etype currE;
  typedef Halfedge Rtype;
  Rtype getCurrent() const;
};

// Adjacent interior halfedge
struct EdgeAdjacentFaceNavigator {
  void advance();
  bool isValid() const;
  typedef Halfedge Etype;
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
  typedef std::pair<Halfedge, Halfedge>
      Etype; // first is current halfedge of this face, second halfedge adjacent to face to return
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

template<> inline size_t nElements<surface::Vertex       >(surface::SurfaceMesh* mesh); 
template<> inline size_t nElements<surface::Face         >(surface::SurfaceMesh* mesh); 
template<> inline size_t nElements<surface::Edge         >(surface::SurfaceMesh* mesh); 
template<> inline size_t nElements<surface::Halfedge     >(surface::SurfaceMesh* mesh); 
template<> inline size_t nElements<surface::Corner       >(surface::SurfaceMesh* mesh); 
template<> inline size_t nElements<surface::BoundaryLoop >(surface::SurfaceMesh* mesh); 

template<> inline size_t elementCapacity<surface::Vertex      >(surface::SurfaceMesh* mesh);
template<> inline size_t elementCapacity<surface::Face        >(surface::SurfaceMesh* mesh);
template<> inline size_t elementCapacity<surface::Edge        >(surface::SurfaceMesh* mesh);
template<> inline size_t elementCapacity<surface::Halfedge    >(surface::SurfaceMesh* mesh);
template<> inline size_t elementCapacity<surface::Corner      >(surface::SurfaceMesh* mesh);
template<> inline size_t elementCapacity<surface::BoundaryLoop>(surface::SurfaceMesh* mesh);

template<> inline size_t dataIndexOfElement<surface::Vertex          >(surface::SurfaceMesh* mesh, surface::Vertex e           );
template<> inline size_t dataIndexOfElement<surface::Face            >(surface::SurfaceMesh* mesh, surface::Face e             );
template<> inline size_t dataIndexOfElement<surface::Edge            >(surface::SurfaceMesh* mesh, surface::Edge e             );
template<> inline size_t dataIndexOfElement<surface::Halfedge        >(surface::SurfaceMesh* mesh, surface::Halfedge e         );
template<> inline size_t dataIndexOfElement<surface::Corner          >(surface::SurfaceMesh* mesh, surface::Corner e           );
template<> inline size_t dataIndexOfElement<surface::BoundaryLoop    >(surface::SurfaceMesh* mesh, surface::BoundaryLoop e     );

template<> struct ElementSetType<surface::Vertex        >   { typedef surface::VertexSet       type; };
template<> struct ElementSetType<surface::Face          >   { typedef surface::FaceSet         type; };
template<> struct ElementSetType<surface::Edge          >   { typedef surface::EdgeSet         type; };
template<> struct ElementSetType<surface::Halfedge      >   { typedef surface::HalfedgeSet     type; };
template<> struct ElementSetType<surface::Corner        >   { typedef surface::CornerSet       type; };
template<> struct ElementSetType<surface::BoundaryLoop  >   { typedef surface::BoundaryLoopSet type; };

template<> inline surface::VertexSet         iterateElements<surface::Vertex      >(surface::SurfaceMesh* mesh);
template<> inline surface::HalfedgeSet       iterateElements<surface::Halfedge    >(surface::SurfaceMesh* mesh);
template<> inline surface::CornerSet         iterateElements<surface::Corner      >(surface::SurfaceMesh* mesh);
template<> inline surface::EdgeSet           iterateElements<surface::Edge        >(surface::SurfaceMesh* mesh);
template<> inline surface::FaceSet           iterateElements<surface::Face        >(surface::SurfaceMesh* mesh);
template<> inline surface::BoundaryLoopSet   iterateElements<surface::BoundaryLoop>(surface::SurfaceMesh* mesh);

template<> inline std::list<std::function<void(size_t)>>& getExpandCallbackList<surface::Vertex      >(surface::SurfaceMesh* mesh);
template<> inline std::list<std::function<void(size_t)>>& getExpandCallbackList<surface::Halfedge    >(surface::SurfaceMesh* mesh);
template<> inline std::list<std::function<void(size_t)>>& getExpandCallbackList<surface::Corner      >(surface::SurfaceMesh* mesh);
template<> inline std::list<std::function<void(size_t)>>& getExpandCallbackList<surface::Edge        >(surface::SurfaceMesh* mesh);
template<> inline std::list<std::function<void(size_t)>>& getExpandCallbackList<surface::Face        >(surface::SurfaceMesh* mesh);
template<> inline std::list<std::function<void(size_t)>>& getExpandCallbackList<surface::BoundaryLoop>(surface::SurfaceMesh* mesh);

template<> inline std::list<std::function<void(const std::vector<size_t>&)>>& getPermuteCallbackList<surface::Vertex       >(surface::SurfaceMesh* mesh);
template<> inline std::list<std::function<void(const std::vector<size_t>&)>>& getPermuteCallbackList<surface::Halfedge     >(surface::SurfaceMesh* mesh);
template<> inline std::list<std::function<void(const std::vector<size_t>&)>>& getPermuteCallbackList<surface::Corner       >(surface::SurfaceMesh* mesh);
template<> inline std::list<std::function<void(const std::vector<size_t>&)>>& getPermuteCallbackList<surface::Edge         >(surface::SurfaceMesh* mesh);
template<> inline std::list<std::function<void(const std::vector<size_t>&)>>& getPermuteCallbackList<surface::Face         >(surface::SurfaceMesh* mesh);
template<> inline std::list<std::function<void(const std::vector<size_t>&)>>& getPermuteCallbackList<surface::BoundaryLoop >(surface::SurfaceMesh* mesh);

template<> inline std::string typeShortName<surface::Vertex       >();
template<> inline std::string typeShortName<surface::Halfedge     >();
template<> inline std::string typeShortName<surface::Corner       >();
template<> inline std::string typeShortName<surface::Edge         >();
template<> inline std::string typeShortName<surface::Face         >();
template<> inline std::string typeShortName<surface::BoundaryLoop >();

// clang-format on

} // namespace geometrycentral

