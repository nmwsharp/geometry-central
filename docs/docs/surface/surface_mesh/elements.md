## Introduction

This section covers the elements types for our surface mesh, as well as the traversal and utility functions that they offer.

!!! note

    In the most proper sense, these element types are really "handles" to the underlying element. They refer to a particular element, but the `Vertex` variable in memory is not really the mesh element itself, just a temporary reference to it.
    
    For instance, it is possible (and common) to have multiple `Vertex` variables which actually refer to the same vertex, and allowing a `Vertex` variable to go out of scope certainly does not delete the vertex in the mesh.

    However, the semantics are very natural, so for the sake of brevity we call the type simply `Vertex`, rather than `VertexHandle` (etc).

Additionally, see [navigation](navigation.md) for iterators to travese adjacent elements, like `for(Vertex v : face.adjacentVertices())`.

#### Construction 

Element types do not have constructors which should be called by the user. Instead, the element will always be created for you, via one of several methods, including:

- Iterating through the mesh `for(Vertex v : mesh.vertices())`
- Traversing from a neighbor element `Face f = halfedge.face()`
- Iterating around an element `for(Halfedge he : vertex.outgoingHalfedges())`

Adding a new element to a mesh is covered in the [mutation section](mutation.md).

#### Comparison & Hashing

All mesh elements support:

- **equality checks** (`==`, `!=`)
- **comparions** (`<`, `>`, `<=`, `>=`, according to the iteration order of the elements)
- **hashing** (so they can be used in a `std::unordered_map`)

---

## Vertex

A vertex is a 0-dimensional point which serves as a node in the mesh.

**Traversal:**

??? func "`#!cpp Halfedge Vertex::halfedge()`"

    Returns one of the halfedges whose tail is incident on this vertex. 
    
    If the vertex is a boundary vertex, then it is guaranteed that the returned halfedge will be the unique interior halfedge along the boundary. That is the unique halfedge such that `vertex.halfedge().twin().isInterior() == false`.


??? func "`#!cpp Corner Vertex::corner()`"

    Returns one of the corners incident on the vertex.


**Utility:**

??? func "`#!cpp bool Vertex::isBoundary()`"

    Returns true if the vertex is along the boundary of the mesh. 
    
    See [boundaries](boundaries.md) for more information.


??? func "`#!cpp bool Vertex::isManifold()`"

    Returns true if the vertex neighborhood is manifold.


??? func "`#!cpp bool Vertex::isManifoldAndOriented()`"

    Returns true if the vertex neighborhood is manifold, and all faces have consistent orientations.


??? func "`#!cpp size_t Vertex::degree()`"

    The degree of the vertex, i.e. the number of edges incident on the vertex.


??? func "`#!cpp size_t Vertex::faceDegree()`"

    The face-degree of the vertex, i.e. the number of faces incident on the vertex. On the interior of a mesh, this will be equal to `Vertex::degree()`, at the boundary it will be smaller by one.


---

## Halfedge

A halfedge is a the basic building block of a halfedge mesh. As its name suggests, the halfedge is *half* of an *edge*, connecting two vertices and sitting on on side of an edge in some face. The halfedge is directed, from its _tail_, to its _tip_. Our halfedges have a counter-clockwise orientation: the halfedges with in a face will always point in the counter-clockwise direction. On a `ManifoldSurfaceMesh`, a halfedge and its _twin_ (the neighbor across an edge) will point in opposite directions.

**Traversal:**

??? func "`#!cpp Halfedge Halfedge::twin()`"

    Returns the halfedge's _twin_, a neighbor across an edge. See notes in `sibling()`.
    
    On a `ManifoldSurfaceMesh`, it will point in the opposite direction.  Calling twin twice will always return to the initial halfedge: `halfedge.twin().twin() == halfedge`.


??? func "`#!cpp Halfedge Halfedge::sibling()`"

    A synonym for `twin()`, which makes it clear that on a general `SurfaceMesh` which is not manifold there may be many halfedges incident on an edge. Calling `sibling()` repeatedly orbits around these halfedges.


??? func "`#!cpp Halfedge Halfedge::next()`"

    Returns the _next_ halfedge in the same face as this halfedge, according to the counter-clockwise ordering.


??? func "`#!cpp Vertex Halfedge::vertex()`"

    Returns the vertex at the tail of this halfedge.

??? func "`#!cpp Vertex Halfedge::tailVertex()`"

    Returns the vertex at the tail of this halfedge.

??? func "`#!cpp Vertex Halfedge::tipVertex()`"

    Returns the vertex at the tip of this halfedge.


??? func "`#!cpp Edge Halfedge::edge()`"

    Returns the edge that this halfedge sits along.


??? func "`#!cpp Face Halfedge::face()`"

    Returns the face that this halfedge sits inside.

    Note that in the case of a mesh with boundary, if the halfedge is exterior the result of this function will really be a boundary loop. See [boundaries](boundaries.md) for more information.


??? func "`#!cpp Corner Halfedge::corner()`"

    Returns the corner at the tail of this halfedge.


**Fancy Traversal:**

??? func "`#!cpp Halfedge Halfedge::prevOrbitFace()`"

    Returns the _previous_ halfedge, that is the halfedge such that `he.next() == *this`. This result is found by orbiting around the shared face.
    
    Because our halfedge mesh is singly-connected, this is not a simple $O(1)$ lookup, but must be computed by orbiting around the face. Be careful: calling `he.prevOrbitFace()` on each exterior halfedge can easily lead to $O(N^2)$ algorithm complexity, as each call walks all the way around a a boundary loop.
    
    Generally this operation can (and should) be avoided with proper code structure.


**Utility:**

??? func "`#!cpp bool Halfedge::isInterior()`"

    Returns true if the edge is _interior_, and false if it is _exterior_ (i.e., incident on a boundary loop).
    
    See [boundaries](boundaries.md) for more information.


??? func "`#!cpp bool Halfedge::orientation()`"

    Returns true if the halfedge points in the same direction as the arbitrary orientation defined at each edge.

    For instance, for two halfedges incident on the same edge, you can test if they point in the same direction by checking `heA.orientation() == heB.orientation()`. Recall that two halfedges will point in the same direction along an edge when their respective faces have _opposite_ orientation.

---



## Edge

An _edge_ is a 1-dimensional element that connects two vertices in the mesh.


**Traversal:**

??? func "`#!cpp Halfedge Edge::halfedge()`"

    Returns one of the two halfedges incident on this edge. If the edge is a boundary edge, it is guaranteed that the returned edge will be the interior one.


??? func "`#!cpp Vertex Edge::otherVertex(Vertex v)`"

    Of the two vertices incident on the edge, returns the one which is not `v`. If both incident vertices are `v`, returns `v`.


??? func "`#!cpp Vertex Edge::firstVertex()`"

    Returns one of the two vertices incident on the edge. In particular, returns `edge.halfedge().tailVertex()`.

    See also: `Edge::secondVertex()`.

??? func "`#!cpp Vertex Edge::secondVertex()`"

    Returns one of the two vertices incident on the edge. In particular, returns `edge.halfedge().tipVertex()`.
    
    See also: `Edge::firstVertex()`.
   

**Utility:**

??? func "`#!cpp bool Edge::isBoundary()`"

    Returns true if the edge is along the boundary of the mesh. Note that edges which merely touch the boundary at one endpoint are not considered to be boundary edges.
    
    See [boundaries](boundaries.md) for more information.


??? func "`#!cpp bool Edge::isOriented()`"

    Returns true if the faces on either side have matching orientations. A nonmanifold edge is never oriented.

??? func "`#!cpp size_t Edge::degree()`"

    Returns the number of faces indcident on the edge.

--- 


## Face

A _face_ is a 2-dimensional element formed by a loop of 3 or more edges. In general, our faces can be polygonal with $d \ge 3$ edges, though many of the routines in geometry central are only valid on triangular meshes.

**Traversal:**

??? func "`#!cpp Halfedge Face::halfedge()`"

    Returns any one of the halfedges inside of this face.

??? func "`#!cpp BoundaryLoop Face::asBoundaryLoop()`"

    Reinterprets this element as a boundary loop. Only valid if the face is, in fact, a boundary loop.  See [boundaries](boundaries.md) for more information.

**Utility:**
  
??? func "`#!cpp bool Face::isBoundaryLoop()`"

    Returns true if the face is really a boundary loop.  See [boundaries](boundaries.md) for more information.

??? func "`#!cpp bool Face::isTriangle()`"

    Returns true if the face has three sides.  

??? func "`#!cpp size_t Face::degree()`"

    Returns the number of sides in the face. Complexity $O(d)$, where $d$ is the resulting degree.


--- 


## Boundary Loop

A _boundary loop_ is a special face-like element used to represent holes in the mesh due to surface boundary.  These are only defined on a `ManifoldSurfaceMesh`, which necessarily has coherent boundary loops.  See [boundaries](boundaries.md) for more information.


**Traversal:**

??? func "`#!cpp Halfedge BoundaryLoop::halfedge()`"

    Returns any one of the halfedges inside of the boundary loop.


**Utility:**

??? func "`#!cpp size_t BoundaryLoop::degree()`"

    Returns the number of sides in the boundary loop. Complexity $O(d)$, where $d$ is the resulting degree.


--- 


## Corner

A _corner_ is a convenience type referring to a corner inside of a face. Tracking corners as a separate type is useful, because one often logically represents data defined at corners.


**Traversal:**

??? func "`#!cpp Halfedge Corner::halfedge()`"

    Returns the halfedge whose tail touches this corner. That is to say, `corner.halfedge().vertex() == corner.vertex()`.


??? func "`#!cpp Vertex Corner::vertex()`"

    Returns the vertex which this corner is incident on.

??? func "`#!cpp Face Corner::face()`"

    Returns the face that this corner sits inside of.

