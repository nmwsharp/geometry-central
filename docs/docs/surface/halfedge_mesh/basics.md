# Halfedge meshes

The halfedge mesh is a powerful and flexible data structure for representing oriented, manifold polygonal meshes, and is the core data structure in geometry-central.

The halfedge mesh has several key advantages over other data structures, most notably that all adjacent-neighborhood traversals can be implemented in constant time, without the use of any variably-sized neighbor lists. Furthermore, common mutation operations like edge splits and vertex insertions can be performed in constant time.  This halfedge mesh implementation furthermore stores all elements in contiguous buffers of memory, which makes it fast (see [internals](internals.md) for implementation details).

As the name suggests, the primary type in a halfedge mesh is a _halfedge_, in addition to the usual _vertex_, _edge_ and _face_ types. A halfedge is a directed edge incident on a face, as shown below. Two halfedges, oriented in opposite directions, make up each edge in the mesh. Each halfedge has relationships with five adjacent elements: 

- `Halfedge::twin()` the other halfedge across the incident edge
- `Halfedge::next()` the next halfedge in clockwise order around the incident face
- `Halfedge::vertex()` the vertex at the tail (back) of the halfedge
- `Halfedge::edge()` the incident edge
- `Halfedge::face()` the incident face

![halfedge pointers](../../media/halfedge_pointers.png)

Each vertex, edge, and face need just one relationship:

- `Vertex::halfedge()` _any_ of the incident halfedges (which point outward from the vertex)
- `Edge::halfedge()` _any_ of the incident halfedges
- `Face::halfedge()` _any_ of the incident halfedges

In fact, this fixed set of relationships is sufficient to implement pretty much _any_ local traversal. Geometry central provides a wide range of convience iterators which wrap these relationships to traverse neighborhoods, such as the example below.
```cpp
for(Edge e : vertex.adjacentEdges()) {
  // do science
}
```
See [navigation](navigation.md) for more information on traversals and convenience iterators.

Notice that the lightweight `Halfedge` (etc) types serve simply as logical references, or "handles" to a mesh element. Deleting one of these handles does not delete the underlying element, and one may have multiple handles to the same element `Vertex a; Vertex b; a == b;`.

## Manifold, Oriented Surfaces

The basic halfedge mesh imposes two requirements: manifoldness and orientability. 

Manifoldness means that our surface must locally look like a plane in any neighborhood. This disallows structures such as three faces meeting at an edge, or two cones of faces meeting at a single vertex like an hourglass. 

Furthermore the halfedge mesh implies a combinatorial _orientation_ of the surface, indicated by the clockwise ordering of halfedges around each face (see figure below). Because the halfedge mesh implies an orientation, it cannot represent non-orientable surfaces, like a Klein bottle.

![halfedge orientation](../../media/halfedge_orientation.png)

These properties are invariants which always hold for any meaningful halfedge mesh; in practice we check them during construction and ensure that all operations preserve them.

Note that our halfedge mesh _does not_ require that faces be triangles or quads; arbitrary faces with degree >= 3 are supported, and faces of different degree may be intermingled. However, many operations are only defined for triangle meshes and will throw errors if invoked on other meshes.


## Basic API


### Constructors


??? func "`#!cpp HalfedgeMesh(const std::vector<std::vector<size_t>>& polygons, bool verbose = false)`"
    Constructs a halfedge mesh from a face-index list.

    - `polygons` a list of faces, each holding the indices of the vertices incident on that face, zero-indexed and in counter-clockwise order.
    - `verbose` if true, prints some statistics to `std::cout` during construction.

### Element counts

??? func "`#!cpp size_t HalfedgeMesh::nVertices()`"
    Returns the number of vertices. 

??? func "`#!cpp size_t HalfedgeMesh::nInteriorVertices()`"
    Returns the number of vertices not incident on the boundary.

??? func "`#!cpp size_t HalfedgeMesh::nBoundaryVertices()`"
    Returns the number of vertices incident on the boundary.

??? func "`#!cpp size_t HalfedgeMesh::nEdges()`"
    Returns the number of edges. 

??? func "`#!cpp size_t HalfedgeMesh::nFaces()`"
    Returns the number of faces in the mesh.

??? func "`#!cpp size_t HalfedgeMesh::nHalfedges()`"
    Returns the number of halfedges, including both interior halfedges and any exterior halfedges incident on boundary loops. Always exactly twice the number of edges.

??? func "`#!cpp size_t HalfedgeMesh::nInterioHalfedges()`"
    Returns the number of interior halfedges, which are incident on faces of the mesh. Always equal to the sum of the number of sides of all faces.

??? func "`#!cpp size_t HalfedgeMesh::nExteriorHalfedges()`"
    Returns the number of exterior halfedges, which are opposite boundary faces. 

??? func "`#!cpp size_t HalfedgeMesh::nBoundaryLoops()`"
    Returns the number of distinct boundary loops in the mesh, each identified as an fictional face closing a boundary loop in the mesh.


### Properties

??? func "`#!cpp bool HalfedgeMesh::hasBoundary()`"
    Returns true if the mesh has boundary, that is if it is not _closed_.
    
    Complexity $\mathcal{O}(1)$.

??? func "`#!cpp int HalfedgeMesh::eulerCharacteristic()`"
    Returns the Euler characteristic of the surface. Computed in O(1) from element counts. 
    
    **Note:** always computed by naively applying [Euler's polyhedron formula](https://en.wikipedia.org/wiki/Euler_characteristic#Polyhedra), which might not do what you want in the case of multiple-connected components.

??? func "`#!cpp int HalfedgeMesh::genus()`"
    Returns the genus of the surface. Computed in O(1) from element counts.
    
    **Note:** always computed by naively applying [Euler's polyhedron formula](https://en.wikipedia.org/wiki/Euler_characteristic#Polyhedra), which might not do what you want in the case of multiple connected components.

??? func "`#!cpp bool HalfedgeMesh::isTriangular()`"
    Returns true if all faces in the mesh have 3 sides. 
    
    Complexity $\mathcal{O}(n)$, do not call in a tight loop.

??? func "`#!cpp size_t HalfedgeMesh::nConnectedComponents()`"
    Returns the number of distinct connected components of the mesh. 
    
    Complexity $\mathcal{O}(n)$, do not call in a tight loop.


### Utility functions   

??? func "`#!cpp std::unique_ptr<Halfedgemesh> copy()`"
    Constructs a copy of the mesh.


