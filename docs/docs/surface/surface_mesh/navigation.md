### Collection Iterators 

Use these routines to iterate over all of the elements in the mesh.

**Note:** Generally, modifying the mesh while iterating is allowed, but the new elements may or may not be iterated over, and previous elements might even appear again later later in the iteration after modifying.

---

??? func "`#!cpp SurfaceMesh::vertices()`"
    Iterate over the vertices in a mesh.
    ```cpp
    for(Vertex v : mesh.vertices()) {
      // do science here
    }
    ```

??? func "`#!cpp SurfaceMesh::halfedges()`"
    Iterate over all of the halfedges in a mesh (both interior and exterior, if the mesh has boundary).
    ```cpp
    for(Halfedge he : mesh.halfedges()) {
      // do science here
    }
    ```

??? func "`#!cpp SurfaceMesh::interiorHalfedges()`"
    Iterate over the interior halfedges in a mesh.
    ```cpp
    for(Halfedge he : mesh.interiorHalfedges()) {
      // do science here
    }
    ```
    Note that on a boundary edge between vertices `i <--> j`, this set will only include a halfedge from `i --> j`, but not from `j --> i` (or vice versa).

??? func "`#!cpp SurfaceMesh::exteriorHalfedges()`"
    Iterate over the exterior halfedges in a mesh. (only useful on `ManifoldSurfaceMesh`)
    ```cpp
    for(Halfedge he : mesh.exteriorHalfedges()) {
      // do science here
    }
    ```
    Note that on a boundary edge between vertices `i <--> j`, this set will only include a halfedge from `i --> j`, but not from `j --> i` (or vice versa).

??? func "`#!cpp SurfaceMesh::edges()`"
    Iterate over the edges in a mesh.
    ```cpp
    for(Edge e : mesh.edges()) {
      // do science here
    }
    ```

??? func "`#!cpp SurfaceMesh::faces()`"
    Iterate over the faces in a mesh.
    ```cpp
    for(Face f : mesh.faces()) {
      // do science here
    }
    ```

??? func "`#!cpp SurfaceMesh::boundaryLoops()`"
    Iterate over the boundary loops for a mesh.

    Remember that only `ManifoldSurfaceMesh`s have well-defined boundary loops.

    ```cpp
    for(BoundaryLoop bl : mesh.boundaryLoops()) {
      // do science here
    }
    ```


## Neighborhood Iterators 

Use these routines to iterate over the neighbors of a mesh element.


??? note "Note: neighborhoods on $\Delta$-complexes"
    The iterators in this section may have unexpected behavior in the advanced case of a $\Delta$-complex, when there are (e.g.) self-edges, or multiple edges between a pair of vertices. Essentially, these iterators always naively traverse the local neighborhood, even if that neighborhood might include duplicate elements. 
    
    For instance, if a $\Delta$-complex has multiple edges connecting vertex `va` to vertex `vb`, then iterating `va.adjacentVertices()` will return `vb` multiple times.
    
    Of course, for ordinary triangle mesh they will behave as expected. See the [Delta complex](delta_complex.md) section for more information.

---

### Around a vertex

??? func "`#!cpp Vertex::outgoingHalfedges()`"
    Iterate over the halfedges which point outward from a vertex.
    ```cpp
    for(Halfedge he : vert.outgoingHalfedges()) {
      assert(he.vertex() == vert); // true
      // do science here
    }
    ```

??? func "`#!cpp Vertex::incomingHalfedges()`"
    Iterate over the halfedges which point inward at a vertex.
    ```cpp
    for(Halfedge he : vert.incomingHalfedges()) {
      assert(he.twin().vertex() == vert); // true
      // do science here
    }
    ```

??? func "`#!cpp Vertex::adjacentVertices()`"

    Iterate over the vertices edge-connected to this vertex.
    ```cpp
    for(Vertex v : vert.adjacentVertices()) {
      // do science here
    }
    ```

??? func "`#!cpp Vertex::adjacentEdges()`"

    Iterate over the edges incident on this vertex.
    ```cpp
    for(Edge e : vert.adjacentEdges()) {
      // do science here
    }
    ```

??? func "`#!cpp Vertex::adjacentFaces()`"

    Iterate over the faces incident on this vertex.
    ```cpp
    for(Face f : vert.adjacentFaces()) {
      // do science here
    }
    ```

### Around an edge

??? func "`#!cpp Edge::adjacentHalfedges()`"

    Iterate over the halfedges incident on this edge.
    ```cpp
    for(Halfedge he : edge.adjacentHalfedges()) {
      // do science here
    }
    ```

??? func "`#!cpp Edge::adjacentFaces()`"

    Iterate over the (one or two) faces incident on this edge.
    ```cpp
    for(Face f : edge.adjacentFaces()) {
      // do science here
    }
    ```


??? func "`#!cpp Edge::adjacentVertices()`"

    Iterate over the (two) vertices which are endpoints of the edge.
    ```cpp
    for (Vertex v : edge.adjacentVertices()) {
      // do science here
    }
    ```

    Note: unlike most navigators, this routine actually returns a fixed-size array, so you can alternately write things like:
    ```cpp 
    std::array<Vertex, 2> verts = edge.adjacentVertices();
    ```


??? func "`#!cpp Edge::diamondBoundary()`"

    Iterate over the four halfedges bounding the diamond with this edge as its center diagonal.

    More precisely, for an interior edge on a manifold triangle mesh, this returns
    ```cpp
    Halfedge he = edge.halfedge();
    return {he.next(), he.next().next(), he.twin().next(), he.twin().next().next()}
    ```


    Example:
    ```cpp
    for (Halfedge he : edge.diamondBoundary()) {
      // do science here
    }
    ```

    Throws an exception if there are non-triangular faces, the edge is on the boundary, or if the edge is nonmanifold.

    Note: unlike most navigators, this routine actually returns a fixed-size array, so you can alternately write things like:
    ```cpp 
    std::array<Halfedge, 4> halfedges = edge.diamondBoundary();
    ```


### Around a face

??? func "`#!cpp Face::adjacentVertices()`"
    Iterate over the vertices adjacent to a face.
    ```cpp
    for(Vertex v : face.adjacentVertices()) {
      // do science here
    }
    ```

??? func "`#!cpp Face::adjacentHalfedges()`"
    Iterate over the halfedges incident on a face.
    ```cpp
    for(Halfedge he : face.adjacentHalfedges()) {
      // do science here
    }
    ```

??? func "`#!cpp Face::adjacentEdges()`"
    Iterate over the edges on the boundary of a face.
    ```cpp
    for(Edge e : face.adjacentEdges()) {
      // do science here
    }
    ```

??? func "`#!cpp Face::adjacentFaces()`"
    Iterate over the faces adjacent to a face, across each edge.
    ```cpp
    for(Face f : face.adjacentFaces()) {
      // do science here
    }
    ```


### Around a boundary loop

??? func "`#!cpp BoundaryLoop::adjacentVertices()`"
    Iterate over the vertices adjacent to a boundary loop.
    ```cpp
    for(Vertex v : boundaryLoop.adjacentVertices()) {
      // do science here
    }
    ```

??? func "`#!cpp BoundaryLoop::adjacentHalfedges()`"
    Iterate over the (exterior) halfedges incident on a boundary loop.
    ```cpp
    for(Halfedge he : boundaryLoop.adjacentHalfedges()) {
      // do science here
    }
    ```

??? func "`#!cpp BoundaryLoop::adjacentEdges()`"
    Iterate over the edges on the boundary of a boundary loop.
    ```cpp
    for(Edge e : boundaryLoop.adjacentEdges()) {
      // do science here
    }
    ```


## Accessors 

Use these routines to access elements of the mesh by their index.

!!! warning
    The indexing routines in the section are only valid when the mesh is [compressed](mutation.md#compressed-mode).

---

??? func "`#!cpp Halfedge SurfaceMesh::halfedge(size_t index)`"
    Constructs a reference to the i'th halfedge in the mesh. `0 <= index < nHalfedges()`.
    
??? func "`#!cpp Vertex SurfaceMesh::vertex(size_t index)`"
    Constructs a reference to the i'th vertex in the mesh. `0 <= index < nVertices()`.
    
??? func "`#!cpp Face SurfaceMesh::face(size_t index)`"
    Constructs a reference to the i'th face in the mesh. `0 <= index < nFaces()`.
    
??? func "`#!cpp Edge SurfaceMesh::edge(size_t index)`"
    Constructs a reference to the i'th edge in the mesh. `0 <= index < nEdges()`.
    
??? func "`#!cpp Face SurfaceMesh::face(size_t index)`"
    Constructs a reference to the i'th face in the mesh. `0 <= index < nFaces()`.

??? func "`#!cpp Face SurfaceMesh::boundaryLoop(size_t index)`"
    Constructs a reference to the i'th boundary loop in the mesh. `0 <= index < nBoundaryLoops()`.
