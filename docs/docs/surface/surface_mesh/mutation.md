These routines allow modification of the mesh connectivity and insertion/deletion of elements.

Geometry-central is designed from the ground up to have good support for mesh mutation. The underlying `SurfaceMesh` data structure is index-based, with lazy expansion and deletion, so all operations run in (amortized) constant time with respect to the number of mesh elements, and usually do not incur any memory allocations. [Containers](containers.md) automatically update after mesh operations.

As much as possible, these routines will check for validity before executing and throw an exception if something isn't right. The `NGC_SAFETY_CHECKS` define disables this behavior for a modest increase in performance, but checks are enabled by default even in release builds.

Note that aggressive use of these routines may reduce a mesh from a _simplicial complex_ to a _$\Delta$-complex_. For instance, flipping enough edges in a mesh might create self-edges, which connect a vertex to itself. See the [$\Delta$-complex](delta_complex.md) section for details, and an explanation of why these complexes are important.

!!! note "General `SurfaceMesh`"
    
    Many of these operations could be defined for a general `SurfaceMesh`, but are currently only implemented on the `ManifoldSurfaceMesh`. Future releases of geometry-central will gradually port them over.

### Compressed mode

Internally, the halfedge mesh is represented by dense arrays of indices which are lazily expanded (see [interals](internals.md) for details). To support fast deletion operations, we simply mark elements as deleted, without re-packing the index space. We say that the mesh is _compressed_ if the index space is dense, with no such marked elements. When a mesh is not compressed, the `index` of a mesh element no longer serves as a proper enumeration from `[0,N)`, but merely as a unique ID.

There are two consequences to being non-compressed:

  - Some operations cannot be implemented efficiently/correctly (e.g., random access of the i'th vertex), as noted in the documentation. Calling such function will generally throw an error in Debug mode, but may fail silently in Release.
  - Storage space is wasted by deleted elements.


All meshes are compressed after construction, and only become non-compressed if the user performs an insertion or deletion operation. The `compress()` function can be called to re-index the elements of the mesh as a proper enumeration from `[0,N)`.

The `compress()` function invalidates pointers, and incurs an update of existing containers. As such, it is recommended to be called sporadically, after a sequence of operations is completed.


??? func "`#!cpp bool SurfaceMesh::isCompressed()`"

    Returns true if the mesh is compressed.

??? func "`#!cpp void SurfaceMesh::compress()`"

    Re-index the elements of the mesh to yield a dense enumeration. Invalidates all `Vertex`, `Edge` (etc) objects. All `VertexData<>`, `FaceData<>`, etc containers are automatically resized and re-indexed.

    Does nothing if the mesh is already compressed.


!!! note "Preserving notable elements"

    In some rare situations, you might want to manually keep track of a significant mesh elements (vertices, faces, etc) through a call to `compress()` (which invalidates all element references).

    One way to do this is to leverage the `MeshData<>` containers, which automatically stay valid through updates:
    

    ```cpp
    SurfaceMesh& mesh; // our mesh
    Vertex vA, vB; // two special vertices we want to keep track of

    // Label the special vertices
    VertexData<int> specialVerts(mesh, 0);
    specialVerts[vA] = 1;
    specialVerts[vB] = 2;

    // Compress the mesh
    // (invalidating the Vertex objects)
    mesh.compress();   
    // specialVerts is automatically maintained through the compression

    // Find the interesting vertices by label
    for(Vertex v : mesh.vertices()) {
      if(specialVerts[v] == 1) vA = v;
      if(specialVerts[v] == 2) vB = v;
    }
    ```

    This is an intentionally simplistic example, but generally speaking the `MeshData<>` arrays can be used to track mesh elements and other data through a `compress()`.



## In-place modifications

These routines modify a mesh, but do not require inserting or deleting elements.

??? func "`#!cpp bool SurfaceMesh::flip(Edge e, bool preventSelfEdges = true)`"

    Flip an edge by rotating counter-clockwise.

    An edge cannot be combinatorially flipped if it is:

      - a boundary edge
      - incident on a degree-1 vertex.

    If `true` is passed as the optional argument `preventSelfEdges`, then the edge will also not be flipped if it would result in both endpoints of the edge becoming the same vertex. 

    **Return:** true if the edge was actually flipped 


## Insertions

These routines modify a mesh by inserting new elements. Elements (like `Vertex`) and containers (like `VertexData<>`) will remain valid through insertions and automatically resize themselves to accommodate the new elements. However, the mesh will no longer be [compressed](#compressed-mode).

Note that some operations my re-use existing elements to create their output. For instance, `splitEdge()` turns a single edge in to two; the input edge will be re-used as one of the two output edges, and data along that edge will be unchanged in any containers.

---

??? func "`#!cpp Halfedge ManifoldSurfaceMesh::insertVertexAlongEdge(Edge e)`"

    Adds a degree 2 vertex along an edge. Unlike `splitEdge()`, this _does not_ triangulate the adjacent faces; the degree of adjacent faces will be increased by 1. Works as expected on boundary edges.

    Returns a halfedge `he` along the newly created edge, which points in the same direction as `e.halfedge()`, and such that `he.vertex()` is the newly inserted vertex.

    Preserves canonical direction of edge.halfedge() for both halves of new edge. The original edge is repurposed as one of the two new edges (same for original halfedges).


??? func "`#!cpp Halfedge ManifoldSurfaceMesh::splitEdgeTriangular(Edge e)`"

    Inserts a vertex along an edge, and triangulates the adjacent faces. On a triangle mesh, the newly inserted vertex will be a degree 4 vertex.  Works as expected on boundary edges.

    Returns a halfedge `he` along the newly created edge, which points in the same direction as `e.halfedge()`, and such that `he.vertex()` is the newly inserted vertex.

    Preserves canonical direction of edge.halfedge() for both halves of new edge. The original edge is repurposed as one of the new edges (same for original halfedges).
    

??? func "`#!cpp Vertex ManifoldSurfaceMesh::insertVertex(Face f)`"
    
    Add vertex inside face and triangulate. Returns new vertex.


??? func "`#!cpp Halfedge connectVertices(Halfedge heA, Halfedge heB)`"

    Call to add an edge to a face, splitting it to two faces.

    Creates a new edge connecting `heA.vertex()` to `heB.vertex()`. The initial shared face will be repurposed as one of the two resulting faces.
    
    `heA` and `heB` must be distinct halfedges in the same face, and their vertices must not already be adjacent in that face.

    Returns new halfedge with `heA.vertex()` at tail, and `he.twin().face()` is the new face.


??? func "`#!cpp std::vector<Face> ManifoldSurfaceMesh::triangulate(Face face)`"

    Triangulate a face in the mesh, returning all of the resulting faces.
    
    One of the returned faces will be the input face, repurposed as a face in the triangulation.


## Deletions

These routines delete mesh elements. Elements (like `Vertex`) and containers (like `VertexData<>`) will remain valid through deletions. However, performing any deletion will cause the mesh to no longer be [compressed](#compressed-mode).

??? func "`#!cpp Vertex ManifoldSurfaceMesh::collapseEdgeTriangular(Edge e)`"

    Collapse an edge. Returns the vertex adjacent to that edge which still exists. Returns Vertex() if not
    collapsible.

??? func "`#!cpp bool ManifoldSurfaceMesh::removeFaceAlongBoundary(Face f)`"

    Remove a face which is adjacent to the boundary of the mesh (along with its edge on the boundary).
    Face must have exactly one boundary edge. Returns true if could remove.


