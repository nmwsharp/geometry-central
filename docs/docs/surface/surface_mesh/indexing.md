Several methods can be used to obtain integer indices associated with mesh elements.  Such indices are often needed for, e.g., building matrices.  These methods index elements from 0 through the number of elements of the same kind.  I.e., vertices are indexed 0, 1, 2, ..., up through the number of vertices.  In general nothing can/should be assumed about this ordering, except that it will remain fixed as long as the mesh connectivity is not changed.

Note that in some cases, the most straightforward thing may be to just assign your own indices---especially if the indexing scheme needed to build your matrix is different from a straightforward enumeration of mesh elements.  For instance, in the following example a boolean vertex attribute `skip` is provided, which determines whether each vertex is used in a linear system:

```cpp
VertexData<bool> skip; // assume skip[v] is true if v is not used in the linear system
VertexData<size_t> vertexIndex( mesh );
int i = 0;
for( Vertex v : mesh.vertices() ) {
   if( !skip[v] ) {
      vertexIndex[v] = i;
      i++;
   }
}
```

## Uncached index arrays

These methods build and return an `ElementData<size_t>` array that indexes the mesh elements. Note that these arrays are not cached, and arrays are re-built each time the method is invoked.

??? func "`#!cpp HalfedgeData<size_t> SurfaceMesh::getHalfedgeIndices()`"

    Returns an array that provides indices for all halfedges, with values between 0 and _H_-1, where _H_ is the number of halfedges in the mesh.

??? func "`#!cpp CornerData<size_t> SurfaceMesh::getCornerIndices()`"

    Returns an array that provides indices for all corners, with values between 0 and _C_-1, where _C_ is the number of corners in the mesh.

??? func "`#!cpp VertexData<size_t> SurfaceMesh::getVertexIndices()`"

    Returns an array that provides indices for all vertices, with values between 0 and _V_-1, where _V_ is the number of vertices in the mesh.

??? func "`#!cpp VertexData<size_t> SurfaceMesh::getInteriorVertexIndices()`"

    Returns an array that provides indices for all interior vertices, with values between 0 and _I_-1, where _I_ is the number of interior vertices in the mesh.

??? func "`#!cpp EdgeData<size_t> SurfaceMesh::getEdgeIndices()`"

    Returns an array that provides indices for all edges, with values between 0 and _E_-1, where _E_ is the number of edges in the mesh.

??? func "`#!cpp FaceData<size_t> SurfaceMesh::getFaceIndices()`"

    Returns an array that provides indices for all faces, with values between 0 and _F_-1, where _F_ is the number of faces in the mesh.

??? func "`#!cpp BoundaryLoopData<size_t> SurfaceMesh::getBoundaryLoopIndices()`"

    Returns an array that provides indices for all boundary loops, with values between 0 and _B_-1, where _B_ is the number of boundary loops in the mesh.

## Cached index arrays

These methods build the same arrays as the uncached methods, but cache them in a `Geometry` object. Note that indices of course do not relate to the actual geometry (i.e., shape and size) of the mesh---however, this mechanism allows indices to be tracked with the geometry classes, so they can participate in the buffer management and dependency management.  Basic example usage:

```cpp
SurfaceMesh mesh;
VertexPositionGeometry geometry( mesh );
geometry->requireVertexIndices();
for( Vertex v : mesh->vertices() ) {
   int i = geometry->vertexIndices[v];
   // do science
}
```

??? func "`#!cpp void BaseGeometryInterface::requireHalfedgeIndices()`"

    Builds and maintains the array `BaseGeometryInterface::halfedgeIndices` that provides indices for all halfedges, with values between 0 and _H_-1, where _H_ is the number of halfedges in the mesh.

??? func "`#!cpp void BaseGeometryInterface::requireCornerIndices()`"

    Builds and maintains the array `BaseGeometryInterface::cornerIndices` that provides indices for all corners, with values between 0 and _C_-1, where _C_ is the number of halfedges in the mesh.

??? func "`#!cpp void BaseGeometryInterface::requireVertexIndices()`"

    Builds and maintains the array `BaseGeometryInterface::vertexIndices` that provides indices for all vertices, with values between 0 and _V_-1, where _V_ is the number of vertices in the mesh.

??? func "`#!cpp void BaseGeometryInterface::requireInteriorVertexIndices()`"

    Builds and maintains the array `BaseGeometryInterface::interiorVertexIndices` that provides indices for all interior vertices, with values between 0 and _I_-1, where _I_ is the number of interior vertices in the mesh.

??? func "`#!cpp void BaseGeometryInterface::requireEdgeIndices()`"

    Builds and maintains the array `BaseGeometryInterface::edgeIndices` that provides indices for all edges, with values between 0 and _E_-1, where _E_ is the number of edges in the mesh.

??? func "`#!cpp void BaseGeometryInterface::requireFaceIndices()`"

    Builds and maintains the array `BaseGeometryInterface::faceIndices` that provides indices for all faces, with values between 0 and _F_-1, where _F_ is the number of faces in the mesh.

??? func "`#!cpp void BaseGeometryInterface::requireBoundaryLoopIndices()`"

    Builds and maintains the array `BaseGeometryInterface::boundaryLoopIndices` that provides indices for all boundary loops, with values between 0 and _B_-1, where _B_ is the number of boundary loops in the mesh.

## Per element indices

!!! warning
    The routines in this section are valid only when the mesh is [compressed](/surface/surface_mesh/mutation/#compressed-mode).

These methods access the index directly from the mesh element itself.  For example:

```cpp
Edge e;
size_t i = e.getIndex();
```

??? func "`#!cpp size_t Halfedge::getIndex() const`"

    Returns an index between 0 and _H_-1, where _H_ is the number of halfedges in the mesh.

??? func "`#!cpp size_t Corner::getIndex() const`"

    Returns an index between 0 and _C_-1, where _C_ is the number of corners in the mesh.

??? func "`#!cpp size_t Vertex::getIndex() const`"

    Returns an index between 0 and _V_-1, where _V_ is the number of vertices in the mesh.

??? func "`#!cpp size_t Edge::getIndex() const`"

    Returns an index between 0 and _E_-1, where _E_ is the number of edges in the mesh.

??? func "`#!cpp size_t Face::getIndex() const`"

    Returns an index between 0 and _F_-1, where _F_ is the number of faces in the mesh.

## Setting custom indices
In most cases, if you need custom element indices you can just maintain you own outside of the mesh data structure (as shown [above](indexing.md)).

However, there are methods to set the element indices within a mesh. The indices must start at 0 and go up through the number of elements (so you cannot, for example, 1-index your vertices).

Since [mesh element types are handles to the underlying elements,](elements/#introduction) stored elements change meaning or become invalidated after reindexing. For example, if `Vertex v = mesh.vertex(0)`, then after changing vertex indices `v` would refer to the vertex whose new index is 0, which may be a different vertex. If the mesh was not [compressed](mutation/#compressed-mode) before reindexing, then `Vertex` variables may become invalid.

On the other hand, [Containers](containers.md) automatically update when you set custom indices. Here is an example of both behaviors:
```cpp
Vertex v0 = mesh.vertex(0);
Vertex v1 = mesh.vertex(1);

VertexData<size_t> oldVertexIndices = mesh.getVertexIndices();
VertexData<size_t> newVertexIndices = oldVertexIndices;
std::swap(newVertexIndices[v0], newVertexIndices[v1]); // swap indices of v0 and v1

mesh.setVertexIndices(newVertexIndices);
// now v0 refers to the mesh with new index 0, and v1 refers to the vertex with new index v1

std::cout << v0 << " has new index " << newVertexIndices[v0]
                << " and old index " << oldVertexIndices[v0] << std::endl;
// v_0 has new index 0 and old index 1

std::cout << v1 << " has new index " << newVertexIndices[v0]
                << " and old index " << oldVertexIndices[v0] << std::endl;
// v_1 has new index 1 and old index 0

```

This is also discussed in the page on [compression](mutation/#compressed-mode), where the same sort of invalidation happens.

??? func "`#!cpp void SurfaceMesh::setVertexIndices(const VertexData<size_t>& newIndices)`"

    Sets vertex indices to `newIndices`.

??? func "`#!cpp void SurfaceMesh::setEdgeIndices(const EdgeData<size_t>& newIndices)`"

    Sets edge indices to `newIndices`.

??? func "`#!cpp void SurfaceMesh::setFaceIndices(const FaceData<size_t>& newIndices)`"

    Sets face indices to `newIndices`.

### Setting custom edge orientations

!!! warning
    In a `ManifoldSurfaceMesh` edge orientations are set by reindexing the halfedges, so all of the caveats about [reindexing possibly invalidating elements](#setting-custom-indices) apply

??? func "`#!cpp void SurfaceMesh::setEdgeOrientations(const EdgeData<Halfedge>& newOrientations)`"

    Sets `e.halfedge()` to `newOrientation[e]` for every edge `e` in the mesh.
    `newOrientations` must obey two conditions:
    
    - `newOrientations[e].edge() == e` for every edge, since `e.halfedge().edge()` must equal `e`
    - `newOrientations[e].isInterior()`, meaning that for any boundary edge `e`, `newOrientations[e]` must be its interior halfedge and not the exterior halfedge.
