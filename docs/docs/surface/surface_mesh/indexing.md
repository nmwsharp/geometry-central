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
