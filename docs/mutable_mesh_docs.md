# Mutating the Halfedge Mesh

The primary `HalfedgeMesh` type in geometry central (as well as the associated `MeshData<T>` types) stores elements in densely packed arrays, which yields significant performance benefits. However, unlike most dense representations, this mesh is also mutable, meaning you can insert/delete mesh elements and the packed representation will be maintained under the hood. These modifications are implemented in O(1) time via amortized resizing of the underlying buffers (a la std::vector), as well as tombstones to represent deleted elements.

Most of the complexity is hidden under the hood, and the net effect is pretty great: you can utilize almost all of the functionality of geometry-central even while mutating a mesh, and there is no performance penalty compared to an unmodifiable mesh implmentation. However, there are a few rules you need to follow as a user-- I'm documenting these here for now.

Beware: As of Oct 2018, these routines are still recently implemented and not well-tested. They may very well hold obnoxious lurking bugs.

## Modifications and Invalidation

`HalfedgeMesh.h` contains a variety of the usual mesh modification routines like `triangulate()`, `insertVertex()`, and `collapseEdge()`. You should not modify the mesh except via these functions. Whenever any of these mutation routines is invoked, the following things happen:

- Any element pointers (`Vertex`, etc) are immediately invalidated.
- Any element iterators are immediately invalidated (as in `for(Vertex v : mesh->vertices())`, or `for(Face n : f.adjacentFaces())`.
- If the operation involved deletion, the mesh is no longer _compressed_ (see below).
- The mesh is no longer _canonical_ (see below).

Note that the invalidation rules here are similar to `std::vector`-- the invalidation only _actually_ happens rarely, when an internal buffer needs to be resized. However, as the user you don't know when this will happen, so you must assume it happens every time.

### Dynamic Pointers

It's a bit inconvenient that mesh pointers get invalidated after mutation calls.  However, the mutation routines return a normal, valid pointer to some relevant mesh element after the modifcation, and with that I've found that you can quite often get away without actually needing to preserve any pointers across a call.

Nonetheless, sometimes you will need to track mesh elements across mutations. For that purpose, there are equivalent _dynamic pointer_ classes (`DynamicVertex` etc), which work just like the usual mesh pointers, except that they are not invalidated by mutation. Of course, dynamic pointers can be converted to traditional pointers, and vice-versa. Use these pointers to keep track of particular elements across mutations.

We don't recommended using dynamic pointers for everything, because they are noticebly more expensive in practice, as they use callbacks to respond to buffer resizes and permutations (though they are still amortized O(1)). 

TODO as of Sep 29, 2018 dynamic pointers are NOT preserved through compress() events, this still needs to be implemented.

### Compressing the mesh

As deletions occur, the mesh will become sparse. Most outward-facing routines (like `mesh->nVertices()`, or iterating over all elements) handle this internally. However, there are a few operations that require the mesh to be compressed, such as accessing elements by index, or converting `MeshData<T>.toVector()`. Furthermore, if many, many deletions have been performed, the mesh structure may be occupying significant extra memory.

To compress the mesh, call `mesh->compress()`. The invalidation rules mentioned above apply to this operation, and `MeshData<>` containers will be automatically updated.

### Canonicalizing the mesh

After constructing an initial halfedge mesh from polygon soup, there is a useful canonical structure in the way elements are ordered. That is, the ordering of halfedges is the same as if you traversed the faces in order, and traversed the halfedges within each face. The ordering of the edges corresponds to the first time one of the adjacent halfedges appears in the halfedge ordering (TODO reference some docs on this).

As we modify the mesh, this canonical ordering is not necessarily preserved. The mesh is not "broken" in any sense, but you might like to make use of this canonical ordering in your code. For instance, the visualization library Polyscope uses this ordering to abstract over halfedge mesh data structures, so a mesh must be canonical when used with Polyscope.

To permute the elements in to the canonical ordering, call `mesh->canonicalize()`. The invalidation rules mentioned above apply to this operation, and `MeshData<>` containers will be automatically updated.

## MeshData Objects

`MeshData<T>` contains are automatically updated as their corresponding mesh is modified. They can _always_ be indexed with a `Vertex` (etc), even if the corresponding mesh is not compressed, or is not canonical. This is very useful, because it allows you to write algorithms that track data while modifying a mesh without any extra effort (or storing data on element structs).

Indexing a `MeshData<T>` container with an integer index is only well-defined if the underlying mesh is compressed. Same for calling `MeshData<T>.toVector()`.

What happens when you create a `VertexData<T>` container, put some stuff in it, then insert a new vertex in to the mesh? What value does `vData[newV]` hold? The answer is the default value for the container. This value can specified at construction (`VertexData<int>(mesh, 42)`). If it is not specified, it is the default initialization for that type.


## Misc Notes

- `Edge.flip()` is special, and not subject to nearly any of the restrictions mentioned here, because edge flips neither insert nor delete elements. However, flips do invalidate the canonical ordering.
- Under the hood, containers and pointers are maintained by registering callbacks with the mesh, which the mesh invokes upon buffer resize operations. To implmement your own container obeying these semantics, check out `mesh->vertexExpandCallbackList` (etc).
