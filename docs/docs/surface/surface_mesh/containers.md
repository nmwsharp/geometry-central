The halfedge mesh class is equipped with a system of containers for associating data with mesh vertices, halfedges, edges, and faces. For instance, to represent a scalar value at vertices, or a vector value at faces, one can use 

```cpp
// on vertices
VertexData<double> myVertexScalar(mesh);
Vertex v = /* some vertex */;
myVertexScalar[v] = 42.;

// on faces
FaceData<Vector3> myFaceVector(mesh);
Face f = /* some face */;
myFaceVector[f] = Vector3{1., 2., 3.};
```
and so on.


A key feature of the `MeshData<>` containers is that they **automatically adapt to mutation of the underlying mesh**. All existing `MeshData<>` containers will remain valid during any sequence of mesh element insertions and deletions, adaptively and efficiently resizing themselves as needed. These containers can also be [automatically written to file](/surface/utilities/io/#rich-surface-mesh-data).


## Mesh data types

The mesh data types are all templated on a common base class: `MeshData<E,T>`, where `E` is an element pointer type (such as `Vertex`) and `T` is a scalar type (such as `double`). The first template argument should usually be omitted in user code; the various element containers are all typedef'd with concise names as follows:

- `VertexData<T>` data at vertices
- `HalfedgeData<T>` data at (interior and exterior) halfedges 
- `CornerData<T>` data at corners
- `EdgeData<T>` data at edges 
- `FaceData<T>` data at faces 
- `BoundaryLoopData<T>` data at boundary loops

Most functionality is identical between all of these classes, so the sections below are written in terms of the generic `MeshData<>` class.

## Construction


??? func "`#!cpp MeshData<E,T>::MeshData<E,T>(SurfaceMesh& mesh)`"
    Construct a new container over a mesh. Elements will be default-initialized with `T()`.

??? func "`#!cpp MeshData<E,T>::MeshData<E,T>(SurfaceMesh& mesh, T initVal)`"
    Construct a new container over a mesh. 
    
    Elements will be initialized with `initVal`, and any newly-created mesh elements will have their default values set to `initVal`.

Additionally, see the vector-based initializers in [vector interoperability](containers.md#vector-interoperability).


## Accessors

??? func "`#!cpp T& MeshData<E,T>::operator[](E ptr)`"

    Access data stored in the container with a reference to a mesh element. A const version also exists; expect semantics like `std::vector<>`.

    For example:
    ```cpp
    // on vertices
    VertexData<double> myVertexScalar(mesh);
    Vertex v = /* some vertex */;
    myVertexScalar[v] = 42.;
    double val = myVertexScalar[v];
    ```

??? func "`#!cpp T& MeshData<E,T>::operator[](size_t ind)`"

    Access data stored in the container by the index of a mesh element. A const version also exists; expect semantics like `std::vector<>`.

    Only valid when the underlying mesh is [compressed](mutation.md#compressed-mode).
    
    Must have `0 <= ind < N`, where `N` is the number of elements of that type.

    For example:
    ```cpp
    // on vertices
    VertexData<double> myVertexScalar(mesh);
    myVertexScalar[11] = 42.;
    double val = myVertexScalar[11];

    // equivalent to:
    double val = myVertexScalar[mesh->vertex(11)];

    ```
    

??? func "`#!cpp void MeshData<E,T>::fill(T fillVal)`"

    Fill all entries in the container with `fillVal`.

??? func "`#!cpp size_t MeshData<E,T>::size()`"

    The size of the underlying buffer for the container. In particular, the largest integer `i` such that `data[i]` is safe.

    Generally on a compressed mesh this is the same as the number of elements of type `E`, e.g. `SurfaceMesh::nVertices()`, but on an uncompressed mesh or in the presence of exterior halfedges it may be larger.

    NOTE: The behavior of this function as changed in recent versions.

??? func "`#!cpp SurfaceMesh* MeshData<E,T>::getMesh() const`"

    The mesh on which the container is defined.


## Arithmetic

`MeshData<>` containers support arithmetic operations with each other, and with scalar values. All arithmetic is applied independently to each value in the container, and is only well-defined for containers defined on the same mesh.

```cpp
// add two vertex datas together
VertexData<double> A(*mesh, 1.); // (sample data, filled with all 1's)
VertexData<double> B(*mesh, 2.); 
VertexData<double> C = A + B;

// multiply times a scalar
FaceData<double> vals(*mesh, 1.);
vals *= 12.0;

// types do not need to be the same, as long as the operation
// is well-defined
VertexData<float> scales(*mesh, 2.);
VertexData<Vector3> vecs(*mesh, Vector3{1., 2., 3.});
VertexData<Vector3> scaledVecs = scales * vecs ;
```

The binary operators `+,-,*,/,%,&,|,^,<< ,>>,&&,||` and the unary operators `+,-,!,~` are all supported, along with the matching assignment operators like `+=`. Of course, the underlying container entry types must support the operation, and the result of the operation must be compatible with the destination container.



## Vector interoperability

To support easy common-case linear algebra operations, `MeshData<>` containers support conversion to and from Eigen vector types.

The corresponding vectors are indexed according to the indices of the underlying mesh elements, or by a user-supplied index map which maps each elements to a dense set of zero-based indices.


**Construct from a vector:**

??? func "`#!cpp MeshData<E,T>::MeshData<E,T>(SurfaceMesh& mesh Eigen::Matrix<T, Eigen::Dynamic, 1> vec)`"

    Construct a new container over a mesh, with the contents of `vec`.
  

??? func "`#!cpp MeshData<E,T>::MeshData<E,T>(SurfaceMesh& mesh Eigen::Matrix<T, Eigen::Dynamic, 1> vec, MeshData<E, size_t>& indexer)`"

    Construct a new container over a mesh, with the contents of `vec`, indexed according to `indexer`.


**Fill from a vector:**

??? func "`#!cpp void MeshData<E,T>::fromVector(Eigen::Matrix<T, Eigen::Dynamic, 1> vec)`"

    Fill this container with the contents of `vec`.
    

??? func "`#!cpp void MeshData<E,T>::fromVector(Eigen::Matrix<T, Eigen::Dynamic, 1> vec, MeshData<E, size_t>& indexer)`"

    Fill this container with the contents of `vec`, indexed according to `indexer`.


**Convert to a vector:**

??? func "`#!cpp Eigen::Matrix<T, Eigen::Dynamic, 1> MeshData<E,T>::toVector()`"

    Return a new `std::vector` which holds the contents of this container.

    Detail: this vector will always be a dense listing of values per-element, regardless of whether the mesh is compressed, etc. Therefore, the contents of this vector are _not_ necessarily always identical to the raw underlying buffer via `raw()`. Even in the case of a compressed mesh, for `CornerData<>` the resulting vector will omit implicit indices for exterior "outside" corners which may exist on meshes with boundary.
    

??? func "`#!cpp Eigen::Matrix<T, Eigen::Dynamic, 1> MeshData<E,T>::toVector(MeshData<E, size_t>& indexer)`"

    Return a new vector which holds the contents of this container, indexed according to `indexer`.

    See `toVector()` for more details.
    

    
## Default values

All containers track a default value for their elements, which can optionally be set at construction; if not set it is simply `T()`. After construction this value is significant because it will be used as the value for any newly-created mesh elements if the underlying mesh is mutated. The getter and setter below allow you to modify the default value for an existing container.

??? func "`#!cpp void MeshData<E,T>::setDefault(T newDefault)`"

    Sets a new default value for the container. 
    
    Does not modify any existing data in the container.


??? func "`#!cpp T MeshData<E,T>::getDefault() const`"

    Get the current default value for the container. 


## Transferring data

`MeshData<>` containers are defined with respect to a particular mesh object. Sometimes one may need to transfer data defined on one mesh to another, for instance after making a copy of a mesh, or when reading data from file.

??? func "`#!cpp MeshData<E,T> MeshData<E,T>::reinterpretTo(SurfaceMesh& target)`"

    Map data defined on one halfedge mesh to another. The meshes must have the same number of elements, and data will be naively transferred between elements with the same index.

    Requires that both meshes be [compressed](mutation.md#compressed-mode).

    Example usage:
    ```cpp
    SurfaceMesh meshA = /* something */;
    std::unique_ptr<SurfaceMesh> meshB = meshA.copy();

    FaceData<Vector3> myDataOnA(meshA);
    /* fill myDataOnA with interesting values */

    FaceData<Vector3> myDataOnB = myDataOnA.reinterpretTo(*meshB);
    ```


## Advanced features

Under the hood, all `MeshData<>` types use a `Eigen::Matrix<T>` to store their values. However, the size and indexing indexing are carefully managed in conjunction with the underlying mesh. This vector will only be a dense listing if the mesh is [compressed](mutation.md#compressed-mode).

??? func "`#!cpp Eigen::Matrix<T, Eigen::Dynamic, 1>& MeshData<E,T>::raw()`"

    Access the raw underlying Eigen vector of storage.

??? func "`#!cpp const Eigen::Matrix<T, Eigen::Dynamic, 1>& MeshData<E,T>::raw() const`"

    Access the raw underlying Eigen vector of storage (const).


