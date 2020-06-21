`#include "geometrycentral/surface/meshio.h"`


## Reading meshes

Construct a surface mesh from a file on disk, or more generally any `std::istream`.

Example usage:
```cpp
#include "geometrycentral/surface/meshio.h"
using namespace geometrycentral::surface;

// Load a surface mesh which is required to be manifold
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;
std::tie(mesh, geometry) = readManifoldSurfaceMesh("spot.obj"); 


// Load a more general surface mesh, which might be nonmanifold
std::unique_ptr<SurfaceMesh> mesh2;
std::unique_ptr<VertexPositionGeometry> geometry2;
std::tie(mesh2, geometry2) = readSurfaceMesh("spot_messy.obj"); 
```


??? note "Why use `std::unique_ptr<>`?"

    The mesh loader, like many functions in geometry-central, returns constructed objects via a `unique_ptr`. Unique pointers are an important tool for memory management in modern C++; if you haven't used them before, we suggest you give them a try!

    In most ways, a `unique_ptr` acts just like a normal C++ pointer. You can dereference it with `*uPtr`, and access its members and function like `uPtr->function()`. However, the `unique_ptr` helps prevent common memory-management mistakes, and communicates the programmer's intent about object lifetime. This is accomplished with two properties:

      - You don't need to call `delete` on a `unique_ptr`, it happens automatically when the pointer is destructed, e.g. when it goes out of scope at the end of a function, or when the object it is a member of gets deleted. This helps prevent memory leaks where you forget to deallocate the object.

      - You cannot copy the `unique_ptr`; hence it is "unique"! You can still pass around references, or `std::move()` the pointer, which are sufficient for most reasonable uses. This helps prevent you from creating a copy, and then accidentally deleting the pointer twice.

    ---

    The general paradigm in geometry-central (and a recommended style in modern C++) is to return long-lived, allocated objects with a `unique_ptr`, and pass these objects in to functions and dependent classes by reference.

    For instance, we might write a function which takes a mesh as an argument like

    ```cpp
    void processMesh(SurfaceMesh& inputMesh) { /* do stuff */}
    ```
   
    and call it by dereferencing the unique pointer to pass a reference

    ```cpp
    std::unique_ptr<SurfaceMesh> mesh;
    std::unique_ptr<VertexPositionGeometry> geometry;
    std::tie(mesh, geometry) = loadMesh("spot.obj"); 
    
    processMesh(*mesh);
    ```

    For more details about unique pointers, see the [language documentation](https://en.cppreference.com/w/cpp/memory/unique_ptr), or many tutorials around the web.

    --- 

    If you really don't want to use unique pointers, you can simply release the unique pointer to an ordinary pointer:

    ```cpp
    std::unique_ptr<SurfaceMesh> mesh /* populated as above */;
    SurfaceMesh* meshPtr = mesh.release();
    ```
    
    The `meshPtr` now points the mesh object, and you are responsible for eventually deleting this pointer. After calling `release()`, the unique pointer points to nothing and will no longer deallocate the object.


??? func "`#!cpp std::tuple<std::unique_ptr<SurfaceMesh>, std::unique_ptr<VertexPositionGeometry>> readSurfaceMesh(std::string filename, std::string type = "")`"

    Load a general, not-necessarily-manifold mesh from file. Returns both a `SurfaceMesh` representing the connectivity, and a `Geometry` representing the geometry. See the example above to concisely unpack.

    If the file includes vertices which do not appear in any face, they will be stripped from the vertex listing and ignored.

    The `type` parameter determines the type of file to load. For example, `type="ply"` will attempt to read the target file as a .ply file. If no type is given, the type will be inferred from the file extension. 

    See the matrix below for all supported file types.

??? func "`#!cpp std::tuple<std::unique_ptr<SurfaceMesh>, std::unique_ptr<VertexPositionGeometry>> readSurfaceMesh(std::istream& in, std::string type)`"

    Same as above, but reads from a general `istream` object rather than a filename. 
    
    When reading from a stream, the type _must_ be specified and cannot be automatically inferred.


??? func "`#!cpp std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>> readManifoldSurfaceMesh(std::string filename, std::string type = "")`"
    
    Load a manifold mesh from file. Returns both a `ManifoldSurfaceMesh` representing the connectivity, and a `Geometry` representing the geometry. See the example above to concisely unpack. If the mesh to be loaded is not manifold, an exception will be thrown.

    If the file includes vertices which do not appear in any face, they will be stripped from the vertex listing and ignored.

    The `type` parameter determines the type of file to load. For example, `type="ply"` will attempt to read the target file as a .ply file. If no type is given, the type will be inferred from the file extension. 

    See the matrix below for all supported file types.


??? func "`#!cpp std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>> readManifoldSurfaceMesh(std::istream& in, std::string type)`"

    Same as above, but reads from a general `istream` object rather than a filename. 

    When reading from a stream, the type _must_ be specified and cannot be automatically inferred.



## Writing meshes

Write a mesh to file, or more generally any `std::ostream`. 

Example usage:
```cpp
#include "geometrycentral/surface/meshio.h"
using namespace geometrycentral::surface;

// Write a surface mesh
// Recall that ManifoldSurfaceMesh is a subclass of SurfaceMesh, so
// it can be used just the same
std::unique_ptr<SurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;
writeSurfaceMesh(*mesh, *geometry, "my_mesh.obj"); 
```


??? func "`#!cpp void writeSurfaceMesh(SurfaceMesh& mesh, EmbeddedGeometryInterface& geometry, std::string filename, std::string type = "")`"
    
    Write a mesh to file. 
    
    The `type` parameter determines the type of file to write. For example, `type="obj"` will write the target file as a .obj file. If no type is given, the type will be inferred from the file extension. 

    See the matrix below for all supported file types.


??? func "`#!cpp void writeSurfaceMesh(SurfaceMesh& mesh, EmbeddedGeometryInterface& geometry, CornerData<Vector2>& texCoords, std::string filename, std::string type = "") `"

    Write a mesh to file as above, except corner per-corner coordinates are also written, if supported by the file format.

??? func "`#!cpp void writeSurfaceMesh(SurfaceMesh& mesh, EmbeddedGeometryInterface& geometry, std::ostream& out, std::string type) `"


    Write a mesh to file as above, to a general `std::ostream` rather than directly to a named file.
    
    When writing to a stream, the type _must_ be specified and cannot be automatically inferred.
    

??? func "`#!cpp void writeSurfaceMesh(SurfaceMesh& mesh, EmbeddedGeometryInterface& geometry, CornerData<Vector2>& texCoords, std::ostream& out, std::string type) `"

    Write a mesh to file as above, to a general `std::ostream` and with the specified texture coordinates, if the file format supports it.

    When writing to a stream, the type _must_ be specified and cannot be automatically inferred.


### Packing scalar data

A handy trick for transferring data between programs (e.g., to create visualizations) is to pack scalar data in to texture coordinates when writing a mesh to file. To make this easier, the helper `packToParam()` stores data as corner coordinates which can be passed to the writers.

```cpp
#include "geometrycentral/surface/meshio.h"
using namespace geometrycentral::surface;

std::unique_ptr<SurfaceMesh> mesh; // some mesh
std::unique_ptr<VertexPositionGeometry> geometry; // some geometry
VertexData<double> vals; // some values at vertices

writeSurfaceMesh(*mesh, *geometry, packToParam(*mesh, vals), "my_mesh.obj"); 
```


??? func "`#!cpp CornerData<Vector2> packToParam(SurfaceMesh& mesh, VertexData<double>& vals)`"

    Create a `CornerData` with the value from `vals` stored in the first coordinate at each corner. The second coordinate will hold `0`s.

??? func "`#!cpp CornerData<Vector2> packToParam(SurfaceMesh& mesh, VertexData<double>& valsX, VertexData<double>& valsY)`"
    
    Create a `CornerData` with the value from `valsX` stored in the first coordinate at each corner, and `valsY` stored in the second coordinate.


## Supported file types

| key | reading | writing | tex coords | notes                                                |
|-----|:-------:|:-------:|:----------:|------------------------------------------------------|
| `obj` |    ✅    |    ✅    |      ✅     |                                                      |
| `ply` |    ✅    |         |            |                                                      |
| `off` |    ✅    |         |            |                                                      |
| `stl` |    ✅    |         |            | Exactly coincident vertices are automatically merged |


## Serializing containers 

Data stored in `MeshData<>` containers can be automatically written and loaded from file. Internally, data is stored as additional custom fields of a `.ply` file.  Here, we're using the `.ply` format as a general container for structured data---other software may not automatically understand the additional fields in these files.

The `PlyHalfedgeMeshData` class is used to read and write these souped-up `.ply` files, and is distinct from the simple mesh-loading `.ply` interface above.

`#include "geometrycentral/surface/ply_halfedge_mesh_data.h"`

Example usage:

TODO

```cpp
#include "geometrycentral/surface/ply_halfedge_mesh_data.h"
using namespace geometrycentral::surface;

// Open a file and load the mesh therein
std::unique_ptr<HalfedgeMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;
std::unique_ptr<PlyHalfedgeMeshData> plyData;
std::tie(mesh, geometry) = PlyHalfedgeMeshData::loadMeshAndData("archive.ply"); 

// Read a stored value
FaceData<double> faceValues = plyData->getFaceProperty<double>("name_a");

// Do some science
CornerData<double> cornerValues = /* something important */
EdgeData<int> edgeValues = /* something else important */

// Add these values to the reader/writer
// note: alternately could create a new record like
// PlyHalfedgeMeshData newData(*mesh);
plyData->addCornerProperty("name_b", cornerValues);
plyData->addEdgeProperty("name_c", edgeValues);

// Write the data to file
plyData->write("new_archive.ply")
```

??? func "`#!cpp PlyHalfedgeMeshData::PlyHalfedgeMeshData(HalfedgeMesh& mesh, std::string filename, bool verbose = false)`"

    Open a `ply` file, and interpret its fields as living on the existing halfedge mesh `mesh`. Any read properties will be returned in containers defined on `mesh`.

??? func "`#!cpp PlyHalfedgeMeshData::PlyHalfedgeMeshData(HalfedgeMesh& mesh_, bool verbose = false)`"

    Construct from an existing mesh. The mesh connectivity will be included when writing the file.

??? func "`#!cpp static std::tuple<std::unique_ptr<HalfedgeMesh>, std::unique_ptr<PlyHalfedgeMeshData>> PlyHalfedgeMeshData::loadMeshAndData(std::string filename, bool verbose = false)`"

    Convenience factory function to open a `.ply` file and load the mesh contained within, as well as creating a `PlyHalfedgeMeshData` reader/writer to access any other properties stored in the file.

??? func "`#!cpp void PlyHalfedgeMeshData::write(std::string filename)`"

    Write the object to file. The binary/ascii writing mode is determined by the `PlyHalfedgeMeshData::outputFormat` option.
    
    Note that if this reader/writer was created by loading a file, and is later written using `write()` all fields from the initial file will be automatically written out.

### Writing properties

These methods add properties to the `PlyHalfedgeMeshData` object, which will be written when `write()` is called. The set of scalar types supported is the same as the [.ply](https://github.com/nmwsharp/happly) format, including list types.  For instance, a property of type `double` on vertices could written to a new `ply` file with.

```cpp
PlyHalfedgeMeshData data(mesh);
VertexData<double> values = /* incredibly important data */
data.addVertexProperty("important_values", values);
data.write("my_file.ply");
```

??? func "`#!cpp void PlyHalfedgeMeshData::addVertexProperty<>(std::string name, const VertexData<>& data)`"
    Add a property for writing. 

    - `name` A user-defined name with which the property will be written to file. Must be unique within the file.
    - `data` The data to be written, such a `VertexData<double>`.

??? func "`#!cpp void PlyHalfedgeMeshData::addHalfedgeProperty<>(std::string name, const HalfedgeData<>& data)`"
    Add a property for writing. 

    - `name` A user-defined name with which the property will be written to file. Must be unique within the file.
    - `data` The data to be written, such a `HalfedgeData<double>`. 

??? func "`#!cpp void PlyCornerMeshData::addCornerProperty<>(std::string name, const CornerData<>& data)`"
    Add a property for writing. 

    - `name` A user-defined name with which the property will be written to file. Must be unique within the file.
    - `data` The data to be written, such a `CornerData<double>`. 

??? func "`#!cpp void PlyHalfedgeMeshData::addEdgeProperty<>(std::string name, const EdgeData<>& data)`"
    Add a property for writing. 

    - `name` A user-defined name with which the property will be written to file. Must be unique within the file.
    - `data` The data to be written, such a `EdgeData<double>`.

??? func "`#!cpp void PlyHalfedgeMeshData::addFaceProperty<>(std::string name, const FaceData<>& data)`"
    Add a property for writing. 

    - `name` A user-defined name with which the property will be written to file. Must be unique within the file.
    - `data` The data to be written, such a `FaceData<double>`. 

??? func "`#!cpp void PlyHalfedgeMeshData::addBoundaryLoopProperty<>(std::string name, const BoundaryLoopData<>& data)`"
    Add a property for writing. 

    - `name` A user-defined name with which the property will be written to file. Must be unique within the file.
    - `data` The data to be written, such a `BoundaryLoopData<double>`. 

??? func "`#!cpp void PlyHalfedgeMeshData::addGeometry(const Geometry<Euclidean>& geometry)`"
    Add geometry to the record, which will written as `double` vertex coordinates properties named "x", "y", and "z".

### Reading properties

These methods read properties from the `PlyHalfedgeMeshData` object, which exist either because they were read from an opened file, or because they were previously added with the `add___()` functions above.

The template argument to this function will likely be necessary to resolve the expected type of the data. For instance, a property of type `double` on vertices could be accessed with.

```cpp
PlyHalfedgeMeshData data(mesh, "my_file.ply");
VertexData<double> values = data.getVertexProperty<double>("important_values");
```

The automatic type promotion in [hapPLY](https://github.com/nmwsharp/happly) gives some flexibility in specifying the type of the read data--- for instance if property `"propName"` in the example above was stored as a `float`, it could still be read as a `double`. See the documentation there for details.

??? func "`#!cpp VertexData<T> PlyHalfedgeMeshData::getVertexProperty<T>(std::string name)`"
    Read a property from a loaded file. 

    - `name` A user-defined name with which the property will be read from the file. Throws if no such property exists.
    - **Return:** The requested container.

??? func "`#!cpp HalfedgeData<T> PlyHalfedgeMeshData::getHalfedgeProperty<T>(std::string name)`"
    Read a property from a loaded file. 

    - `name` A user-defined name with which the property will be read from the file. Throws if no such property exists.
    - **Return:** The requested container.

??? func "`#!cpp CornerData<T> PlyCornerMeshData::getCornerProperty<T>(std::string name)`"
    Read a property from a loaded file. 

    - `name` A user-defined name with which the property will be read from the file. Throws if no such property exists.
    - **Return:** The requested container.

??? func "`#!cpp EdgeData<T> PlyHalfedgeMeshData::getEdgeProperty<T>(std::string name)`"
    Read a property from a loaded file. 

    - `name` A user-defined name with which the property will be read from the file. Throws if no such property exists.
    - **Return:** The requested container.

??? func "`#!cpp FaceData<T> PlyHalfedgeMeshData::getFaceProperty<T>(std::string name)`"
    Read a property from a loaded file. 

    - `name` A user-defined name with which the property will be read from the file. Throws if no such property exists.
    - **Return:** The requested container.

??? func "`#!cpp BoundaryLoopData<T> PlyHalfedgeMeshData::getBoundaryLoopProperty<T>(std::string name)`"
    Read a property from a loaded file. 

    - `name` A user-defined name with which the property will be read from the file. Throws if no such property exists.
    - **Return:** The requested container.
  
  
??? func "`#!cpp std::unique_ptr<Geometry<Euclidean>> PlyHalfedgeMeshData::getGeometry()`"
    Build a geometry object from vertex postions stored in a file.


## Storing Delta-complexes

Most mesh file formats store connectivity via a face-vertex list; this format used by default in all IO functions above. However, this format is insufficient for representing more general $\Delta$-complexes. To support IO for $\Delta$-complexes, connectivity can instead be encoded via halfedge adjacency indices as described in the [Internals](internals.md) section. This representation has the additional advantage that loading halfedge meshes will be very fast, as no connectivity needs to be detected.

The `.ply` readers automatically support reading this format. The option below enables writing `.ply` files in this format via the `PlyHalfedgeMeshData` class.


??? func "`#!cpp bool PlyHalfedgeMeshData::useHalfedgeAdjacency`"
    If true, writing will produce a `.ply` file which stores connectivity using haflfedge permutation indices rather than the usual face-vertex list.

    Default value: `false`.
