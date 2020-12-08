`#include "geometrycentral/surface/meshio.h"`


## Reading meshes

Construct a surface mesh from a file on disk, or more generally any `std::istream`.

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
    std::tie(mesh, geometry) = readSurfaceMesh("spot.obj"); 
    
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


??? func "`#!cpp std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>, std::unique_ptr<CornerData<Vector2>>> readParameterizedManifoldSurfaceMesh(std::string filename, std::string type="")`"

    Loads a manifold surface mesh plus UV (texture) coordinates from a file.  See other variants for details.

    Currently only OBJ files are supported.

??? func "`#!cpp std::tuple<std::unique_ptr<SurfaceMesh>, std::unique_ptr<VertexPositionGeometry>, std::unique_ptr<CornerData<Vector2>>> readParameterizedSurfaceMesh(std::string filename, std::string type="")`"

    Loads a general surface mesh plus UV (texture) coordinates from a file.  See other variants for details.

    Currently only OBJ files are supported.

## Writing meshes

Write a mesh to file, or more generally any `std::ostream`. 

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
| `obj` |    ✅    |    ✅    |      ✅     |                                                 |
| `ply` |    ✅    |         |            |                                                   |
| `off` |    ✅    |         |            |                                                   |
| `stl` |    ✅    |         |            | Exactly colocated vertices are automatically merged |



## Rich Surface Mesh Data

The `RichSurfaceMeshData` offers advanced IO which interoperates directly with the geometry-central mesh data structures. In particular, it has two useful features:

- Data stored in `MeshData<>` containers (`VertexData<>`, `EdgeData<>`, etc.) can be automatically written to and read from file
- The internal representation used by `SurfaceMesh` is directly written to file, improving performance and supporting meshes which cannot be represented as face index lists.

The written files can hold mesh's connectivity, associated geometry and/or any number of properties. These files can be used to load/save a whole mesh from file, or to load/save individual properties associated with some existing mesh.

Internally, data is stored as additional custom fields of a `.ply` file.  Here, we're using the `.ply` format as a general container for structured data---other software will not automatically understand the additional fields in these files.  A benefit of the `ply` format is that data can be written in both efficient binary or plaintext ascii formats---the default is binary.  The `RichSurfaceMeshData` class is used to read and write these souped-up `.ply` files, and is distinct from the simple mesh-loading `.ply` interface above. 

`#include "geometrycentral/surface/rich_surface_mesh_data.h"`

Example: loading and saving data on a surface

```cpp
#include "geometrycentral/surface/rich_surface_mesh_data.h"
using namespace geometrycentral::surface;

// Open a file and load the mesh therein
std::unique_ptr<SurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;
std::unique_ptr<RichSurfaceMeshData> richData;
std::tie(mesh, richData) = RichSurfaceMeshData::readMeshAndData("archive.ply"); 
geometry = richData->getGeometry();

// Read a stored value
FaceData<double> faceValues = richData->getFaceProperty<double>("name_a");

// Do some science
CornerData<double> cornerValues = /* something important */
EdgeData<int> edgeValues = /* something else important */

// Add these values to the reader/writer
// Note: alternately you could create a new record like
// RichSurfaceMeshData newData(*mesh);
richData->addCornerProperty("name_b", cornerValues);
richData->addEdgeProperty("name_c", edgeValues);

// Write the data to file
richData->write("new_archive.ply")
```

Example: saving and loading a surface along with some properties

```cpp
#include "geometrycentral/surface/rich_surface_mesh_data.h"
using namespace geometrycentral::surface;

// Your existing mesh and data
std::unique_ptr<SurfaceMesh> mesh; // your mesh
std::unique_ptr<VertexPositionGeometry> geometry; // your geometry
EdgeData<double> data; // your data

// Store data
RichSurfaceMeshData richData(*mesh);
richData.addMeshConnectivity();
richData.addGeometry(*geometry);
richData.addEdgeProperty("my prop", data);

// Write to file
richData.write("file.ply");
```

later, when loading...

```cpp
// Load the mesh and the data from file
std::unique_ptr<SurfaceMesh> mesh; 
std::unique_ptr<RichSurfaceMeshData> richData;
std::tie(mesh, richData) = RichSurfaceMeshData::readMeshAndData("file.ply");  

std::unique_ptr<VertexPositionGeometry> geometry = richData->getGeometry();

EdgeData<double> data = richData->getEdgeProperty<double>("my prop");
```

### Writing rich data

These methods add properties to the `RichSurfaceMeshData` object, which will be written when `write()` is called. The set of scalar types supported is the same as the [.ply](https://github.com/nmwsharp/happly) format, including list types.  For instance, a property of type `double` on vertices could written to a new `ply` file with.

```cpp
RichSurfaceMeshData richData(*mesh);
VertexData<double> values = /* incredibly important data */
richData.addVertexProperty("important_values", values);
richData.write("my_file.ply");
```

The following routines can create the data object and write it to file:

??? func "`#!cpp RichSurfaceMeshData::RichSurfaceMeshData(SurfaceMesh& mesh)`"

    Construct a data object from an existing mesh.


??? func "`#!cpp void RichSurfaceMeshData::write(std::string filename)`"

    Write the object to file. The binary/ascii writing mode is determined by the `RichSurfaceMeshData::outputFormat` option.
    
    Note that if this reader/writer was created by loading a file, and is later written using `write()` all fields from the initial file will be automatically written out.

??? func "`#!cpp void RichSurfaceMeshData::write(std::ostream& out)`"

    Write the object to stream. The binary/ascii writing mode is determined by the `RichSurfaceMeshData::outputFormat` option.
    
    Note that if this reader/writer was created by loading a file, and is later written using `write()` all fields from the initial file will be automatically written out.

To store the mesh connectivity itself in the file, call `addMeshConnectivity()`---this is required if you want to load the mesh from the file later. Similarly, the `addGeometry()` helpers will store geometry as vertex positions or edge lengths.

??? func "`#!cpp void RichSurfaceMeshData::addMeshConnectivity()`"

    Store the meshes connectivity in the file.

    This routine always stores both the internal `SurfaceMesh` representation (for use with geometry-central), and a traditional face-index list representation (for use with other software).

??? func "`#!cpp void RichSurfaceMeshData::addGeometry(const EmbeddedGeometryInterface& geometry)`"

    Add geometry to the record, which will written as `double` vertex coordinates properties named "x", "y", and "z".

??? func "`#!cpp void RichSurfaceMeshData::addIntrinsicGeometry(const IntrinsicGeometryInterface& geometry)`"

    Add geometry to the record, which will written as a `double` edge propery called `intrinsic_edge_lengths`.

General properties can then be written as:

??? func "`#!cpp void RichSurfaceMeshData::addVertexProperty<>(std::string name, const VertexData<>& data)`"
    Add a property for writing. 

    - `name` A user-defined name with which the property will be written to file. Must be unique within the file.
    - `data` The data to be written, such a `VertexData<double>`.

??? func "`#!cpp void RichSurfaceMeshData::addHalfedgeProperty<>(std::string name, const HalfedgeData<>& data)`"
    Add a property for writing. 

    - `name` A user-defined name with which the property will be written to file. Must be unique within the file.
    - `data` The data to be written, such a `HalfedgeData<double>`. 

??? func "`#!cpp void PlyCornerMeshData::addCornerProperty<>(std::string name, const CornerData<>& data)`"
    Add a property for writing. 

    - `name` A user-defined name with which the property will be written to file. Must be unique within the file.
    - `data` The data to be written, such a `CornerData<double>`. 

??? func "`#!cpp void RichSurfaceMeshData::addEdgeProperty<>(std::string name, const EdgeData<>& data)`"
    Add a property for writing. 

    - `name` A user-defined name with which the property will be written to file. Must be unique within the file.
    - `data` The data to be written, such a `EdgeData<double>`.

??? func "`#!cpp void RichSurfaceMeshData::addFaceProperty<>(std::string name, const FaceData<>& data)`"
    Add a property for writing. 

    - `name` A user-defined name with which the property will be written to file. Must be unique within the file.
    - `data` The data to be written, such a `FaceData<double>`. 

??? func "`#!cpp void RichSurfaceMeshData::addBoundaryLoopProperty<>(std::string name, const BoundaryLoopData<>& data)`"
    Add a property for writing. 

    - `name` A user-defined name with which the property will be written to file. Must be unique within the file.
    - `data` The data to be written, such a `BoundaryLoopData<double>`. 


### Reading rich data

When reading a `RichSurfaceMeshData` object, you have two options:

**Option A** Open the object _on_ an existing mesh, like `RichSurfaceMeshData(*mesh, "file.ply")`. All resulting properties will be defined on this mesh. The existing mesh must be the "same" as the one the file was saved from, with the same number of elements in the same semantic ordering.

??? func "`#!cpp RichSurfaceMeshData::RichSurfaceMeshData(SurfaceMesh& mesh, std::string filename)`"

    Open a `ply` file, and interpret its fields as living on the existing halfedge mesh `mesh`. Any read properties will be returned in containers defined on `mesh`.

??? func "`#!cpp RichSurfaceMeshData::RichSurfaceMeshData(SurfaceMesh& mesh, std::istream& in)`"

    Same as above, loading from a general `istream`.

**Option B** Simultaneously construct a new mesh from the file, and open the file _on_ that mesh, via `readMeshAndData(...)`. The file must have been saved with mesh connectivity included by calling `addMeshConnectivity()`.

??? func "`#!cpp static std::tuple<std::unique_ptr<SurfaceMesh>, std::unique_ptr<RichSurfaceMeshData>> RichSurfaceMeshData::readMeshAndData(std::string filename)`"

    Convenience factory function to open a rich `.ply` file and load the mesh contained within, as well as creating a `RichSurfaceMeshData` reader/writer to access any other properties stored in the file.

    The base class of the created `SurfaceMesh` will match the mesh from which is was created.


??? func "`#!cpp static std::tuple<std::unique_ptr<SurfaceMesh>, std::unique_ptr<RichSurfaceMeshData>> RichSurfaceMeshData::readMeshAndData(std::istream& in)`"

    Same as above, loading from a stream.


??? func "`#!cpp static std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<RichSurfaceMeshData>> RichSurfaceMeshData::readManifoldMeshAndData(std::string filename)`"

    Convenience factory function to open a rich `.ply` file and load the mesh contained within, as well as creating a `RichSurfaceMeshData` reader/writer to access any other properties stored in the file.

    Returns a _manifold_ mesh, erroring out out if the file was not saved from a `ManifoldSurfaceMesh`.

??? func "`#!cpp static std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<RichSurfaceMeshData>> RichSurfaceMeshData::readManifoldMeshAndData(std::istream& in)`"
    
    Same as above, loading from a stream.


Once the `RichSurfaceMeshData` has been opened, properties can then be read from the file like `getVertexProperty<double>(name)`, etc.  The template argument to this function will likely be necessary to resolve the expected type of the data. For instance, a property of type `double` on vertices could be accessed with.

```cpp
RichSurfaceMeshData data(mesh, "my_file.ply");
VertexData<double> values = data.getVertexProperty<double>("important_values");
```

The automatic type promotion in [hapPLY](https://github.com/nmwsharp/happly) gives some flexibility in specifying the type of the read data--- for instance if property `"propName"` in the example above was stored as a `float`, it could still be read as a `double`. See the documentation there for details.

Properties can then be read as:

??? func "`#!cpp VertexData<T> RichSurfaceMeshData::getVertexProperty<T>(std::string name)`"
    Read a property from a loaded file. 

    - `name` A user-defined name with which the property will be read from the file. Throws if no such property exists.
    - **Return:** The requested container.

??? func "`#!cpp HalfedgeData<T> RichSurfaceMeshData::getHalfedgeProperty<T>(std::string name)`"
    Read a property from a loaded file. 

    - `name` A user-defined name with which the property will be read from the file. Throws if no such property exists.
    - **Return:** The requested container.

??? func "`#!cpp CornerData<T> PlyCornerMeshData::getCornerProperty<T>(std::string name)`"
    Read a property from a loaded file. 

    - `name` A user-defined name with which the property will be read from the file. Throws if no such property exists.
    - **Return:** The requested container.

??? func "`#!cpp EdgeData<T> RichSurfaceMeshData::getEdgeProperty<T>(std::string name)`"
    Read a property from a loaded file. 

    - `name` A user-defined name with which the property will be read from the file. Throws if no such property exists.
    - **Return:** The requested container.

??? func "`#!cpp FaceData<T> RichSurfaceMeshData::getFaceProperty<T>(std::string name)`"
    Read a property from a loaded file. 

    - `name` A user-defined name with which the property will be read from the file. Throws if no such property exists.
    - **Return:** The requested container.

??? func "`#!cpp BoundaryLoopData<T> RichSurfaceMeshData::getBoundaryLoopProperty<T>(std::string name)`"
    Read a property from a loaded file. 

    - `name` A user-defined name with which the property will be read from the file. Throws if no such property exists.
    - **Return:** The requested container.
  
  
??? func "`#!cpp std::unique_ptr<VertexPositionGeometry> RichSurfaceMeshData::getGeometry()`"

    Build a new geometry object from vertex postions stored in a file (by `addGeometry()`).

??? func "`#!cpp std::unique_ptr<EdgeLengthGeometry> RichSurfaceMeshData::getIntrinsicGeometry()`"

    Build a new geometry object from edge lengths stored in a file (by `addIntrinsicGeometry()`).


## Factory constructors

  These simultaneously construct the connectivity and geometry of a mesh, and are used internally in many of the subroutines above.

  `#include "geometrycentral/surface/surface_mesh_factories.h"`

  Construct a mesh and geometry from a list of polygons and vertex positions:

??? func "`#!cpp std::tuple<std::unique_ptr<SurfaceMesh>, std::unique_ptr<VertexPositionGeometry>> makeSurfaceMeshAndGeometry(const std::vector<std::vector<size_t>>& polygons, const std::vector<Vector3> vertexPositions)`"

    Construct a new surface mesh and geometry from a list of face indices and a list of vertex positions. See the constructors of `SurfaceMesh` and `VertexPositionGeometry` for more details about their semantics.
    
    - `polygons` is a nested list of zero-indexed face indices 
    - `vertexPositions` is a list of 3D vertex positions

??? func "`#!cpp std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>> makeManifoldSurfaceMeshAndGeometry(const std::vector<std::vector<size_t>>& polygons, const std::vector<Vector3> vertexPositions)`"

    Same as above, but the result is a `ManifoldSurfaceMesh` (and thus the connectivity must describe a manifold mesh).


  Construct a mesh and geometry from face indices and vertex positions, stored in dense (Eigen) matrices:

  ```cpp
  #include "geometrycentral/surface/surface_mesh_factories.h"

  // matrices describing mesh (populated somehow)
  Eigen::MatrixXd vMat = /* ... */;  // V x 3 array of vertex positions
  Eigen::MatrixXi fMat = /* ... */;  // F x 3 array of zero-indexed face indices 
                                     // (or F x 4 for quads, etc)

  // construct geometry-central mesh types
  std::unique_ptr<SurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geometry;
  std::tie(mesh, geometry) = makeSurfaceMeshAndGeometry(vMat, fMat);

  // OR to construct a mesh which must be manifold
  std::unique_ptr<ManifoldSurfaceMesh> meshManifold;
  std::unique_ptr<VertexPositionGeometry> geometryManifold;
  std::tie(meshManifold, geometryManifold) = makeManifoldSurfaceMeshAndGeometry(vMat, fMat);
  ```

??? func "`#!cpp std::tuple<std::unique_ptr<SurfaceMesh>, std::unique_ptr<VertexPositionGeometry>> makeSurfaceMeshAndGeometry(const Eigen::MatrixBase<Scalar_V>& vMat, const Eigen::MatrixBase<Scalar_F>& fMat)`"

    Construct a new surface mesh and geometry from a list of face indices and a list of vertex positions. See the constructors of `SurfaceMesh` and `VertexPositionGeometry` for more details about their semantics. Note that these arguments are ordered `V,F` to match MATLAB & friends conventions.

    - `vMat` is a `Vx3` floating-point valued matrix of vertex positions. Any floating point type can be used.
    - `fMat` is an `FxD` index valued matrix of zero-indexed face indices (e.g. an Fx3 array of triangles, of Fx4 array of quads). The index type can be any integer type (like `size_t` or `int`).
    
    The `Eigen:MatrixBase<T>` type is just a general type which accepts most Eigen matrix types as input, including geometry-central's nicely-named wrapper `DenseMatrix<T>`.

??? func "`#!cpp std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>> makeManifoldSurfaceMeshAndGeometry(const Eigen::MatrixBase<Scalar_V>& vMat, const Eigen::MatrixBase<Scalar_F>& fMat)`"

    Same as above, but the result is a `ManifoldSurfaceMesh` (and thus the connectivity must describe a manifold mesh).








