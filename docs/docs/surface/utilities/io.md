## Reading meshes

Construct a halfedge mesh from a file on disk.

`#include "geometrycentral/surface/meshio.h"`

Example usage:
```cpp
#include "geometrycentral/surface/meshio.h"
using namespace geometrycentral::surface;

std::unique_ptr<HalfedgeMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;
std::tie(mesh, geometry) = loadMesh("spot.obj"); 
```

??? func "`#!cpp std::tuple<std::unique_ptr<HalfedgeMesh>,std::unique_ptr<VertexPositionGeometry>> loadMesh(std::string filename, std::string type="")`"

    Load a mesh from file. Returns both a `HalfedgeMesh` representing the connectivity, and a `Geometry` representing the geometry. See example below to concisely unpack.

    If the file includes vertices which do not appear in any face, they will be stripped from the vertex listing and ignored.

    The `type` parameter determines the type of file to load. For example, `type="ply"` will attempt to read the target file as a .ply file. If no type is given, the type will be inferred from the file name. 

    Currently the following types are supported:
    
    - `obj`
    - `ply` (using [hapPLY](https://github.com/nmwsharp/happly))

## Writing meshes

TODO

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
