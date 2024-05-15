`SimplePolygonMesh` is small mesh helper class which stores vertex positions, and indexed list of faces, and (optionally) texture coordinates. Used mainly as an intermediate representation for input and output.

`#include "geometrycentral/surface/simple_polygon_mesh.h"`

!!! note "This is not the mesh you're looking for"

    `SurfaceMesh` is the [main mesh class in geometry-central](../../surface_mesh/basics/). It supports a wide range of traversals, geometry, and other operations, and is used in all higher-level algorithms. 

    This class, `SimplePolygonMesh`, is a small utility used mainly for input/output and occasionally as an implementation detail for other algorithms.

    If you are implementing any nontrivial algorithm, you should almost certainly be using `SurfaceMesh` instead.

!!! warning "PolygonSoupMesh deprecation"

    In a previous version of geometry-central, this class was called `PolygonSoupMesh`---that name was potentially misleading. To avoid breaking existing code, `polygon_soup_mesh.h` now simply typedefs `PolygonSoupMesh` as `SimplePolygonMesh`. Please use `SimplePolygonMesh` in all new code.

### Members

  - `std::vector<Vector3> vertexCoordinates` 3D positions for each vertex in the mesh. 
  - `std::vector<std::vector<size_t>> polygons` The list of polygonal faces comprising the mesh. Each inner vector is a face, given by the 0-based vertex indices in to the `vertexCoordinates` array. The ordering of these indices is interpreted as the orientation of the face, via a counter-clockwise ordering of the vertices.
  - `std::vector<std::vector<Vector2>> paramCoordinates` (optional) 2D parameterization coordinates associated with each corner of each face. If non-empty, the dimensions of this array should be exactly the same as `polygons`; each coordinate corresponds to the matching polygon corner in `polygons`.
   

### Constructors

??? func "`#!cpp SimplePolygonMesh::SimplePolygonMesh()`"

    Create an empty mesh.

??? func "`#!cpp SimplePolygonMesh::SimplePolygonMesh(std::string filename, std::string type = "")`"

    Load a mesh from file. See documentation for `readMeshFromFile()` below.

??? func "`#!cpp SimplePolygonMesh::SimplePolygonMesh(std::istream& in, std::string type)`"

    Load a mesh from file. See documentation for `readMeshFromFile()` below.

??? func "`#!cpp SimplePolygonMesh::SimplePolygonMesh(const std::vector<std::vector<size_t>>& polygons_, const std::vector<Vector3>& vertexCoordinates_)`"

    Construct a mesh from a list of polygons and vertex coordinates.

??? func "`#!cpp SimplePolygonMesh::SimplePolygonMesh(const std::vector<std::vector<size_t>>& polygons_, const std::vector<Vector3>& vertexCoordinates_, const std::vector<std::vector<Vector2>>& paramCoordinates_)`"

    Construct a mesh from a list of polygons, vertex coordinates, and parameterization coordinates.

??? func "`#!cpp std::unique_ptr<SimplePolygonMesh> unionMeshes(const std::vector<SimplePolygonMesh>& meshes)`"

    Union a collection of polygon meshes in to a single mesh.

### Input & Output
    
Currently the following mesh types are supported for reading:
    
  - `obj`
  - `ply` (using [hapPLY](https://github.com/nmwsharp/happly))
  - `stl`
  - `off`

Currently, UV coordinates are only read from file for the `obj` format.

For writing, only `obj` is supported. Note that `RichSurfaceMeshData` [provides additional ply writing routines](../io/) for surface meshes, including storing buffers of data.

The desired file type can be passed as the `type` to any of the routines below as an unpunctuated lower-case string, like `#!cpp readMeshFromFile(myFilename, "stl")`. 

??? func "`#!cpp void SimplePolygonMesh::readMeshFromFile(std::istream& in, std::string type);`"

    Read a mesh from an input stream. Type must be specificed (see above). 
    
    Clears any existing data from the file before loading.

??? func "`#!cpp void SimplePolygonMesh::readMeshFromFile(std::string filename, std::string type = "");`"

    Read a mesh from a file. `filename` should be the full path to the file. The type can be manually specified (see above), or given as the empty string (`""`) to attempt to auto-detect from the filename extension.

    Clears any existing data from the file before loading.

??? func "`#!cpp void SimplePolygonMesh::writeMesh(std::ostream& out, std::string type);`"

    Write a mesh to an input stream. Type must be specificed (see above).

??? func "`#!cpp void SimplePolygonMesh::writeMesh(std::string filename, std::string type = "");`"

    Write a mesh to file. `filename` should be the full path to the file. The type can be manually specified (see above), or given as the empty string (`""`) to attempt to auto-detect from the filename extension.

### Accessors


??? func "`#!cpp size_t SimplePolygonMesh::nFaces()`"

    The number of faces in the mesh.

??? func "`#!cpp size_t SimplePolygonMesh::nVertices()`"

    The number of vertices in the mesh.

??? func "`#!cpp bool SimplePolygonMesh::hasParameterization()`"

    True if the mesh has a 2D corner parameterization in the `paramCoordinates` member.
 

### Modification

??? func "`#!cpp void SimplePolygonMesh::mergeIdenticalVertices()`"

    Vertices with identical coordinates are merged to be a single vertex entry, and the face indices are updated accordingly.

    Identity is tested using a simple exact floating-point comparison test, no radius or threshold is supported.


??? func "`#!cpp std::vector<size_t> SimplePolygonMesh::stripUnusedVertices()`"

    Remove vertices from `vertexCoordinates` which do not appear in any face. Face indices are updated accordingly. 
  
    Returns an index translation vector mapping old indices to new, such that `vec[ind_old] == ind_new`. Holds `INVALID_IND` for any old indices which were stripped.


??? func "`#!cpp void SimplePolygonMesh::stripFacesWithDuplicateVertices()`"

    Remove any faces from `polygons` for which some vertex index appears multiple times.

??? func "`#!cpp void SimplePolygonMesh::stripInvalidFaces()`"

    Remove any faces from `polygons` which:
        - have < 3 vertices
        - OR have out of bounds vertex indices
        - OR have repated vertices

??? func "`#!cpp void SimplePolygonMesh::stripDuplicateFaces()`"

    Remove repeated faces from `polygons`, defined as having identical vertex sets even if the permutation of the vertices may be different.

??? func "`#!cpp void SimplePolygonMesh::consistentlyOrientFaces(bool outwardOrient=true)`"

    Greedily re-orient faces by flipping them to have consistent orientation within each manifold subpatch.

    The inward/outward orientation is arbitrary, unless `outwardOrient=true`. In this case, an (approximate) geometric raycasting check will be performed to attempt to orient the patches to face outward. This optional check only supports triangular meshes.
    
    In the case of non-orientable patches, there will be a not-matching boundary where greedy outward orientation fails.

??? func "`#!cpp void SimplePolygonMesh::triangulate()`"

    All faces with more than 3 sides are triangulated, replacing the original face with multiple triangular faces. Currently uses a naive fan triangulation strategy.

    An error is thrown for any faces with < 3 sides.

??? func "`#!cpp void SimplePolygonMesh::clear()`"

    Empty all data arrays for the mesh.
