`#include "geometrycentral/surface/common_subdivision.h"`

Intuitively, the common subdivision of meshes $T_A$ and $T_B$ is the polygon mesh that you obtain if you cut $T_B$ along the edges of $T_A$. The vertices of the common subdivision correspond to intersections of simplices of $T_A$ with simplices of $T_B$ (e.g. intersections between edges). 
![common subdivision of two meshes](/media/common_subdivision.svg){: style="max-height: 10em; display: block; margin-left: auto; margin-right: auto;"}

The `CommonSubdivision` class represents the common subdivision of two meshes, referred to as `meshA` and `meshB`. 
!!! note 
    In the common case of intrinsic triangulations, the **input triangulation** is `meshA`, and the **intrinsic triangulation** drawn along its surface is `meshB`.

!!! warning 
    Every vertex of `meshA` must also be a vertex of `meshB`.

The common subdivision is encoded by a list of [CommonSubdivisionPoints](#common-subdivision-points), representing its vertices.
In addition, it provides ordered lists indicating which points appear along each edge of `meshA` and `meshB`.
Optionally, an explicit mesh of the common subdivision can be constructed via the `constructMesh()` function.

**Example:** construct an intrinsic triangulation, an build its common subdivision as a mesh.

```cpp
#include "geometrycentral/surface/common_subdivision.h"
#include "geometrycentral/surface/intrinsic_triangulation.h"
#include "geometrycentral/surface/integer_coordinates_intrinsic_triangulation.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

ManifoldSurfaceMesh& mesh;
VertexPositionGeometry& geometry; // your mesh and geometry

// Build an interesting intrinsic triangulation
IntegerCoordinatesIntrinsicTriangulation intTri(mesh, geometry);
intTri.delaunayRefine();

// Get the common subdivision
CommonSubdivision& cs = intTri.getCommonSubdivision();

// Build a mesh of the common subdivision, and interpolate vertex positions for it
cs.constructMesh();
ManifoldSurfaceMesh& csMesh = *cs.mesh;
VertexData<Vector3> csPositions = cs.interpolateAcrossA(geometry->vertexPositions);
```

## Members

??? func "`#!cpp ManifoldSurfaceMesh& CommonSubdivision::meshA`<br/>`#!cpp ManifoldSurfaceMesh& CommonSubdivision::meshB`"

    The two meshes whose common subdivision this object encodes.

??? func "`#!cpp std::deque<CommonSubdivisionPoint> CommonSubdivision::subdivisionPoints`"

    A list of the vertices (subdivision points) of the common subdivision. See documentation of `CommonSubdivisionPoint` below.

??? func "`#!cpp EdgeData<std::vector<CommonSubdivisionPoint*>> CommonSubdivision::pointsAlongA`<br/>`#!cpp EdgeData<std::vector<CommonSubdivisionPoint*>> CommonSubdivision::pointsAlongB`"

    Stores the list of the vertices of the common subdivision along each edge of `meshA` (resp. `meshB`), represented as pointers to elements of `subdivisionPoints`. These lists are ordered according to the edge's orientation and include the start and end points of the edges.
    
    If edge `eA` of `meshA` and `eB` of `meshB` run parallel to each other, then `pointsAlongA[eA]` will store an "intersection" of type `CSIntersectionType::EDGE_PARALLEL` to record the fact that it runs parallel to edge `eB`, and vice versa.

## Queries and Accessors

??? func "`#!cpp size_t CommonSubdivision::nVertices() const`"

    Counts the number of vertices of the common subdivision.
    Works even if an explicit mesh of the common subdivision has not been constructed yet.
    
??? func "`#!cpp std::tuple<size_t, size_t, size_t> CommonSubdivision::elementCounts() const`"

    Counts the number of vertices, edges, and faces of the common subdivision.
    Works even if an explicit mesh of the common subdivision has not been constructed yet.

??? func "`#!cpp size_t CommonSubdivision::intersectionsA(Edge eA) const`<br/>`#!cpp size_t CommonSubdivison::intersectionsB(Edge eB) const`"

    Counts the number of common subdivision vertices (not including endpoints) along edge `eA` of `meshA` (resp. `eB` of `meshB`). These intersection counts can be thought of as _normal coordinates_.
    

??? func "`#!cpp std::unique_ptr<SimplePolygonMesh> CommonSubdivision::buildSimpleMesh()`"
  
    Construct and return a `SimplePolygonMesh` of the common subdivision.
    
    In cases where there are some errors in intersection data defining the common subvidision, it may not be possible to construct the `ManifoldSurfaceMesh` with `constructMesh()`, and this is the only option.

## Common Subdivision Points
A `CommonSubdivisionPoint` is a struct representing a vertex of the common subdivision. Each vertex of the common subdivision is the intersection of a mesh element from `meshA` with a mesh element of `meshB`.

??? func "`#!cpp CSIntersectionType CommonSubdivisionPoint::intersectionType`"
    The type of elements intersecting (e.g. edge-edge intersection, or face-vertex intersection).

??? func "`#!cpp SurfacePoint CommonSubdivisionPoint::posA`<br/>`#!cpp SurfacePoint CommonSubdivisionPoint::posB`"
    The location of this intersection on `meshA` (resp. `meshB`), represented as a [SurfacePoint](/surface/utilities/surface_point).

??? func "`#!cpp bool CommonSubdivisionPoint::orientation`"
    For edge-edge intersections, this stores the orientation of the intersection. It is set to `true` for positive intersections, meaning that when traveling along the edge of `meshA`, the intersecting edge of `meshB` points to the left.

### Intersection Types
```cpp
enum class CSIntersectionType {
  VERTEX_VERTEX,
  EDGE_TRANSVERSE,
  EDGE_PARALLEL,
  FACE_VERTEX, // Face of mesh A, Vertex of mesh B
  EDGE_VERTEX  // Edge of mesh A, Vertex of mesh B
};
```

???+ info "'Parallel' intersections"
    In addition to storing points representing transverse intersections of mesh elements, we also store points which represent edges of `meshA`  and `meshB` that run parallel to each other. These points are not really vertices of the common subdivision, but provide a convenient way of encoding the fact that edges run parallel to each other. These points are tagged with type `CSIntersectionType::EDGE_PARALLEL`.


## Mesh connectivity

Optionally, the raw intersections stored in the common subdivision can be used to explicitly construct a `SurfaceMesh` (with all the usual halfedge connectivity) of the common subdivision; this mesh can then be used for many higher-level geometric operations.

**The routines and members in this section all require that constructMesh() has been called first.**

Note that in cases where there are some errors in intersection data defining the common subdivision, it may not be possible to construct a manifold mesh of the common subdivision. `constructMesh()` will fail with an exception in this case. Note that `buildSimpleMesh()` above can be used to build a plain old vertex-face adjacency list representation of the mesh, although it cannot be used for the various routines in this section.

By default, the common subdivision mesh is connectivity-only; there is no geometry associated with. Geometry can be recovered by either interpolating vertex positions from a mesh sitting in space (like `interpolateAcrossA(geometry->vertexPositions)`), or intrinsically via (like `interpolationEdgeLengthsA()`).

### Mesh Members

??? func "`#!cpp std::unique_ptr<ManifoldSurfaceMesh> CommonSubdivision::mesh`"

    An explicit mesh of the common subdivision.

??? func "`#!cpp VertexData<CommonSubdivisionPoint*> CommonSubdivision::sourcePoints`"

    Defined per-vertex of the mesh `CommonSubdivision::mesh`, associates each vertex of `mesh` with the corresponding element of `subdivisionPoints`.

??? func "`#!cpp FaceData<Face> CommonSubdivision::sourceFaceA` <br/> `#!cpp FaceData<Face> CommonSubdivision::sourceFaceB`"

    Defined per-face of the mesh `CommonSubdivision::mesh`, associates each face of `mesh` with the corresponding face of `meshA` (resp. `meshB`) which contains it.

### Mesh Constructors

??? func "`#!cpp void CommonSubdivision::constructMesh(bool triangulate = true, bool skipIfAlreadyConstructed = true)`"

    Compute an explicit mesh of the common subdivision, storing it in the `mesh` member.

    Initially, the faces of the common subdivision are polygons. If `triangulate` is `true`, then these faces are immediately triangulated internally.

    If `skipIfAlreadyConstructed` is `true`, this function does nothing when called multiple times. Otherwise, it deletes and reconstructs on subsequent calls.
    
??? func "`#!cpp void CommonSubdivision::triangulateMesh() const`"

    Triangulate the `mesh` member. Internally, the `sourceFaceA` and `sourceFaceB` are also updated to reflect the triangulation.
    
    This is already called be default in `constructMesh()` if `triangulate=true`.


### Mesh Utilities

Interpolate values defined at vertices from either `meshA` or `meshB` to the vertices of the common subdivision mesh.

??? func "`#!cpp VertexData<T> CommonSubdivision::interpolateAcrossA(const VertexData<T>& dataA)`<br/>`#!cpp VertexData<T> CommonSubdivision::interpolateAcrossB(const VertexData<T>& dataB)`"

    Linearly interpolates data on `meshA` (resp. `meshB`) to the common subdivision.


??? func "`#!cpp FaceData<T> CommonSubdivision::copyFromA(const FaceData<T>& dataA)`<br/>`#!cpp FaceData<T> CommonSubdivision::copyFromB(const FaceData<T>& dataB)`"

    Copy data at faces from one of the meshes to the common subdivision. Each face of the common subdivision gets the value from the face which contains it. The return value is defined per-face of the common subdivision mesh.

??? func "`#!cpp SparseMatrix<double> CommonSubdivison::interpolationMatrixA()`<br/>`#!cpp SparseMatrix<double> CommonSubdivision::interpolationMatrixB()`"

    Yields a `|V| x |V_A|` matrix (resp. `|V| x |V_B|`) which linearly interpolates data on `meshA` (resp. `meshB`) to the common subdivision. Here `|V|` denotes the number of vertices in the common subdivision.

    
??? func "`#!cpp EdgeData<double> CommonSubdivision::interpolateEdgeLengthsA(const EdgeData<double>& lengthA)`<br/>`#!cpp EdgeData<double> CommonSubdivision::interpolateEdgeLengthsB(const EdgeData<double>& lengthB)`"

    Takes in edge lengths for `meshA` (resp. `meshB`) and computes the edge lengths of the common subdivision.
    
    Note that in the standard case of an intrinsic triangulation with Euclidean metric-preserving edge flips, calling either of these methods with the edge lengths from the respective triangulation will produce identical outputs (up to floating-point error).

??? func "`#!cpp SparseMatrix<double> CommonSubdivision::vertexGalerkinMassMatrixFromPositionsA(const VertexData<Vector3>& positionsA)`<br/>`#!cpp SparseMatrix<double> CommonSubdivision::vertexGalerkinMassMatrixFromPositionsB(const VertexData<Vector3>& positionsB)`"

    Takes in vertex positions for `meshA` (resp. `meshB`) and computes the Galerkin mass matrix of the common subdivision.

??? func "`#!cpp SparseMatrix<double> CommonSubdivision::vertexGalerkinMassMatrixFromLengthsA(const EdgeData<Vector3>& lengthsA)`<br/>`#!cpp SparseMatrix<double> CommonSubdivision::vertexGalerkinMassMatrixFromLengthsB(const EdgeData<Vector3>& lengthsB)`"

    Takes in edge lengths for `meshA` (resp. `meshB`) and computes the Galerkin mass matrix of the common subdivision.
    
    

### Additional Utilities

??? func "`#!cpp FaceData<double> niceColors(ManifoldSurfaceMesh& mesh, int kColors=7)`"

    Take `kColors` evenly spaced values on [0,1] and treat them as categorical labels. These labels are assigned to mesh faces in the sense of a graph coloring, with further heuristics to try to avoid neighbors-of-neighbors at vertices.

    This function is useful for generating a nice coloring of an intrinsic triangulation, defined along the common subdivision, for visualization.
