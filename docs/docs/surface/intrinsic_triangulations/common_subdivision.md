`#include "geometrycentral/surface/common_subdivision.h"`

Intuitively, the common subdivision of meshes $T_A$ and $T_B$ is the polygon mesh that you obtain if you cut $T_B$ along the edges of $T_A$. The vertices of the common subdivision correspond to intersections of simplices of $T_A$ with simplices of $T_B$ (e.g. intersections between edges). 
![common subdivision of two meshes](/media/common_subdivision.svg){: style="max-height: 10em; display: block; margin-left: auto; margin-right: auto;"}

The `CommonSubdivision` class represents the common subdivision of two meshes, referred to as `meshA` and `meshB`.
??? warning "Every vertex of `meshA` must also be a vertex of `meshB`."

The common subdivision is encoded by a list of [CommonSubdivisionPoints](#common-subdivision-points), representing its vertices.
In addition, it provides ordered lists indicating which points appear along each edge of `meshA` and `meshB`.
Optionally, an explicit mesh of the common subdivision can be constructed via the `constructMesh()` function.

## Members

??? func "`#!cpp ManifoldSurfaceMesh& CommonSubidivision::meshA`<br/>`#!cpp ManifoldSurfaceMesh& CommonSubdivision::meshB`"

    The two meshes whose common subdivision this object encodes.

??? func "`#!cpp std::deque<CommonSubdivisionPoint> CommonSubdivision::subdivisionPoints`"

    A list of the vertices of the common subdivision.

??? func "`#!cpp EdgeData<std::vector<CommonSubdivisionPoint*>> CommonSubdivision::pointsAlongA`<br/>`#!cpp EdgeData<std::vector<CommonSubdivisionPoint*>> CommonSubdivision::pointsAlongB`"

    Stores the list of the vertices of the common subdivision along each edge of `meshA` (resp. `meshB`), represented as pointers to elements of `subdivisionPoints`. These lists are ordered according to the edge's orientation and include the start and end points of the edges.
    
    If edge `eA` of `meshA` and `eB` of `meshB` run parallel to each other, then `pointsAlongA[eA]` will store an "intersection" of type `CSIntersectionType::EDGE_PARALLEL` to record the fact that it runs parallel to edge `eB`, and vice versa.

??? func "`#!cpp std::unique_ptr<ManifoldSurfaceMesh> CommonSubdivision::mesh`"

    An explicit mesh of the common subdivision.

??? func "`#!cpp VertexData<CommonSubdivisionPoint*> CommonSubdivision::sourcePoints`"

    Associates each vertex of `mesh` with the corresponding element of `subdivisionPoints`.

## Queries and Accessors

??? func "`#!cpp VertexData<T> CommonSubdivision::interpolateAcrossA(const VertexData<T>& dataA)`<br/>`#!cpp VertexData<T> CommonSubdivision::interpolateAcrossB(const VertexData<T>& dataB)`"

    Linearly interpolates data on `meshA` (resp. `meshB`) to the common subdivision.

??? func "`#!cpp SparseMatrix<double> CommonSubdivison::interpolationMatrixA()`<br/>`#!cpp SparseMatrix<double> CommonSubdivision::interpolationMatrixB()`"

    Yields a `|V| x |V_A|` matrix (resp. `|V| x |V_B|`) which linearly interpolates data on `meshA` (resp. `meshB`) to the common subdivision. Here `|V|` denotes the number of vertices in the common subdivision.
    
??? func "`#!cpp EdgeData<double> CommonSubdivision::interpolationEdgeLengthsA(const EdgeData<double>& lengthA)`<br/>`#!cpp EdgeData<double> CommonSubdivision::interpolationEdgeLengthsB(const EdgeData<double>& lengthB)`"

    Takes in edge lengths for `meshA` (resp. `meshB`) and computes the edge lengths of the common subdivision.

??? func "`#!cpp SparseMatrix<double> CommonSubdivision::vertexGalerkinMassMatrixFromPositionsA(const VertexData<Vector3>& positionsA)`<br/>`#!cpp SparseMatrix<double> CommonSubdivision::vertexGalerkinMassMatrixFromPositionsB(const VertexData<Vector3>& positionsB)`"

    Takes in vertex positions for `meshA` (resp. `meshB`) and computes the Galerkin mass matrix of the common subdivision.

??? func "`#!cpp SparseMatrix<double> CommonSubdivision::vertexGalerkinMassMatrixFromLengthsA(const EdgeData<Vector3>& lengthsA)`<br/>`#!cpp SparseMatrix<double> CommonSubdivision::vertexGalerkinMassMatrixFromLengthsB(const EdgeData<Vector3>& lengthsB)`"

    Takes in edge lengths for `meshA` (resp. `meshB`) and computes the Galerkin mass matrix of the common subdivision.
    
??? func "`#!cpp size_t CommonSubdivision::nVertices() const`"

    Counts the number of vertices of the common subdivision.
    Works even if an explicit mesh of the common subdivision has not been constructed yet.
    
??? func "`#!cpp std::tuple<size_t, size_t, size_t> CommonSubdivision::elementCounts() const`"

    Counts the number of vertices, edges, and faces of the common subdivision.
    Works even if an explicit mesh of the common subdivision has not been constructed yet.
    
??? func "`#!cpp size_t CommonSubdivision::intersectionA(Edge eA) const`<br/>`#!cpp size_t CommonSubdivison::intersectionB(Edge eB) const`"

    Counts the number of common subdivision vertices (not including endpoints) along edge `eA` of `meshA` (resp. `eB` of `meshB`). These intersection counts can be thought of as _normal coordinates_.
    
## Mutators
??? func "`#!cpp void CommonSubdivision::constructMesh() const`"

    Compute an explicit mesh of the common subdivision, storing it in the `mesh` member.
    
??? func "`#!cpp void CommonSubdivision::triangulateMesh() const`"

    Triangulate the `mesh` member (by default, `mesh` is a polygon mesh).
   
## Misc
??? func "`#!cpp FaceData<double> niceColors(ManifoldSurfaceMesh& mesh, int kColors=7)`"

    Take `kColors` evenly spaced values on [0,1] and treat them as categorical labels. These labels are assigned to mesh faces in the sense of a graph coloring, with further heuristics to try to avoid neighbors-of-neighbors at vertices.

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

??? info "'Parallel' intersections"
    In addition to storing points representing transverse intersections of mesh elements, we also store points which represent edges of `meshA`  and `meshB` that run parallel to each other. These points are not really vertices of the common subdivision, but provide a convenient way of encoding the fact that edges run parallel to each other. These points are tagged with type `CSIntersectionType::EDGE_PARALLEL`.
