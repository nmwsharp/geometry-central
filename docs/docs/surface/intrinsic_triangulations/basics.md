Intrinsic triangulations begin with a simple idea: rather than representing a triangle mesh's geometry via vertex positions, instead store a length associated with each edge. This turns out to be just the right data to compute many geometric quantities, including areas, angles, Laplacians, geodesic distance, etc; from these one can evaluate many algorithms directly on an intrinsic triangulation.

The general workflow with intrinsic triangulations is to start out with some initial mesh having 3D vertex positions, and construct a new intrinsic triangulation which "sits on top" of the initial mesh, and can then be modified and improved.  The benefit of working with intrinsic triangulations is that they enable a family of simple & powerful mesh operations to robustly compute with a surface while preserving its shape exactly---these operations are not possible when working directly with vertex positions. For a general introduction to intrinsic triangulations and the technique used here, [see this course](https://nmwsharp.com/media/papers/int-tri-course/int_tri_course.pdf).

![intrinsic triangulation teaser](/media/int_tri_teaser.jpg)
The above image shows an intrinsic triangulation of a poorly-tessellated surface. The original mesh is given by the black wireframe, while the intrinsic triangulation (after Delaunay refinement) is given by the colored triangles. Computing with the intrinsic triangulation yields significant benefits, while still implicitly encoding the original shape exactly.

The tricky part of working with intrinsic triangulations like this is encoding the correspondence to the original mesh. Geometry-central contains two data structures for encoding this correspondence, which allow you to e.g. transfer data between to original mesh and the intrinsic triangulation, or map a point in barycentric coordinates between the triangulations. Additionally, these data structures offer routines for standard operations such as generating a high-quality intrinsic triangulation of a given surface.

The main interface for working with intrinsic triangulations is `IntrinsicTriangulation`. This interface is realized by the `SignpostIntrinsicTriangulation` and `IntegerCoordinatesIntrinsicTriangulation` classes (see [here](/surface/geometry/geometry/#geometry-hierarchy) for a brief discussion of interfaces, realization, and polymorphism in C++). The routines in this section show how to initialize a triangulation, improve its quality, and run algorithms on it. In particular, this intrinsic triangulation satisfies the [IntrinsicGeometryInterface](/surface/geometry/geometry/#intrinsic-geometry), and thus can be used directly as input to many algorithms in geometry-central.

[This repository](https://github.com/nmwsharp/navigating-intrinsic-triangulations-demo) contains a simple demo application with a GUI for manipulating intrinsic triangulations, and demonstrates much of the fuctionality documented here.

Headers:
```cpp
#include "geometrycentral/surface/intrinsic_triangulation.h"
#include "geometrycentral/surface/signpost_intrinsic_triangulation.h"
#include "geometrycentral/surface/integer_coordinates_intrinsic_triangulation.h"
```


## Example

```cpp
#include "geometrycentral/surface/signpost_intrinsic_triangulation.h"
#include "geometrycentral/surface/direction_fields.h"
#include "geometrycentral/surface/meshio.h"
using namespace geometrycentral;
using namespace surface;

// Load a mesh
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;
std::tie(mesh, geometry) = loadMesh(filename);

// Construct an intrinsic triangulation
// (here, using the signpost data structure)
std::unique_ptr<IntrinsicTriangulation> intTri(
  new SignpostIntrinsicTriangulation(*mesh, *geometry));


// == Improve the quality

// flip edges to get the intrinsic Delaunay triangulation
intTri->flipToDelaunay(); 

// optional: go even further and insert new vertices guarantee 
// a minimum angle bound (default: 25 deg) in the underlying 
// intrinsic triangulation
intTri->delaunayRefine();


// == Do some computation on the intrinsic triangulation 

// Compute a smooth direction field
VertexData<Vector2> directionsOnInt = computeSmoothestVertexDirectionField(*intTri);

// Copy the result back to the vertices of the input mesh
VertexData<Vector2> directionsOnOrig = intTri->restrictToInput(directionsOnInt);
```

## Intrinsic Triangulation

This class tracks the geometry and connectivity of an intrinsic triangulation sitting atop a triangle mesh. It "is-a" `IntrinsicGeometryInterface`, so it can be used directly as a geometry input to many geometry-central subroutines. Furthermore, we can compute [all sorts of geometric quantities](/surface/geometry/quantities) on it just as we normally would with a geometry object.

`IntrinsicTriangulation` is an interface and cannot be constructed directly: to use it, you must construct a `SignpostIntrinsicTriangulation` or `IntegerCoordinatesIntrinsicTriangulation`.

Example: compute some quantities
```cpp
// Construct an intrinsic triangulation
// (here, using the signpost data structure)
std::unique_ptr<IntrinsicTriangulation> intTri(
  new SignpostIntrinsicTriangulation(*mesh, *geometry));

intTri->flipToDelaunay(); 

// edge lengths
intTri->requireEdgeLengths();
for(Edge e : intTri->intrinsicMesh.edges()) {
  // print them, as an example
  std::cout << "Length of edge " << e << " is " << intTri->edgeLengths[e] << std::endl;
}

// vertex dual area (mass)
intTri->requireVertexDualAreas();
for(Vertex v : intTri->intrinsicMesh.vertices()) {
  std::cout << "Area of vertex " << v << " is " << intTri->vertexDualAreas[v] << std::endl;
}

// Laplace matrix
// (this is this intrinsic Delaunay Laplace matrix, since we flipped to Delaunay)
intTri->requireCotanLaplacian();
SparseMatrix<double> L = intTri->cotanLaplacian; // an Eigen sparse matrix
```


## Basic API
Note that some additional functions and members can be found in the `IntrinsicTriangulation` class at 
`geometrycentral/surface/intrinsic_triangulation.h`.

### Members

??? func "`#!cpp ManifoldSurfaceMesh& inputMesh`"
  
    The underlying mesh which the intrinsic triangulation sits on top of. It must not be modified during the intrinsic triangulation's lifetime.


??? func "`#!cpp IntrinsicGeometryInterface& inputGeom`"
  
    The underlying geometry corresponding to `inputMesh`. It must not be modified during the intrinsic triangulation's lifetime.


??? func "`#!cpp std::unique_ptr<ManifoldSurfaceMesh>& intrinsicMesh`"

    A triangle mesh (distinct from `inputMesh`), which represents the connectivity of the intrinsic triangulation. Right after construction, it will be a copy of `inputMesh`, but may be modified by subsequent mutations like edge flips.

    Note that that because the `IntrinsicTriangulation` is-a `IntrinsicGeometryInterface`, there is also a  `IntrinsicTriangulation::mesh` member reference which refers to this same mesh object, inherited from the geometry interface.


??? func "`#!cpp VertexData<SurfacePoint>& vertexLocations`"
  
    A position on `inputMesh` for each vertex in the intrinsic triangulation. Initially, these positions will all corresponding vertices on the input mesh, but as e.g. refinements are performed, some intrinsic vertices will have positions inside faces of the original triangulation.
    

### Queries and Accessors

??? func "`#!cpp EdgeData<std::vector<SurfacePoint>> IntrinsicTriangulation::traceAllIntrinsicEdgesAlongInput()`"
  
    Traces out the path that each intrinsic edge takes along the surface, as a sequence of [SurfacePoint](/surface/utilities/surface_point)s beginning with the vertex at the tail of the halfedge and ending with the vertex at the tip of the halfedge.
    
??? func "`#!cpp std::vector<SurfacePoint> IntrinsicTriangulation::traceIntrinsicHalfedgeAlongInput(Halfedge intrinsicHe)`"

    Traces out the path that an intrinsic edge takes along the surface, as a sequence of [SurfacePoint](/surface/utilities/surface_point)s beginning with the vertex at the tail of the halfedge and ending with the vertex at the tip of the halfedge.
    
    Nore: When using an `IntegerCoordinatesIntrinsicTriangulation`, calling `traceEdges()` is significantly more efficient than calling `traceHalfedge(he)` on each edge.
    
??? func "`#!cpp EdgeData<std::vector<SurfacePoint>> IntrinsicTriangulation::traceAllInputEdgesAlongIntrinsic()`"
  
    Traces out the path that each input edge takes along the surface, as a sequence of [SurfacePoint](/surface/utilities/surface_point)s beginning with the vertex at the tail of the halfedge and ending with the vertex at the tip of the halfedge.
    
??? func "`#!cpp std::vector<SurfacePoint> IntrinsicTriangulation::traceInputHalfedgeAlongIntrinsic(Halfedge inputHe)`"

    Traces out the path that an input edge takes along the surface, as a sequence of [SurfacePoint](/surface/utilities/surface_point)s beginning with the vertex at the tail of the halfedge and ending with the vertex at the tip of the halfedge.

??? func "`#!cpp CommonSubdivision& IntrinsicTriangulation::getCommonSubdivision()`"
  
    Returns the [common subdivison](/surface/intrinsic_triangulations/common_subdivision) of the input and intrisnic meshes. May construct it from scratch if this is the first time it is needed. The intrinsic triangulation manages the lifetime of the subdivision---it will be deallocated if
    
    (a) this object is deleted, or
    
    (b) the triangulation is mutated, invalidating the common subdivision.
    
    Be sure to copy it if you want to retain it through those operations.
    
??? func "`#!cpp SurfacePoint IntrinsicTriangulation::equivalentPointOnIntrinsic(const SurfacePoint& pointOnInput)`"
  
    Given a point on the input triangulation, returns the corresponding point on the current intrinsic triangulation.
    
    The input and output are given as a [SurfacePoint](/surface/utilities/surface_point), which may be a vertex, point along an edge, or point within some face.
    
??? func "`#!cpp SurfacePoint IntrinsicTriangulation::equivalentPointOnInput(const SurfacePoint& pointOnIntrinsic)`"
  
    Given a point on the current intrinsic triangulation, returns the corresponding point on the input triangulation.

    The input and output are given as a [SurfacePoint](/surface/utilities/surface_point), which may be a vertex, point along an edge, or point within some face.
    
??? func "`#!cpp VertexData<T> IntrinsicTriangulation::sampleFromInput(const VertexData<T>& dataOnInput)`"
  
    Given data defined on the vertices of the input triangulation, samples it to the vertices of the intrinsic triangulation.
    
??? func "`#!cpp VertexData<T> IntrinsicTriangulation::restrictToInput(const VertexData<T>& dataOnIntrinsic)`"
  
    Given data defined on the vertices of the intrinsic triangulation, restrict it to the vertices of the input triangulation.
    
??? func "`#!cpp bool IntrinsicTriangulation::isDelaunay()`"

    Returns true if all edges in the intrinsic triangulation satisfy the _intrinsic Delaunay property_ (aka have nonnegative cotangent weights).

??? func "`#!cpp double IntrinsicTriangulation::minAngleDegrees()`"

    Returns the smallest corner angle in the intrinsic triangulation, in degrees.


    
### High-Level Mutators

??? func "`#!cpp void IntrinsicTriangulation::flipToDelaunay()`"
  
    Flips edges in the intrinsic triangulation until is satisfies the intrinsic Delaunay criterion.
    
??? func "`#!cpp void IntrinsicTriangulation::delaunayRefine(double angleThreshDegrees = 25, double circumradiusThresh = inf, size_t maxInsertions = inf)`"

    Applies Chew's 2nd algorithm to the intrinsic triangulation, flipping edges and inserting vertices until the triangulation simultaneously:
    
     - satisfies the intrinsic Delaunay criterion
     
     - has no angles smaller than `angleThreshDegrees` (values > 30 degrees may not terminate)
     
     - has no triangles larger than `circumradiusThresh`
     
    Terminates no matter what after `maxInsertions` insertions (infinite by default)

    The algorithm converges with angle threshold settings up to 30 degrees (away from ultra-skinny needle vertices and boundary angles which cannot be improved).
  
    
??? func "`#!cpp void IntrinsicTriangulation::delaunayRefine(cosnt std::function<bool(Face)>& shouldRefine, size_t maxInsertions = inf)`"
    
    General version of intrinsic Delaunay refinement, taking a function which will be called to determine if a triangle should be refined. Will return only when all triangles pass this function, or `maxInsertions` is exceeded, so be sure to chose arguments such that the function terminates.
    
??? func "`#!cpp void IntrinsicTriangulation::setMarkedEdges(const EdgeData<bool>& markedEdges)`"

    Set a subset of the edges which are special marked edges. If true, the edge is fixed, and will not be flipped (e.g. in `flipToDelaunay()`. 
 
    The class automatically internally hanldes updates to this array as edge splits are performed, so if a marked edge is split the two resulting edges will be marked.
    
### Low-Level Mutators


!!! note "Refreshing quantities"

    The intrinsic triangulation is-a [geometry](/surface/geometry/geometry/) object, which means that one may `require()` quantities from it. However, for efficiency reasons, these quantities are not automatically updated after each low-level mutataion. Call `refreshQuantities()` after a sequence of mutations to update.
  

??? func "`#!cpp bool IntrinsicTriangulation::flipEdgeIfNotDelaunay(Edge e)`"

    If the edge is not (intrinsic) Delaunay, flip it. Returns true if flipped.
    
??? func "`#!cpp bool IntrinsicTriangulation::flipEdgeIfPossible(Edge e)`"

    If the edge can be flipped, flip it. Returns true if flipped.

    Edges cannot be flipped if:
    - the are boundary edges
    - either endpoint is a degree 1 vertex
    - when laid out in the plane, the diamond containing the edge is a nonconvex quadrilateral
    
??? func "`#!cpp Vertex IntrinsicTriangulation::insertVertex(SurfacePoint newPositionOnIntrinsic)`"

    Inserts a new vertex in the intrinsic triangulation. The adjacent faces will be triangulated to maintain a triangle mesh. In degenerate situations, this procedure may fail and return `Vertex()`.
    
??? func "`#!cpp Vertex IntrinsicTriangulation::insertCircumcenter(Face f)`"

    Insert an new vertex in to the intrinsic triangulation, at the intrinsic circumcenter of some face. Note that for an obtuse face, the circumcenter is outside of the triangle; the resulting location will be found by tracing. Adjacent faces will be triangulated to maintain a triangle mesh. In degenerate situations, this procedure may fail and return `Vertex()`.
    
??? func "`#!cpp Vertex IntrinsicTriangulation::insertBarycenter(Face f)`"

    Insert a new vertex in the intrinsic triangulation, at the intrinsic barycenter of face `f`.
    
??? func "`#!cpp Face IntrinsicTriangulation::removeInsertedVertex(Vertex v)`"

    Removes an inserted vertex from the triangulation. Returns `Face()` if the vertex cannot be removed.
    
??? func "`#!cpp Halfedge IntrinsicTriangulation::splitEdge(Halfedge he, double tSplit)`"

    Splits an oriented edge. Returns the halfedge which starts at the inserted vertex and points in the same direction as the input halfedge.
    
Additionally, callbacks can be registered by inserting them in to the following lists. These callbacks will be invoked whenever the corresponding mesh operation is performed within the intrinsic triangulation. This can be useful e.g. for maintaining values as the mesh is modified.

??? func "`#!cpp std::list<std::function<void(Edge)>> edgeFlipCallbackList`"

    Member variable, a list of callbacks which are invoked whenever an edge is flipped.

    Example: construct a lambda callback and register it
    ```cpp
    SignpostIntrinsicTriangulation intTri = /* constructed somehow */;

    auto updateOnFlip = [&](Edge flipE) {
      // do something
    };
    auto callbackRef = intTri->edgeFlipCallbackList.insert(std::end(intTri->edgeFlipCallbackList), updateOnFlip);
    ```

??? func "`#!cpp std::list<std::function<void(Face, Vertex)>> faceInsertionCallbackList`"
    
    Member variable, a list of callbacks which are invoked whenever vertex is inserted in to a face. The first argument is the old face, and the second is the new vertex.  

    Example: construct a lambda callback and register it
    ```cpp
    SignpostIntrinsicTriangulation intTri = /* constructed somehow */;

    auto updateOnInsert = [&](Face oldF, Vertex newV) {
      // do something
    };
    auto callbackRef = intTri->faceInsertionCallbackList.insert(std::end(intTri->faceInsertionCallbackList), updateOnInsert);
    ```


??? func "`#!cpp std::list<std::function<void(Edge, Halfedge, Halfedge)>> edgeSplitCallbackList`"
    
    Member variable, a list of callbacks which are invoked whenever an edge is split. The first edge is the old edge which was just split, and the two halfedges are along the new edges which were created by the split.

    Example: construct a lambda callback and register it
    ```cpp
    SignpostIntrinsicTriangulation intTri = /* constructed somehow */;

    auto updateOnSplit = [&](Edge oldE, Halfedge newHe1, Halfedge newHe2) {
      // do something
    };
    auto callbackRef = intTri->edgeSplitCallbackList.insert(std::end(intTri->edgeSplitCallbackList), updateOnSplit);
    ```

## Citations

The above data structures are described in the following works:

```bib
@article{Sharp:2021:GPI,
  author = {Sharp, Nicholas and Gillespie, Mark and Crane, Keenan},
  title = {Geometry Processing with Intrinsic Triangulations},
  booktitle = {ACM SIGGRAPH 2021 courses},
  series = {SIGGRAPH '21},
  year = {2021},
  publisher = {ACM},
  address = {New York, NY, USA},
}
```

```bib
@article{sharp2019navigating,
  title={Navigating intrinsic triangulations},
  author={Sharp, Nicholas and Soliman, Yousuf and Crane, Keenan},
  journal={ACM Transactions on Graphics (TOG)},
  volume={38},
  number={4},
  pages={55},
  year={2019},
  publisher={ACM}
}
```

```bib
@article{gillespie2021integer,
  title={Integer Coordinates for Intrinsic Geometry Processing},
  author={Gillespie, Mark and Sharp, Nicholas and Crane, Keenan},
  journal={arXiv preprint arXiv:2106.00220},
  year={2021}
}
```
