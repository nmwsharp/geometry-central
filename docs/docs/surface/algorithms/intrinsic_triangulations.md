Intrinsic triangulations begin with a simple idea: rather than representing a triangle mesh's geometry via vertex positions, instead store a length associated with each edge. This turns out to be just the right data to compute many geometric quantities, including areas, angles, Laplacians, geodesic distance, etc; from these one can evaluate many algorithms directly on an intrinsic triangulation.

The general workflow with intrinsic triangulations is to start out with some initial mesh having 3D vertex positions, and construct a new intrinsic triangulation which "sits on top" of the initial mesh, and can then be modified and improved.  The benefit of working with intrinsic triangulations is that they enable a family of simple & powerful mesh operations to robustly compute with a surface while preserving its shape exactly---these operations are not possible when working directly with vertex positions. For a general introduction to intrinsic triangulations and the technique used here, [see this paper](http://www.cs.cmu.edu/~kmcrane/Projects/NavigatingIntrinsicTriangulations/paper.pdf).

![intrinsic triangulation teaser](/media/int_tri_teaser.jpg)
The above image shows an intrinsic triangulation of a poorly-tesselated surface. The original mesh is given by the black wireframe, while the intrinsic triangulation (after Delaunay refinement) is given by the colored triangles. Computing with the intrinsic triangulation yields significant benefits, while still implicitly encoding the original shape exactly.

The main class for working with intrinsic triangulations is the `SignpostIntrinsicTriangulation`. The routines in this section show how to initialize the triangulation, improve its quality, and run algorithms on it. In particular, this intrinsic triangulation satisfies the [IntrinsicGeometryInterface](/surface/geometry/geometry/#intrinsic-geometry), and thus can be used directly as input to many algorithms in geometry-central.

[This repository](https://github.com/nmwsharp/navigating-intrinsic-triangulations-demo) contains a simple demo application with a GUI for manipulating intrinsic triangulations, and demonstrates much of the fuctionality documented here.

`#include "geometrycentral/surface/signpost_intrinsic_triangulation.h"`


## Example

```cpp
#include "geometrycentral/surface/signpost_intrinsic_triangulation.h"
#include "geometrycentral/surface/direction_fields.h"
#include "geometrycentral/surface/meshio.h"
using namespace geometrycentral;
using namespace surface;

// Load a mesh
std::unique_ptr<HalfedgeMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;
std::tie(mesh, geometry) = loadMesh(filename);

// Construct an intrinsic triangulation
SignpostIntrinsicTriangulation intTri(*mesh, *geometry);


// == Improve the quality

// flip edges to get the intrinsic Delaunay triangulation
intTri.flipToDelaunay(); 

// optional: go even further and insert new vertices guarantee 
// a minimum angle bound (default: 25 deg) in the underlying 
// intrinsic triangulation
intTri.delaunayRefine();


// == Do some computation on the intrinsic triangulation 

// (example here: compute a smooth direction field)
VertexData<Vector2> directionsOnInt = computeSmoothestVertexDirectionField(intTri);

// copy the result back to the vertices of the input mesh
VertexData<Vector2> directionsOnOrig = intTri.restrictToInput(directionsOnInt);
```

## Signpost Intrinsic Triangulation

This class tracks the geometry and connectivity of an intrinsic triangulation sitting atop a triangle mesh. It "is-a" `IntrinsicGeometryInterface`, so it can be used directly as a geometry input to many geometry-central subroutines.

Note that some additional functions and members can be found in the `SignpostIntrinsicTriangulation` class at 
`geometrycentral/surface/signpost_intrinsic_triangulation.h`.

### Members

??? func "`#!cpp ManifoldSurfaceMesh& inputMesh`"
  
    The underlying mesh which the intrinsic triangulation sits on top of. It must not be modified during the intrinsic triangulation's lifetime.


??? func "`#!cpp IntrinsicGeometryInterface& inputGeom`"
  
    The underlying geometry corresponding to `inputMesh`. It must not be modified during the intrinsic triangulation's lifetime.


??? func "`#!cpp std::unique_ptr<ManifoldSurfaceMesh>& intrinsicMesh`"

    A triangle mesh (distinct from `inputMesh`), which represents the connectivity of the intrinsic triangulation. Right after construction, it will be a copy of `inputMesh`, but may be modified by subsequent mutations like edge flips.

    Note that that because the `SignpostIntrinsicTriangulation` is-a `IntrinsicGeometryInterface`, there is also a  `SignpostIntrinsicTriangulation::mesh` member reference which refers to this same mesh object, inherited from the geometry interface.


??? func "`#!cpp VertexData<SurfacePoint>& vertexLocations`"
  
    A position on `inputMesh` for each vertex in the intrinsic triangulation. Initially, these positions will all corresponding vertices on the input mesh, but as e.g. refinements are performed, some intrinsic vertices will have positions inside faces of the original triangulation.


Remember, because the `SignpostIntrinsicTriangulation` is-a `IntrinsicGeometryInterface`, we can compute [all sorts of geometric quantities](/surface/geometry/quantities) on it just as we normally would with a geometry object.

Example: compute some quantities
```cpp
SignpostIntrinsicTriangulation intTri(*mesh, *geometry);
intTri.flipToDelaunay(); 

// edge lengths
intTri.requireEdgeLengths();
for(Edge e : intTri->intrinsicMesh.edges()) {
  // print them, as an example
  std::cout << "Length of edge " << e << " is " << intTri.edgeLengths[e] << std::endl;
}

// vertex dual area (mass)
intTri.requireVertexDualAreas();
for(Vertex v : intTri->intrinsicMesh.vertices()) {
  std::cout << "Area of vertex " << v << " is " << intTri.vertexDualAreas[v] << std::endl;
}

// Laplace matrix
// (this is this intrinsic Delaunay Laplace matrix, since we flipped to Delaunay)
intTri.requireCotanLaplacian();
SparseMatrix<double> L = intTri.cotanLaplacian; // an Eigen sparse matrix
```


### Constructors

??? func "`#!cpp SignpostIntrinsicTriangulation::SignpostIntrinsicTriangulation(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& inputGeom)`"

    Initialize an intrinsic triangulation sitting on top of `mesh`. Recall that `IntrinsicGeometryInterface` can be almost any geometry object, including a `VertexPositionGeometry`.

    Initially, the intrinsic triangulation will be identical to the input mesh; it can be modified with the routines below.

### Queries 

??? func "`#!cpp SurfacePoint SignpostIntrinsicTriangulation::equivalentPointOnIntrinsic(const SurfacePoint& pointOnInput)`"

    Given a point on the input triangulation, returns the corresponding point on the current intrinsic triangulation.

    The input and output are given as a [SurfacePoint](/surface/utilities/surface_point), which may be a vertex, point along an edge, or point within some face.

??? func "`#!cpp SurfacePoint SignpostIntrinsicTriangulation::equivalentPointOnInput(const SurfacePoint& pointOnIntrinsic)`"

    Given a point on the current intrinsic triangulation, returns the corresponding point on the input triangulation.

    The input and output are given as a [SurfacePoint](/surface/utilities/surface_point), which may be a vertex, point along an edge, or point within some face.

??? func "`#!cpp std::vector<SurfacePoint> SignpostIntrinsicTriangulation::traceHalfedge(Halfedge he, bool trimEnd=true)`"

    Traces out the path an intrinsic edge takes along the surface, as a sequence of [SurfacePoint](/surface/utilities/surface_point) beginning with the vertex at the tail of the halfedge and ending with the vertex at the tip of the halfedge.

    - `he` the halfedge to trace
    - `bool trimEnd` if true, a simple heuristc is applied to clean up tiny fragments of the trajectory which result from not-quite-perfectly hitting the target vertex in floating point


??? func "`#!cpp EdgeData<std::vector<SurfacePoint>> SignpostIntrinsicTriangulation::traceEdges()`"

    Calls `traceHalfedge()` (above) on each edge of intrinsic triangulation to extract all edge paths.


??? func "`#!cpp bool SignpostIntrinsicTriangulation::isDelaunay()`"

    Returns true if all edges in the intrinsic triangulation satisfy the _intrinsic Delaunay property_ (aka have nonnegative cotangent weights).

??? func "`#!cpp double SignpostIntrinsicTriangulation::minAngleDegrees()`"

    Returns the smallest corner angle in the intrinsic triangulation, in degrees.


### High-level mutators

??? func "`#!cpp void SignpostIntrinsicTriangulation::flipToDelaunay()`"

    Flips edges in the intrinsic triangulation to give it the _intrinsic Delaunay property_. 

??? func "`#!cpp void SignpostIntrinsicTriangulation::delaunayRefine(double angleThreshDegrees = 25., double circumradiusThresh = std::numeric_limits<double>::infinity(), size_t maxInsertions = INVALID_IND)`"

    Applies Chew's 2nd algorithm to the intrinsic triangulation, flipping edges and inserting vertices until the triangulation is _intrinsic Delaunay_, AND satisfies the requested angle boundd and triangle size threshold.

    The algorithm converges with angle threshold settings up to 30 degrees (away from ultra-skinny needle vertices and boundary angles which cannot be improved).

    - `angleThreshDegrees` Remove skinny angles. Inserts vertices until all intrinsic triangle corner angles are larger than this threshold. Should converge for settings 30 degrees (though for values very close to 30 floating point numerics may become an issue).
    - `circumradiusThresh` Refine large triangles if not infinity. Inserts vertices until all triangles have a circumradius smaller than this threshold.
    - `maxInsertions` insert at most this many vertices, overriding the other two thresholds (ignored if `INVALID_IND`).

??? func "`#!cpp void SignpostIntrinsicTriangulation::setMarkedEdges(const EdgeData<bool>& markedEdges)`"

    Set a subset of the edges which are special marked edges. If true, the edge is fixed, and will not be flipped (e.g. in `flipToDelaunay()`. 
 
    The class automatically internally hanldes updates to this array as edge splits are performed, so if a marked edge is split the two resulting edges will be marked.  


### Low-level mutators
  
??? func "`#!cpp bool SignpostIntrinsicTriangulation::flipEdgeIfNotDelaunay(Edge e)`"

    If the edge is not (intrinsic) Delaunay, flip it. Returns true if flipped.

    - `Edge e` Some edge in the intrinsic triangulation to flip.

??? func "`#!cpp bool SignpostIntrinsicTriangulation::flipEdgeIfPossible(Edge e)`"

    If the edge can be flipped, flip it. Returns true if flipped.

    Edges cannot be flipped if:
    - the are boundary edges
    - either endpoint is a degree 1 vertex
    - when laid out in the plane, the diamond containing the edge is a nonconvex quadrilateral

    Arguments:
    - `Edge e` Some edge in the intrinsic triangulation to flip.

??? func "`#!cpp Vertex SignpostIntrinsicTriangulation::insertVertex(SurfacePoint newPositionOnIntrinsic)`"

    Insert an new vertex in to the intrinsic triangulation, inside some face or along some edge. The adjacent faces will be triangulated to maintain a triangle mesh.

??? func "`#!cpp Vertex SignpostIntrinsicTriangulation::insertCircumcenter(Face f)`"

    Insert an new vertex in to the intrinsic triangulation, at the intrinsic circumcenter of some face. Note that for an obtuse face, the circumcenter is outside of the triangle; the resulting location will be found by tracing. Adjacent faces will be triangulated to maintain a triangle mesh.


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


??? func "`#!cpp std::list<std::function<void(Edge, Halfedge, Halfedge)>> SignpostIntrinsicTriangulation::edgeSplitCallbackList`"
    
    Member variable, a list of callbacks which are invoked whenever an edge is split. The first edge is the old edge which was just split, and the two halfedges are along the new edges which were created by the split.

    Example: construct a lambda callback and register it
    ```cpp
    SignpostIntrinsicTriangulation intTri = /* constructed somehow */;

    auto updateOnSplit = [&](Edge oldE, Halfedge newHe1, Halfedge newHe2) {
      // do something
    };
    auto callbackRef = intTri->edgeSplitCallbackList.insert(std::end(intTri->edgeSplitCallbackList), updateOnSplit);
    ```



## Citation

For the signpost data structure and most of the algorithms presented above, please cite

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
