This section describes algorithms for computing distance along a surface, or _geodesic_ distance.

Note that distance depends on the _intrinsic_ geometry of a surface (via the `IntrinsicGeometryInterface`). Therefore, these routines can be run on abstract geometric domains as well as traditional surfaces in 3D.

## Polyhedral Distance

These routines use [Danil Kirsanov's implementation](https://code.google.com/archive/p/geodesic/) of the [MMP algorithm for exact polyhedral geodesic distance](https://www.cs.jhu.edu/~misha/ReadingSeminar/Papers/Mitchell87.pdf) to compute geodesic distance (and geodesic paths) along a surface.

This class supports any (possibly-nonmanifold) triangle mesh as input, and requires only intrinsic geometry (aka edge lengths) to function.

`#include "geometrycentral/surface/exact_geodesics.h"`

### Single Solves

Geodesic distance from a single source vertex can be computed via the following utility function.
More general source data, queries in the interior of triangles, and geodesic path extraction can be handled using the stateful version below.

Example
```cpp
#include "geometrycentral/surface/exact_geodesics.h"
#include "geometrycentral/surface/meshio.h"

// Load a mesh
std::unique_ptr<SurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;
std::tie(mesh, geometry) = readSurfaceMesh(filename);

// Pick a vertex
Vertex sourceVert = /* some vertex */

// Compute distance
VertexData<double> distToSource = exactGeodesicDistance(*mesh, *geometry, sourceVert);
/* do something useful */
```

??? func "`#!cpp VertexData<double> exactGeodesicDistance(SurfaceMesh& mesh, IntrinsicGeometryInterface& geom, Vertex v)`"

    Compute the distance from the source using MMP. See the stateful class below for further options.

### Advanced Queries

The stateful class `GeodesicAlgorithmExact` runs the MMP algorithm to compute geodesic distance from a given set of source points. The resulting distance field can be queried at any point on the input mesh to find the identity of the nearest source point, the distance to the source point, and the shortest path to the source point.

Both the source points and query point are represented as `SurfacePoints` (see [here](../../utilities/surface_point/)), i.e. locations on a surface which may be a vertex, a point along an edge, or a point inside a face.

Note that unlike the [heat method](#heat-method-for-distance), this precomputation does not speed up future distance computations using different sets of source points.

Example:
```cpp
#include "geometrycentral/surface/exact_geodesics.h"
#include "geometrycentral/surface/meshio.h"

// Load a mesh
std::unique_ptr<SurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;
std::tie(mesh, geometry) = readSurfaceMesh(filename);

// Create the GeodesicAlgorithmExact
GeodesicAlgorithmExact mmp(*mesh, *geometry);

// Pick a few points as the source set
std::vector<SurfacePoint> sourcePoints;
Vertex v = /* some vertex */
sourcePoints.push_back(SurfacePoint(v));

Edge e = /* some edge */
double tEdge = /* some coordinate along edge e*/
sourcePoints.push_back(SurfacePoint(e, tEdge));

Face f =  /* some face */;
Vector3 fBary =  /* some barycentric coords in face f */;
sourcePoints.push_back(SurfacePoint(f, fBary));

// Run MMP from these source points
mmp.propagate(sourcePoints);

// Get the distance function at all mesh vertices
VertexData<double> distToSource = mmp.getDistanceFunction();

// Query the distance function at some point
SurfacePoint queryPoint = /* some point on the surface */
double dist = mmp.getDistance(queryPoint);

// Get the geodesic path from a query point to the nearest source
SurfacePoint queryPoint2 = /* some point on the surface */
double pathLength;
std::vector<SurfacePoint> path = mmp.traceBack(queryPoint2, pathLength);

// Get the gradient of the distance function at some point
Vector2 gradDist = mmp.getDistanceGradient(queryPoint);

// Find the tangent vector pointing to a query point from the closest source
Vector2 log_base_source_of_query = mmp.getLog(queryPoint);
```


??? func "`#!cpp GeodesicAlgorithmExact::GeodesicAlgorithmExact(SurfaceMesh& mesh, IntrinsicGeometryInterface& geom)`"

    Creates a new solver, but does not do any computation

??? func "`#!cpp void GeodesicAlgorithmExact::propagate(const std::vector<SurfacePoint>& sources, double max_propagation_distance = GEODESIC_INF, const std::vector<SurfacePoint>& stop_points = {})`"

    Compute the distance field from a set of source vertices. This distance field can then be queried by several functions below.

    - `sources` is a list of points on the surface. The computed distance field gives the distance from any point on the mesh to the closest of these source points.

    - `max_propagation_distance` is a cutoff value. Once `propagate` identifies all vertices within `max_propagation_distance` of the source points, it stops. By default, `max_propagation_distance` is set to infinity so that distances are computed across the whole input surface.

    - `stop_points` is a list of points on the mesh at which we want to know the distance function. By default it is empty. If it is nonempty, then propagate stops once accurate distances at all of the stop points have been computed. This can speed up run time considerably if one is only interested in a small set of points near the source.

??? func "`#!cpp void GeodesicAlgorithmExact::propagate(const std::vector<Vertex>& sources, double max_propagation_distance = GEODESIC_INF, const std::vector<Vertex>& stop_points = {})`"

    Performs the same computation as the first `propagate` function, but takes vertices rather than arbitrary surface points for convenience.

??? func "`#!cpp void GeodesicAlgorithmExact::propagate(const SurfacePoint& source, double max_propagation_distance = GEODESIC_INF, const std::vector<SurfacePoint>& stop_points = {})`"

    Performs the same computation as the first `propagate` function, but takes a single source point for convenience.

??? func "`#!cpp void GeodesicAlgorithmExact::propagate(const Vertex& source, double max_propagation_distance = GEODESIC_INF, const std::vector<Vertex>& stop_points = {})`"

    Performs the same computation as the first `propagate` function, but takes a single source vertex rather than arbitrary surface points for convenience.

??? func "`#!cpp std::vector<SurfacePoint> GeodesicAlgorithmExact::traceBack(const SurfacePoint& point, double& pathLength [optional]) const`"

    Compute the geodesic path from `point` to the closest source. This path is encoded as a list of `SurfacePoints` starting at `point` and ending at the source.
    (Optionally also returns the length of the path via the `pathLength` argument.)

??? func "`#!cpp std::vector<SurfacePoint> GeodesicAlgorithmExact::traceBack(const Vertex& v, double& pathLength [optional]) const`"

    Compute the geodesic path from `v` to the closest source point. This path is encoded as a list of `SurfacePoints` starting at `v` and ending at the source point.
    (Optionally also returns the length of the path via the `pathLength` argument.)

??? func "`#!cpp std::pair<unsigned, double> GeodesicAlgorithmExact::closestSource(const SurfacePoint& point) const`"

    Identify the closest source to `point` and compute the distance to that source. Returns the index of that source in the source list along with the distance.

??? func "`#!cpp std::pair<unsigned, double> GeodesicAlgorithmExact::closestSource(const Vertex& v) const`"

    Identify the closest source to `v` and compute the distance to that source. Returns the index of that source in the source list along with the distance.

??? func "`#!cpp double GeodesicAlgorithmExact::getDistance(const SurfacePoint& point) const`"

    Returns the distance from `point` to the closest source.

??? func "`#!cpp double GeodesicAlgorithmExact::getDistance(const Vertex& v) const`"

    Returns the distance from `v` to the closest source.

??? func "`#!cpp VertexData<double> GeodesicAlgorithmExact::getDistanceFunction() const`"

    Evaluate the distance function at every vertex of the mesh.

??? func "`#!cpp IntervalList GeodesicAlgorithmExact::getEdgeIntervals(Edge e) const`"

    Get the list of windows along edge `e` that MMP uses to represent the distance function.

### Vector-valued Queries

Computing exact geodesic paths also allows one to compute exact [log maps](/surface/algorithms/vector_heat_method/#logarithmic-map), as well as the exact gradient of the distance function.

??? func "`#!cpp Vector2 GeodesicAlgorithmExact::getLog(const SurfacePoint& point) const`"

    Returns the log map at the source closest to `point` (i.e. the tangent vector based at the closest source which points towards `point` whose magnitude is the distance from the source to `point`).

??? func "`#!cpp Vector2 GeodesicAlgorithmExact::getLog(const Vertex& v) const`"

    Returns the log map at the source closest to `v` (i.e. the tangent vector based at the closest source which points towards `v` whose magnitude is the distance from the source to `point`).

??? func "`#!cpp Vector2 GeodesicAlgorithmExact::getDistanceGradient(const SurfacePoint& point) const`"

    Returns the gradient of the distance function at `point` (i.e. the unit tangent vector at `point` which points away from the closest source)

??? func "`#!cpp Vector2 GeodesicAlgorithmExact::getDistanceGradient(const Vertex& v) const`"

    Returns the gradient of the distance function at `v` (i.e. the unit tangent vector at `v` which points away from the closest source)

## Heat Method for Distance

These routines implement the [Heat Method for Geodesic Distance](http://www.cs.cmu.edu/~kmcrane/Projects/HeatMethod/paper.pdf). This algorithm uses short time heat flow to compute distance on surfaces. Because the main burden is simply solving linear systems of equations, it tends to be faster than polyhedral schemes, especially when computing distance multiple times on the same surface.  In the computational geometry sense, this method is an approximation, as the result is not precisely equal to the polyhedral distance on the surface; nonetheless it is fast and well-suited for many applications.

This class also supports any (possibly-nonmanifold) triangle mesh as input, and requires only intrinsic geometry (aka edge lengths) to function. Furthermore, it can optionally use robust operators internally to improve performance on low-quality meshes.

`#include "geometrycentral/surface/heat_method_distance.h"`

### Single Solves

A one-off utility function is provided which computes the distance from a source vertex using the heat method. Repeated solves or more general source data should use the stateful version below.

Example
```cpp
#include "geometrycentral/surface/heat_method_distance.h"
#include "geometrycentral/surface/meshio.h"

// Load a mesh
std::unique_ptr<SurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;
std::tie(mesh, geometry) = loadMesh(filename);

// Pick a vertex
Vertex sourceVert = /* some vertex */

// Compute distance
VertexData<double> distToSource = heatMethodDistance(*geometry, sourceVert);
/* do something useful */
```

??? func "`#!cpp VertexData<double> heatMethodDistance(IntrinsicGeometryInterface& geom, Vertex v)`"

    Compute the distance from the source using the heat method. See the stateful class below for further options.


### Repeated Solves

The stateful class `HeatMethodDistanceSolver` does precomputation when constructed, then allows many distance solves from different source locations to be performed efficiently. This class also exposes options, like changing the internal short-time parameter, or using a robust operators.

The `computeDistance()` method in `HeatMethodDistanceSolver` can also take `SurfacePoint`(s) as the source location(s). A `SurfacePoint` (see [here](../../utilities/surface_point/)) is a location on a surface, which may be a vertex, a point along an edge, or a point inside a face.

Example:
```cpp
#include "geometrycentral/surface/heat_method_distance.h"
#include "geometrycentral/surface/meshio.h"

// Load a mesh
std::unique_ptr<SurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;
std::tie(mesh, geometry) = loadMesh(filename);

// Create the Heat Method solver
HeatMethodDistanceSolver heatSolver(*geometry);

// Alternately, set useRobustLaplacian=true to get a robustified version
// HeatMethodDistanceSolver heatSolver(*geometry, 1.0, true);

// Some vertices as source set
std::vector<Vertex> sourceVerts = /* some interesting vertices */
for(Vertex v : sourceVerts) {
  VertexData<double> distToSource = heatSolver.computeDistance(v);
  /* do something useful */
}


// A point in a face as a source set
Face sourceF =  /* some face */;
Vector3 sourceFBary =  /* some barycentric coords in face */;
SurfacePoint targetP(sourceF, sourceFBary);

VertexData<double> distToSource = heatSolver.computeDistance(targetP);
/* do something useful */
```


??? func "`#!cpp HeatMethodDistanceSolver::HeatMethodDistanceSolver(IntrinsicGeometryInterface& geom, double tCoef=1.0, bool useRobustLaplacian = false)`"

    Create a new solver to compute geodesic distance using the heat method. All precomputation work is performed immediately at construction time.

    - `geom` is the geometry (and hence mesh) on which to compute. Note that nearly any geometry object (`VertexPositionGeometry`, etc) can be passed here.

    - `tCoef` is the time to use for short time heat flow, as a factor `m * h^2`, where `h` is the mean edge length. The default value of `1.0` is almost always sufficient.

    - `useRobustLaplacian` is true, the solver will internally use a robust intrinsic Laplacian, including mollification & tufting for nonmanifold inputs. See "A Laplacian for Nonmanifold Triangle Meshes" [Sharp & Crane 2020 @ SGP] for algorithmic details and citation.

    Algorithm options (like `tCoef`) cannot be changed after construction; create a new solver object with the new settings.


??? func "`#!cpp VertexData<double> HeatMethodDistanceSolver::computeDistance(Vertex v)`"

    Compute the distance from a single source vertex.


??? func "`#!cpp VertexData<double> HeatMethodDistanceSolver::computeDistance(std::vector<Vertex> verts)`"

    Compute the distance from a set of source vertices.


??? func "`#!cpp VertexData<double> HeatMethodDistanceSolver::computeDistance(SurfacePoint p)`"

    Compute the distance from a single source point.


??? func "`#!cpp VertexData<double> HeatMethodDistanceSolver::computeDistance(std::vector<SurfacePoint> points)`"

    Compute the distance from a set of source points.
