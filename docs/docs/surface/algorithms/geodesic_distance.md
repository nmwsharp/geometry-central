This section describes algorithms for computing distance along a surface, or _geodesic_ distance. 

Note that distance depends on the _intrinsic_ geometry of a surface (via the `IntrinsicGeometryInterface`). Therefore, these routines can be run on abstract geometric domains as well as traditional surfaces in 3D.

## Polyhedral Distance

TODO document

## Heat Method for Distance

These routines implement the [Heat Method for Geodesic Distance](http://www.cs.cmu.edu/~kmcrane/Projects/HeatMethod/paper.pdf). This algorithm uses short time heat flow to compute distance on surfaces. Because the main burden is simply solving linear systems of equations, it tends to be faster than polyhedral schemes, especially when computing distance multiple times on the same surface.  In the computational geometry sense, this method is an approximation, as the result is not precisely equal to the polyhedral distance on the surface; nonetheless it is fast and well-suited for many applications.

This class supports any (possibly-nonmanifold) triangle mesh as input, and requires only intrinsic geometry (aka edge lengths) to function. Furthermore, it can optionally use robust operators internally to improve performance on low-quality meshes.

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
