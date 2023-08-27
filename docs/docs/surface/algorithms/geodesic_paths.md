The function `traceGeodesic` allows one to compute straightest paths along a surface (i.e. _geodesic_ paths). 

Note that straightest paths depend only on the _intrinsic_ geometry of a surface (via the `IntrinsicGeometryInterface`). Therefore, these routines can be run on abstract geometric domains as well as traditional surfaces in 3D. However, these routines do assume that the domain is a `ManifoldSurfaceMesh`.

`#include "geometrycentral/surface/trace_geodesic.h"`


??? func "`#!cpp TraceGeodesicResult traceGeodesic(IntrinsicGeometryInterface& geom, SurfacePoint startP, Vector2 traceVec, const TraceOptions& traceOptions = defaultTraceOptions);`"

    Trace a geodesic path along a surface mesh.
    
    - `inputGeom`: the input geometry (as always, a `VertexPositionGeometry` is valid input)
    - `startP`: the point on the surface where the path should start
    - `traceVec`: the direction the path should proceed in, and the distance that it should travel
    - `traceOptions`: options to specify the behavior of `traceGeodesic` in various situations

The function `traceGeodesic` traces out a geodesic path starting at `startP` which proceeds in the direction of `traceVec` and has length equal to `traceVec.norm()`, unless the path intersects a boundary edge in which case it stops there. 
This is also known as the _exponential map_. (As an aside, `geometry-central` also provides procedures for computing the inverse of the exponential map, known as the [logarithmic map](/surface/algorithms/vector_heat_method/#logarithmic-map).)

### Example

```cpp
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/surface/trace_geodesic.h"

// Load a mesh
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;
std::tie(mesh, geometry) = readManifoldSurfaceMesh(filename);

Vertex v = mesh->vertex(0);
Vector2 traceVec = 3 * Vector2::fromAngle(M_PI/6);
SurfacePoint pathEndpoint = traceGeodesic(*geometry, SurfacePoint(v), traceVec).endPoint;
```

## Helper Types
### Options
Options are passed in to `traceGeodesic` via a `TraceOptions` object.

| Field | Default value |Meaning|
|---|---|---|
| `#!cpp bool includePath`| `false` | whether to return the entire path trajectory (as opposed to merely returning the path's endpoint) |
| `#!cpp bool errorOnProblem`| `false` | whether to throw exceptions if the procedure encounters degenerate geometry |
| `#!cpp EdgeData<bool>* barrierEdges`| `nullptr` | if set, paths will stop when they hit barrier edges |
| `#!cpp size_t maxIters`| `INVALID_IND` | if set, paths will stop after traversing through `maxIters` faces |


### Result
The result is returned as a `TraceGeodesicResult`, which has 5 fields:

| Field | Meaning|
|---|---|
| `#!cpp SurfacePoint endPoint`| the point the path ended at |
| `#!cpp std::vector<SurfacePoint> pathPoints`| all points along the path, including start and end |
| `#!cpp Vector2 endingDir`| the incoming direction to the final point, in its tangent space |
| `#!cpp bool hitBoundary`| did the path stop early because we hit a boundary? |
| `#!cpp bool hasPath`| is `pathPoints` populated? |
| `#!cpp double length`| length of the traced path (generally equals norm of `traceVec` unless tracing stopped early due to hitting a boundary/barrier edge or due to the iteration limit `maxIters`) |


## Tangent Spaces
The input `traceVec` is specified as a vector in the _tangent space_ of the starting point. The meaning of this vector depends on whether the starting point is located on a vertex, edge, or face of the mesh. Tangent space in geometry central are discussed in more detail [on the Quantities page](/surface/geometry/quantities/#tangent-vectors-and-transport), but we give a brief overview here.

Given any mesh element (i.e. vertex, edge, or face) `p`, the $x$-axis of the tangent space at `p` points in the direction of `p.halfedge()`. The $y$-axis then points 90 degrees counterclockwise from the $x$-axis. (This is slightly more complicated at vertices, where one must use rescaled corner angles to define these directions. See the discussion of [vertex tangent spaces](/surface/geometry/quantities/#vertex-tangent-spaces)) for more details.
