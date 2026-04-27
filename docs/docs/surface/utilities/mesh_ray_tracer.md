A BVH-accelerated ray-triangle intersection utility for triangle meshes, backed by [nanort](https://github.com/lighttransport/nanort).

`#!cpp #include "geometrycentral/surface/mesh_ray_tracer.h"`

```cpp
#include "geometrycentral/surface/mesh_ray_tracer.h"
using namespace geometrycentral::surface;

// Build the acceleration structure once
MeshRayTracer tracer(geometry);

// Trace a ray — origin and direction, direction need not be normalized
RayHitResult hit = tracer.trace({0.5, 0.5, 10.0}, {0.0, 0.0, -1.0});

if (hit.hit) {
  // distance along ray
  double t = hit.t;

  // which triangle was hit
  size_t faceIdx = hit.faceIndex;

  // barycentric coords in that triangle
  Vector3 bary = hit.faceCoords;
}
```

The mesh must be triangulated. The BVH is built once at construction and can be queried any number of times.

### Result struct

`RayHitResult` is returned by every `trace()` call.

```cpp
struct RayHitResult {
  bool hit;           // true if the ray intersected the mesh
  double t;           // distance along the ray to the hit point
  size_t faceIndex;   // index of the hit face in the mesh
  Vector3 faceCoords; // barycentric coords in the hit face, in the same order as the face's vertices are given
};
```

### Construction

??? func "`#!cpp MeshRayTracer::MeshRayTracer(EmbeddedGeometryInterface& geom)`"

    Build the BVH over the triangles of `geom`'s mesh. The mesh must be fully triangular.

??? func "`#!cpp MeshRayTracer::MeshRayTracer(const SimplePolygonMesh& mesh)`"

    Build the BVH from a `SimplePolygonMesh`. All polygons must be triangles.

### Queries

??? func "`#!cpp RayHitResult MeshRayTracer::trace(Vector3 origin, Vector3 dir) const`"

    Trace a ray from `origin` in direction `dir` (need not be normalized).

### Utilities

??? func "`#!cpp SurfacePoint toSurfacePoint(const RayHitResult& hit, SurfaceMesh& mesh)`"

    Convert a `RayHitResult` to a `SurfacePoint` on `mesh`. 
    
    Only meaningful when the tracer was built from that mesh (or its geometry). Can only be called if the result was a hit.
    
    The returned point is a face point with the hit's barycentric coordinates.
