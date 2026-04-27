#pragma once

#include "geometrycentral/surface/embedded_geometry_interface.h"
#include "geometrycentral/surface/simple_polygon_mesh.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/utilities/vector3.h"

#include <cstddef>
#include <limits>
#include <memory>

namespace geometrycentral {
namespace surface {

struct RayHitResult {
  bool hit = false;
  double t = std::numeric_limits<double>::infinity();    // distance along the ray to the hit
  size_t faceIndex = std::numeric_limits<size_t>::max(); // index of the hit face
  Vector3 faceCoords{0., 0., 0.};                        // barycentric coords in the hit face, in the
                                                         // same order as the face's vertices are given
};

// Acceleration structure for ray-triangle intersection queries against a triangle mesh.
// Backed by a nanort BVH. The mesh must be triangulated.
class MeshRayTracer {
public:
  MeshRayTracer(EmbeddedGeometryInterface& geom);
  MeshRayTracer(const SimplePolygonMesh& mesh);
  ~MeshRayTracer();

  // Trace a ray from `origin` in direction `dir` (need not be normalized).
  // Returns hit info; result.hit is false if no intersection is found.
  RayHitResult trace(Vector3 origin, Vector3 dir) const;

private:
  class Impl;
  std::unique_ptr<Impl> impl;
};

// Convert a RayHitResult to a SurfacePoint on the given mesh.
// Only valid when the tracer was built from that mesh (or its geometry).
SurfacePoint toSurfacePoint(const RayHitResult& hit, SurfaceMesh& mesh);

} // namespace surface
} // namespace geometrycentral
