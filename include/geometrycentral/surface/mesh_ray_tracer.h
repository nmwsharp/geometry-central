#pragma once

#include "geometrycentral/surface/geometry.h"


#include "nanort/nanort.h"

namespace geometrycentral {
namespace surface {

struct RayHitResult {
  bool hit;
  double tHit;
  Face face;
  Vector3 baryCoords;
};

class MeshRayTracer {
public:
  // Creates a new tracer and builds the acceleration structure
  MeshRayTracer(Geometry<Euclidean>* geometry);

  // Build the BVH for the current geometry. Called automatically after
  // construction, re-call if mesh changes.
  void buildBVH();

  // Trace a ray. Note: geometry should be identical to when BVH was constructed
  RayHitResult trace(Vector3 start, Vector3 dir);

private:
  HalfedgeMesh* mesh;
  Geometry<Euclidean>* geometry;

  // Data for the BVH
  std::vector<double> rawPositions;
  std::vector<unsigned int> rawFaces;
  nanort::BVHAccel<double, nanort::TriangleMesh<double>, nanort::TriangleSAHPred<double>,
                   nanort::TriangleIntersector<double>>
      accel;

  double tFar;
};

} // namespace surface
}; // namespace geometrycentral
