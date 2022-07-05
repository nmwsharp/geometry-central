#pragma once

#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/surface/trace_geodesic.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

namespace geometrycentral {
namespace surface {

class PoissonDiskSampler {

public:
  // ===== Constructor
  // TODO: Currently this takes in a ManifoldSurfaceMesh, because the random sampling uses the traceGeodesic() function,
  // which relies on there being a well-defined tangent space at every vertex. However, perhaps for each non-manifold
  // vertex, we could choose a random outgoing halfedge.
  PoissonDiskSampler(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geometry, double rCoef = 1.0,
                     int kCandidates = 30, std::vector<SurfacePoint> pointsToAvoid = std::vector<SurfacePoint>());

  // ===== The main function: Sample the surface using Poisson Disk Sampling.
  std::vector<SurfacePoint> sample();

  // ===== Public members
  double rCoef;    // minimum distance between samples, expressed as a multiple of the mean edge length
  int kCandidates; // number of candidate points chosen from the (r,2r)-annulus around each sample
  std::vector<SurfacePoint> pointsToAvoid;

private:
  // ===== Members
  ManifoldSurfaceMesh& mesh;
  VertexPositionGeometry& geometry;

  double r; // the actual minimum distance between samples

  std::vector<Face> faceFromEachComponent; // holds once face for each connected component in the mesh
  std::vector<SurfacePoint> activeList;    // holds candidate points
  std::vector<SurfacePoint> samples;       // holds finalized samples

  // Spatial lookup structure adapted from Nick Sharp
  typedef std::array<long long int, 3> SpatialKey;
  std::map<SpatialKey, std::vector<Vector3>> spatialBuckets;
  Vector3 mapCenter;

  // ===== Utility and auxiliary functions
  std::tuple<Vector3, Vector3> boundingBox();
  void storeFaceFromEachComponent();
  SurfacePoint generateCandidate(const SurfacePoint& xi) const;
  void addNewSample(const SurfacePoint& sample);
  SpatialKey positionKey(const Vector3& position) const;
  void addPointToSpatialLookup(const Vector3& newPoint);
  double nearestWithinRadius(const SurfacePoint& candidate) const;
  void sampleOnConnectedComponent(const Face& f);
};

} // namespace surface
} // namespace geometrycentral