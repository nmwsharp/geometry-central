#pragma once

#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/surface/trace_geodesic.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

namespace geometrycentral {
namespace surface {

struct PoissonDiskOptions {
  double minDist = 1.0; // minimum distance r between samples
  int kCandidates = 30; // number of candidate points chosen from the (r,2r)-annulus around each sample
  std::vector<SurfacePoint> pointsToAvoid;
  double minDistAvoidance = 1.0; // radius of avoidance
  bool use3DAvoidance = true;
};

class PoissonDiskSampler {

public:
  // ===== Constructor

  // Currently this takes in a ManifoldSurfaceMesh, because the random sampling uses the traceGeodesic() function,
  // which relies on there being a well-defined tangent space at every vertex. However, perhaps for each non-manifold
  // vertex, we could choose a random outgoing halfedge.
  PoissonDiskSampler(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geometry);

  // ===== The main function: Sample the surface using Poisson Disk Sampling.
  std::vector<SurfacePoint> sample(const PoissonDiskOptions& options = PoissonDiskOptions());

private:
  // ===== Members
  ManifoldSurfaceMesh& mesh;
  VertexPositionGeometry& geometry;

  int kCandidates; // number of candidate points chosen from the (r,2r)-annulus around each sample
  std::vector<SurfacePoint> pointsToAvoid;

  double rMinDist;   // the minimum distance between samples
  double sideLength; // side length of each bucket

  std::vector<Face> faceFromEachComponent; // holds one face for each connected component in the mesh
  std::vector<SurfacePoint> activeList;    // holds candidate points
  std::vector<SurfacePoint> samples;       // holds finalized samples

  // Spatial lookup structure
  typedef std::array<long long int, 3> SpatialKey;
  std::map<SpatialKey, Vector3> spatialBuckets;
  Vector3 mapCenter;

  // ===== Utility and auxiliary functions
  std::tuple<Vector3, Vector3> boundingBox();
  void storeFaceFromEachComponent();
  SurfacePoint generateCandidate(const SurfacePoint& xi) const;
  void addNewSample(const SurfacePoint& sample);
  SpatialKey positionKey(const Vector3& position) const;
  void addPointToSpatialLookup(const Vector3& newPos);
  void addPointToSpatialLookupWithRadius(const SurfacePoint& newPoint, double radius = 0, bool use3DAvoidance = true);
  bool isCandidateValid(const SurfacePoint& candidate) const;
  void sampleOnConnectedComponent(const Face& f);
  void clearData();
};

} // namespace surface
} // namespace geometrycentral