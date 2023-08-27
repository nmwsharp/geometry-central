#pragma once

#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/surface/trace_geodesic.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

namespace geometrycentral {
namespace surface {

class PoissonDiskSampler {

public:
  // ===== Constructor

  // Currently this takes in a ManifoldSurfaceMesh, because the random sampling uses the traceGeodesic() function,
  // which relies on there being a well-defined tangent space at every vertex. However, perhaps for each non-manifold
  // vertex, we could choose a random outgoing halfedge.
  PoissonDiskSampler(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geometry);

  // ===== The main function: Sample the surface using Poisson Disk Sampling.
  std::vector<SurfacePoint> sample(double rCoef = 1.0, int kCandidates = 30,
                                   std::vector<SurfacePoint> pointsToAvoid = std::vector<SurfacePoint>(),
                                   int rAvoidance = 1, bool use3DAvoidanceRadius = true);

private:
  // ===== Members
  ManifoldSurfaceMesh& mesh;
  VertexPositionGeometry& geometry;

  double rCoef;    // minimum distance between samples, expressed as a multiple of the mean edge length
  int kCandidates; // number of candidate points chosen from the (r,2r)-annulus around each sample
  std::vector<SurfacePoint> pointsToAvoid;

  double rMinDist;       // the actual minimum distance between samples
  double meanEdgeLength; // the mean edge length
  double sideLength;     // side length of each bucket

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
  void addPointToSpatialLookupWithRadius(const SurfacePoint& newPoint, int R = 0, bool use3DAvoidanceRadius = true);
  bool isCandidateValid(const SurfacePoint& candidate) const;
  void sampleOnConnectedComponent(const Face& f);
  void clearData();
};

} // namespace surface
} // namespace geometrycentral