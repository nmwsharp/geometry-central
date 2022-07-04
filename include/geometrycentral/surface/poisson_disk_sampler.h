#pragma once

#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/surface/trace_geodesic.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

namespace geometrycentral {
namespace surface {

class PoissonDiskSampler {

public:
  // ===== Constructor
  PoissonDiskSampler(SurfaceMesh& mesh, VertexPositionGeometry& geometry, double rCoef = 1.0, int kCandidates = 30,
                     std::vector<SurfacePoint> pointsToAvoid = std::vector<SurfacePoint>());

  // ===== The main function: Sample the surface using Poisson Disk Sampling.
  std::vector<SurfacePoint> sample();

  // ===== Public members
  double rCoef;    // minimum distance between samples, expressed as a multiple of the mean edge length
  int kCandidates; // number of candidate points chosen from the (r,2r)-annulus around each sample
  std::vector<SurfacePoint> pointsToAvoid;

private:
  // ===== Members
  SurfaceMesh& mesh;
  VertexPositionGeometry& geometry;

  double r; // the actual minimum distance between samples

  std::vector<SurfacePoint> activeList;
  std::vector<SurfacePoint> samples;

  // Spatial lookup structure adapted from Nick Sharp
  typedef std::array<long long int, 3> SpatialKey;
  std::map<SpatialKey, std::vector<Vector3>> spatialBuckets; // each bucket can hold at most 1 point
  Vector3 mapCenter;

  // ===== Utility and auxiliary functions
  SurfacePoint generateCandidate(const SurfacePoint& xi) const;
  void addNewSample(const SurfacePoint& sample);
  SpatialKey positionKey(const Vector3& position) const;
  void addPointToSpatialLookup(const Vector3& newPoint);
  double nearestWithinRadius(const SurfacePoint& candidate) const;
  std::vector<SurfacePoint> sampleOnConnectedComponent();
};

} // namespace surface
} // namespace geometrycentral