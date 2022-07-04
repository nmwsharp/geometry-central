#include "geometrycentral/surface/poisson_disk_sampler.h"

namespace geometrycentral {
namespace surface {

PoissonDiskSampler::PoissonDiskSampler(SurfaceMesh& mesh_, VertexPositionGeometry& geometry_, double rCoef_,
                                       int kCandidates_, std::vector<SurfacePoint> pointsToAvoid_)
    : mesh(mesh_), geometry(geometry_), rCoef(rCoef_), kCandidates(kCandidates_), pointsToAvoid(pointsToAvoid_) {

  // Compute mean edge length and set <r>.
  geometry.requireEdgeLengths();
  double meanEdgeLength = 0.;
  for (Edge e : mesh.edges()) {
    meanEdgeLength += geometry.edgeLengths[e];
  }
  meanEdgeLength /= mesh.nEdges();
  r = rCoef * meanEdgeLength;
  geometry.unrequireEdgeLengths();

  // Prepare spatial lookup
  Vector3 bboxMin, bboxMax;
  geometry.boundingBox(bboxMin, bboxMax);
  mapCenter = 0.5 * (bboxMin + bboxMax);
}

/*
 * Given a sample <xi>, generate a candidate point by tracing a geodesic with a randomly-chosen distance and direction.
 */
SurfacePoint PoissonDiskSampler::generateCandidate(const SurfacePoint& xi) const {

  SurfacePoint pathEndpoint;
  TraceGeodesicResult trace;

  Vector2 dir = Vector2::fromAngle(randomReal(0., 2. * PI));
  double dist = std::sqrt(randomReal(r, 2. * r));
  trace = traceGeodesic(geometry, xi, dist * dir);

  pathEndpoint = trace.endPoint;
  return pathEndpoint;
}

void PoissonDiskSampler::addNewSample(const SurfacePoint& sample) {

  samples.push_back(sample);
  activeList.push_back(sample);
  addPointToSpatialLookup(sample.interpolate(geometry.vertexPositions));
}

// ===== Helpers for spatial lookup. This code was pretty much directly taken from Nick Sharp.

PoissonDiskSampler::SpatialKey PoissonDiskSampler::positionKey(const Vector3& position) const {

  // It doesn't really matter that each bucket have sidelength r/sqrt(3), so that it can hold at most 1 sample. We would
  // have to translate between SurfacePoints and their positions more anyway.
  Vector3 coord = (position - mapCenter) / r;
  long long int keyX = static_cast<long long int>(std::floor(coord.x));
  long long int keyY = static_cast<long long int>(std::floor(coord.y));
  long long int keyZ = static_cast<long long int>(std::floor(coord.z));

  return {keyX, keyY, keyZ};
}

void PoissonDiskSampler::addPointToSpatialLookup(const Vector3& newPoint) {

  SpatialKey key = positionKey(newPoint);

  // Make sure the list exists (does nothing if already present)
  std::vector<Vector3> newV;
  spatialBuckets.insert(std::make_pair(key, newV));

  // Insert
  std::vector<Vector3>& bucketList = spatialBuckets[key];
  bucketList.push_back(newPoint);
}

/*
 * Given a candidate point, return the smallest distance to any of the existing samples.
 */
double PoissonDiskSampler::nearestWithinRadius(const SurfacePoint& candidate) const {

  Vector3 queryPoint = candidate.interpolate(geometry.vertexPositions);
  double minDistance = std::numeric_limits<double>::infinity();
  SpatialKey queryKey = positionKey(queryPoint);

  // Check all adjacent buckets
  for (int offsetX = -1; offsetX <= 1; offsetX++) {
    for (int offsetY = -1; offsetY <= 1; offsetY++) {
      for (int offsetZ = -1; offsetZ <= 1; offsetZ++) {

        SpatialKey bucketKey = queryKey;
        bucketKey[0] += offsetX;
        bucketKey[1] += offsetY;
        bucketKey[2] += offsetZ;

        auto bucketIter = spatialBuckets.find(bucketKey);
        if (bucketIter != spatialBuckets.end()) {
          for (Vector3 testP : bucketIter->second) {
            double dist = (queryPoint - testP).norm();
            minDistance = std::min(minDistance, dist);
          }
        }
      }
    }
  }
  return minDistance;
}

/*
 * The workhorse function.
 */
std::vector<SurfacePoint> PoissonDiskSampler::sampleOnConnectedComponent() {

  // TODO: Add points to avoid.
  // TODO: Get initial sample.
  activeList.clear();
  samples.clear();

  while (activeList.size() > 0) {
    // Select a random point from the active list.
    size_t randIdx = randomIndex(activeList.size());
    SurfacePoint& xi = activeList[randIdx];

    // Generate k samples between a distance of r and 2r away from x_i.
    // For each sample, check if it is within r of any points in <samples>.
    // If a sample is further than r away from all existing samples,
    // add it to <samples> and <activeList>.
    bool gotNewSample = false;
    for (size_t i = 0; i < kCandidates; i++) {

      SurfacePoint candidate = generateCandidate(xi);

      double dist = nearestWithinRadius(candidate);
      if (dist > r) {
        addNewSample(candidate);
        gotNewSample = true;
        break;
      }
    }
    if (gotNewSample) continue;

    // If none of the k samples is at least r away from all existing samples, remove xi from active_list.
    // [https://stackoverflow.com/a/4442529]
    if (randIdx != activeList.size() - 1) {
      activeList[randIdx] = std::move(activeList.back());
    }
    activeList.pop_back();
  }
}

/*
 * The final function.
 */
std::vector<SurfacePoint> PoissonDiskSampler::sample() {
  // Carry out sampling process for each connected component.
}

} // namespace surface
} // namespace geometrycentral