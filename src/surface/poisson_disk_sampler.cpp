#include "geometrycentral/surface/poisson_disk_sampler.h"

namespace geometrycentral {
namespace surface {

PoissonDiskSampler::PoissonDiskSampler(ManifoldSurfaceMesh& mesh_, VertexPositionGeometry& geometry_)
    : mesh(mesh_), geometry(geometry_) {

  // Compute mean edge length.
  geometry.requireEdgeLengths();
  meanEdgeLength = 0.;
  for (Edge e : mesh.edges()) {
    meanEdgeLength += geometry.edgeLengths[e];
  }
  meanEdgeLength /= mesh.nEdges();
  geometry.unrequireEdgeLengths();

  // Prepare spatial lookup
  Vector3 bboxMin, bboxMax;
  std::tie(bboxMin, bboxMax) = boundingBox();
  mapCenter = 0.5 * (bboxMin + bboxMax);

  // Determine connected components.
  storeFaceFromEachComponent();
}

/*
 * Return one of the two pairs of furthest-away corners in the axis-aligned bounding box of the mesh.
 *
 * I.e., the distance between bboxMin and bboxMax is the length of the bbox's diagonal, and the length of the bbox in
 * the i-th dimension is |bboxMin[i] - bboxMax[i]|.
 */
std::tuple<Vector3, Vector3> PoissonDiskSampler::boundingBox() {

  const double inf = std::numeric_limits<double>::infinity();
  Vector3 bboxMin = {inf, inf, inf};
  Vector3 bboxMax = {-inf, -inf, -inf};
  for (Vertex v : mesh.vertices()) {
    Vector3& pos = geometry.vertexPositions[v];
    if (pos[0] <= bboxMin[0] && pos[1] <= bboxMin[1] && pos[2] <= bboxMin[2]) {
      bboxMin = pos;
    }
    if (pos[0] >= bboxMax[0] && pos[1] >= bboxMax[1] && pos[2] >= bboxMax[2]) {
      bboxMax = pos;
    }
  }
  return std::make_tuple(bboxMin, bboxMax);
}

/*
 * Get one face from each connected component in the mesh.
 */
void PoissonDiskSampler::storeFaceFromEachComponent() {

  // Do simple depth-first search.
  faceFromEachComponent.clear();
  FaceData<char> marked(mesh, false); // has face been visited?

  for (Face f : mesh.faces()) {
    // Get the next unmarked face.
    if (marked[f]) continue;

    // We've never visited this face before; it must belong to a new component.
    marked[f] = true;
    faceFromEachComponent.push_back(f);

    // Begin gathering each face in the current component.
    std::vector<Face> queue = {f};
    Face currF;
    while (!queue.empty()) {

      currF = queue.back();
      queue.pop_back();

      for (Face g : currF.adjacentFaces()) {
        if (marked[g]) continue;
        marked[g] = true;
        queue.push_back(g);
      }
    }
  }
}

/*
 * Given a sample <xi>, generate a candidate point by tracing a geodesic with a randomly-chosen distance and direction.
 */
SurfacePoint PoissonDiskSampler::generateCandidate(const SurfacePoint& xi) const {

  SurfacePoint pathEndpoint;
  TraceGeodesicResult trace;

  Vector2 dir = Vector2::fromAngle(randomReal(0., 2. * PI));
  double dist = std::sqrt(randomReal(rMinDist, 2. * rMinDist));
  trace = traceGeodesic(geometry, xi, dist * dir);

  pathEndpoint = trace.endPoint;
  return pathEndpoint;
}

void PoissonDiskSampler::addNewSample(const SurfacePoint& sample) {

  samples.push_back(sample);
  activeList.push_back(sample);
  addPointToSpatialLookup(sample.interpolate(geometry.vertexPositions));
}

// ===== Helpers for spatial lookup. This code was pretty much directly taken from Nick Sharp's Poisson sampling code.

PoissonDiskSampler::SpatialKey PoissonDiskSampler::positionKey(const Vector3& position) const {

  // Each bucket has sidelength s = r/sqrt(3) so that it can hold at most 1 sample.
  Vector3 coord = (position - mapCenter) / sideLength;
  long long int keyX = static_cast<long long int>(std::floor(coord.x));
  long long int keyY = static_cast<long long int>(std::floor(coord.y));
  long long int keyZ = static_cast<long long int>(std::floor(coord.z));

  return {keyX, keyY, keyZ};
}

void PoissonDiskSampler::addPointToSpatialLookup(const Vector3& newPos) {

  SpatialKey key = positionKey(newPos);
  spatialBuckets.insert(std::make_pair(key, newPos));
}

/*
 * Mark all buckets with a radius of <R> buckets as occupied, as well.
 */
void PoissonDiskSampler::addPointToSpatialLookupWithRadius(const SurfacePoint& newPoint, int R,
                                                           bool use3DAvoidanceRadius) {

  Vector3 newPos = newPoint.interpolate(geometry.vertexPositions);
  addPointToSpatialLookup(newPos);

  if (R <= 0) return;

  if (!use3DAvoidanceRadius) {
    // This places fictitious points in a metric ball approximately of radius R*r centered at <newPoint>.
    // The solid ball is built by constructing R layers, radially outward; each layer is spaced r apart, and
    // points in each layer are spaced approx. r apart around the circular layer. The idea is that no point can be added
    // within this ball of points without being at most approx. r * sqrt(3)/2 < r away from another point with the ball.
    // In practice, the curvature of the surface may be s.t. geodesics traced in this manner don't end up being spaced
    // evenly the way we want; heuristically, just divide the angle increment by 2, which should work for all but very
    // "spiky" meshes.
    SurfacePoint pathEndpoint;
    TraceGeodesicResult trace;
    for (int rIter = 0; rIter <= R; rIter++) {
      double dist = rIter * rMinDist;
      for (double theta = 0.; theta <= 2. * PI; theta += rMinDist / dist / 2.) {
        trace = traceGeodesic(geometry, newPoint, dist * Vector2::fromAngle(theta));
        pathEndpoint = trace.endPoint;
        addPointToSpatialLookup(pathEndpoint.interpolate(geometry.vertexPositions));
      }
    }
  } else {
    // This places fictitious points in a *solid 3D ball* approximately of radius R*r centered at <newPoint>.
    // The solid ball is built by constructing R layers, radially outward; each layer is spaced r apart, and
    // points in each layer are spaced approx. r apart in a "grid" (a mapping from the cylinder to the sphere that has
    // been scaled s.t. projected points end up being approx. r apart.)
    for (int rIter = 0; rIter <= R; rIter++) {
      double r = rIter * rMinDist;
      for (double z = -0.99; z <= 0.99; z += rMinDist) {
        double coeff = std::sqrt(1. - z * z);
        for (double theta = 0.0; theta <= 2.0 * PI; theta += rMinDist / coeff) {
          Vector3 pos = {r * coeff * std::cos(theta), r * coeff * std::sin(theta), r * z};
          pos += newPos;
          addPointToSpatialLookup(pos);
        }
      }
    }
  }
}

/*
 * Given a candidate point, return "false" if the smallest distance to any of the existing samples is less than r.
 * Return "true" otherwise.
 */
bool PoissonDiskSampler::isCandidateValid(const SurfacePoint& candidate) const {

  Vector3 queryPoint = candidate.interpolate(geometry.vertexPositions);
  SpatialKey queryKey = positionKey(queryPoint);

  // Check all adjacent buckets. Each bucket's sidelength is chosen to be small enough such that if any two adjacent
  // buckets contain points, their points are guaranteed to be within r of each other.
  for (int offsetX = -1; offsetX <= 1; offsetX++) {
    for (int offsetY = -1; offsetY <= 1; offsetY++) {
      for (int offsetZ = -1; offsetZ <= 1; offsetZ++) {

        SpatialKey bucketKey = queryKey;
        bucketKey[0] += offsetX;
        bucketKey[1] += offsetY;
        bucketKey[2] += offsetZ;

        auto bucketIter = spatialBuckets.find(bucketKey);
        if (bucketIter != spatialBuckets.end()) {
          if ((queryPoint - bucketIter->second).norm() < rMinDist) {
            return false;
          }
        }
      }
    }
  }
  return true;
}

/*
 * The workhorse function. Takes as input a face from a connected component.
 */
void PoissonDiskSampler::sampleOnConnectedComponent(const Face& f) {

  // Get initial sample.
  Vector3 baryCoords = {1. / 3., 1. / 3., 1. / 3.};
  SurfacePoint x0(f, baryCoords);
  activeList.clear();
  addNewSample(x0);

  while (activeList.size() > 0) {
    // Select a random point from the active list.
    size_t randIdx = randomIndex(activeList.size());
    SurfacePoint& xi = activeList[randIdx];

    // Generate k samples between a distance of r and 2r away from x_i.
    // For each sample, check if it is within r of any points in <samples>.
    // If a sample is further than r away from all existing samples, add it to <samples> and <activeList>.
    bool gotNewSample = false;
    for (int i = 0; i < kCandidates; i++) {

      SurfacePoint candidate = generateCandidate(xi);

      if (isCandidateValid(candidate)) {
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
 * Reset all data structures in preparation for a new solve.
 */
void PoissonDiskSampler::clearData() {

  samples.clear();
  spatialBuckets.clear();
}

/*
 * The final function.
 */
std::vector<SurfacePoint> PoissonDiskSampler::sample(double rCoef_, int kCandidates_,
                                                     std::vector<SurfacePoint> pointsToAvoid_, int rAvoidance,
                                                     bool use3DAvoidanceRadius) {

  clearData();

  // Set parameters.
  rCoef = rCoef_;
  rMinDist = rCoef * meanEdgeLength;
  sideLength = rMinDist / std::sqrt(3.0);
  kCandidates = kCandidates_;
  pointsToAvoid = pointsToAvoid_;

  // Add points to avoid.
  for (const SurfacePoint& pt : pointsToAvoid) {
    addPointToSpatialLookupWithRadius(pt, rAvoidance - 1, use3DAvoidanceRadius);
  }

  // Carry out sampling process for each connected component.
  for (const Face& f : faceFromEachComponent) {
    sampleOnConnectedComponent(f);
  }

  return samples;
}

} // namespace surface
} // namespace geometrycentral