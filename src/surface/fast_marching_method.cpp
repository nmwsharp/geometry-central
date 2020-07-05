#include "geometrycentral/surface/fast_marching_method.h"

#include <queue>
#include <tuple>


namespace geometrycentral {
namespace surface {

namespace {

// The super fun quadratic distance function in the Fast Marching Method on triangle meshes
// TODO parameter c isn't actually defined in paper, so I guessed that it was an error
double eikonalDistanceSubroutine(double a, double b, double theta, double dA, double dB) {


  if (theta <= PI / 2.0) {
    double u = dB - dA;
    double cTheta = std::cos(theta);
    double sTheta2 = 1.0 - cTheta * cTheta;

    // Quadratic equation
    double quadA = a * a + b * b - 2 * a * b * cTheta;
    double quadB = 2 * b * u * (a * cTheta - b);
    double quadC = b * b * (u * u - a * a * sTheta2);
    double sqrtVal = std::sqrt(quadB * quadB - 4 * quadA * quadC);
    // double tVals[] = {(-quadB + sqrtVal) / (2*quadA),        // seems to always be the first one
    //                   (-quadB - sqrtVal) / (2*quadA)};

    double t = (-quadB + sqrtVal) / (2 * quadA);
    if (u < t && a * cTheta < b * (t - u) / t && b * (t - u) / t < a / cTheta) {
      return dA + t;
    } else {
      return std::min(b + dA, a + dB);
    }

  }
  // Custom by Nick to get acceptable results in obtuse triangles without fancy unfolding
  else {

    double maxDist = std::max(dA, dB); // all points on base are less than this far away, by convexity
    double c = std::sqrt(a * a + b * b - 2 * a * b * std::cos(theta));
    double area = 0.5 * std::sin(theta) * a * b;
    double altitude = 2 * area / c; // distance to base, must be inside triangle since obtuse
    double baseDist = maxDist + altitude;

    return std::min({b + dA, a + dB, baseDist});
  }
}


} // namespace


VertexData<double> FMMDistance(IntrinsicGeometryInterface& geometry,
                               const std::vector<std::pair<Vertex, double>>& initialDistances) {

  typedef std::pair<double, Vertex> Entry;

  SurfaceMesh& mesh = geometry.mesh;
  geometry.requireEdgeLengths();
  geometry.requireCornerAngles();
  
  // TODO this could handle nonmanifold geometry with a few small tweaks
  if(!mesh.isManifold()) {
    throw std::runtime_error("handling of nonmanifold mesh not yet implemented");
  }

  // Initialize
  VertexData<double> distances(mesh, std::numeric_limits<double>::infinity());
  VertexData<char> finalized(mesh, false);

  std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> frontierPQ;
  for (auto& x : initialDistances) {
    frontierPQ.push(std::make_pair(x.second, x.first));
  }
  size_t nFound = 0;
  size_t nVert = mesh.nVertices();

  // Search
  while (nFound < nVert && !frontierPQ.empty()) {

    // Pop the nearest element
    Entry currPair = frontierPQ.top();
    frontierPQ.pop();
    Vertex currV = currPair.second;
    double currDist = currPair.first;


    // Accept it if not stale
    if (finalized[currV]) {
      continue;
    }
    distances[currV] = currDist;
    finalized[currV] = true;
    nFound++;


    // Add any eligible neighbors
    for (Halfedge he : currV.incomingHalfedges()) {
      Vertex neighVert = he.vertex();

      // Add with length
      if (!finalized[neighVert]) {
        double newDist = currDist + geometry.edgeLengths[he.edge()];
        if (newDist < distances[neighVert]) {
          frontierPQ.push(std::make_pair(currDist + geometry.edgeLengths[he.edge()], neighVert));
          distances[neighVert] = newDist;
        }
        continue;
      }

      // Check the third point of the "left" triangle straddling this edge
      if (he.isInterior()) {
        Vertex newVert = he.next().next().vertex();
        if (!finalized[newVert]) {

          // Compute the distance
          double lenB = geometry.edgeLengths[he.next().next().edge()];
          double distB = currDist;
          double lenA = geometry.edgeLengths[he.next().edge()];
          double distA = distances[neighVert];
          double theta = geometry.cornerAngles[he.next().next().corner()];
          double newDist = eikonalDistanceSubroutine(lenA, lenB, theta, distA, distB);

          if (newDist < distances[newVert]) {
            frontierPQ.push(std::make_pair(newDist, newVert));
            distances[newVert] = newDist;
          }
        }
      }

      // Check the third point of the "right" triangle straddling this edge
      Halfedge heT = he.twin();
      if (heT.isInterior()) {
        Vertex newVert = heT.next().next().vertex();
        if (!finalized[newVert]) {

          // Compute the distance
          double lenB = geometry.edgeLengths[heT.next().edge()];
          double distB = currDist;
          double lenA = geometry.edgeLengths[heT.next().next().edge()];
          double distA = distances[neighVert];
          double theta = geometry.cornerAngles[heT.next().next().corner()];
          double newDist = eikonalDistanceSubroutine(lenA, lenB, theta, distA, distB);

          if (newDist < distances[newVert]) {
            frontierPQ.push(std::make_pair(newDist, newVert));
            distances[newVert] = newDist;
          }
        }
      }
    }
  }

  return distances;
}


} // namespace surface
} // namespace geometrycentral
