#include "geometrycentral/surface/fast_marching_method.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <tuple>


namespace geometrycentral {
namespace surface {

namespace {

// The super fun quadratic distance function in the Fast Marching Method on triangle meshes
// TODO parameter c isn't actually defined in paper, so I guessed that it was an error
double eikonalDistanceSubroutine(double a, double b, double theta, double dA, double dB, int sign) {


  if (theta <= PI / 2.0) {
    double u = dB - dA;
    double cTheta = std::cos(theta);
    double sTheta2 = 1.0 - cTheta * cTheta;

    // Quadratic equation
    double quadA = a * a + b * b - 2 * a * b * cTheta;
    double quadB = 2 * b * u * (a * cTheta - b);
    double quadC = b * b * (u * u - a * a * sTheta2);
    double sqrtVal = std::sqrt(quadB * quadB - 4 * quadA * quadC);
    double tVals[] = {(-quadB + sqrtVal) / (2 * quadA), (-quadB - sqrtVal) / (2 * quadA)};
    double t = sign > 0 ? tVals[0] : tVals[1];
    double y = b * (t - u) / t;
    if (u < sign * t && a * cTheta < y && y < a / cTheta) {
      return dA + t;
    } else {
      return std::min(b * sign + dA, a * sign + dB);
    }

  }
  // Custom by Nick to get acceptable results in obtuse triangles without fancy unfolding
  else {
    double maxDist = std::max(sign * dA, sign * dB); // all points on base are less than this far away, by convexity
    double c = std::sqrt(a * a + b * b - 2 * a * b * std::cos(theta));
    double area = 0.5 * std::sin(theta) * a * b;
    double altitude = 2 * area / c; // distance to base, must be inside triangle since obtuse
    double baseDist = maxDist + altitude;

    return std::min({b * sign + dA, a * sign + dB, sign * baseDist});
  }
}


} // namespace

VertexData<double> FMMDistance(IntrinsicGeometryInterface& geometry,
                               const std::vector<std::vector<std::pair<SurfacePoint, double>>>& initialDistances,
                               bool sign) {

  typedef std::tuple<double, Vertex, int> Entry;

  SurfaceMesh& mesh = geometry.mesh;
  geometry.requireEdgeLengths();
  geometry.requireCornerAngles();

  // TODO this could handle nonmanifold geometry with a few small tweaks
  if (!mesh.isManifold()) {
    throw std::runtime_error("handling of nonmanifold mesh not yet implemented");
  }

  // Initialize
  VertexData<double> distances(mesh, std::numeric_limits<double>::infinity());
  VertexData<int> signs(mesh, 1);
  VertexData<char> finalized(mesh, false);
  VertexData<char> isSource(mesh, false);
  VertexData<int> sourceLabels(mesh, -1);
  EdgeData<char> isSourceBaseEdge(mesh, false);
  EdgeData<int> sourceBaseLabels(mesh, -1);
  int nextSourceLabel = 0;

  auto cmp = [&signs](Entry left, Entry right) {
    Vertex leftV = std::get<1>(left);
    Vertex rightV = std::get<1>(right);
    if (signs[leftV] != signs[rightV]) {
      // We're looking at an edge that intersects the source, which may not have initial distance 0.
      // I *think* this can only happen if the positive vertex lies on the source, in which case it's the "closer" one.
      return signs[leftV] != 1;
    }
    return signs[leftV] * std::get<0>(left) > signs[rightV] * std::get<0>(right);
  };
  std::priority_queue<Entry, std::vector<Entry>, decltype(cmp)> frontierPQ(cmp);

  auto candidateImprovesDistance = [&](Vertex v, double candidateDist) {
    if (!std::isfinite(candidateDist)) return false;
    if (!std::isfinite(distances[v])) return true;
    if (signs[v] == 0) return std::abs(candidateDist) < std::abs(distances[v]);
    return signs[v] * candidateDist < signs[v] * distances[v];
  };

  auto addCandidateDistance = [&](Vertex v, double candidateDist, int sourceLabel) {
    if (!candidateImprovesDistance(v, candidateDist)) return;
    distances[v] = candidateDist;
    sourceLabels[v] = sourceLabel;
    frontierPQ.push(std::make_tuple(candidateDist, v, sourceLabel));
  };

  auto markSourceBaseEdge = [&](Edge e, int sourceLabel) {
    isSourceBaseEdge[e] = true;
    if (sourceBaseLabels[e] == -1) sourceBaseLabels[e] = sourceLabel;
  };

  // Initialize signs
  if (sign) {
    for (Vertex v : mesh.vertices()) signs[v] = 0;
    for (auto& curve : initialDistances) {
      size_t nNodes = curve.size();
      for (size_t i = 0; i < nNodes - 1; i++) {
        const SurfacePoint& pA = curve[i].first;
        const SurfacePoint& pB = curve[i + 1].first;
        Edge commonEdge = sharedEdge(pA, pB);
        if (commonEdge != Edge()) {
          // Assign +/- signs to the "third" vertices of each face straddling this edge.
          // These vertices might themselves lie on the curve, in which case we overwrite them below.
          Halfedge he = commonEdge.halfedge();
          signs[he.next().tipVertex()] = (he.vertex() == pA.vertex) ? -1 : 1;
          if (!commonEdge.isBoundary()) signs[he.twin().next().tipVertex()] = -signs[he.next().tipVertex()];
        } else {
          Face commonFace = sharedFace(pA, pB);
          if (commonFace == Face()) {
            throw std::logic_error("For signed fast marching distance, each curve segment must share a common face.");
          }
          BarycentricVector tangent(pA, pB);
          BarycentricVector normal = -tangent.rotate90(geometry);
          for (Vertex v : commonFace.adjacentVertices()) {
            BarycentricVector u(pA, SurfacePoint(v));
            signs[v] = (dot(geometry, normal, u) > 0) ? 1 : -1;
          }
        }
      }
    }
    for (auto& curve : initialDistances) {
      size_t nNodes = curve.size();
      for (size_t i = 0; i < nNodes - 1; i++) {
        Edge commonEdge = sharedEdge(curve[i].first, curve[i + 1].first);
        if (commonEdge != Edge()) markSourceBaseEdge(commonEdge, nextSourceLabel++);
      }
    }
    // Vertices on the curve are always assumed to have positive sign.
    for (auto& curve : initialDistances) {
      size_t nNodes = curve.size();
      for (size_t i = 0; i < nNodes; i++) {
        const SurfacePoint& p = curve[i].first;
        if (p.type != SurfacePointType::Vertex) continue;
        signs[p.vertex] = 1;
      }
    }
    // Fill in the signs of faces around the fan of any vertices.
    for (auto& curve : initialDistances) {
      size_t nNodes = curve.size();
      for (size_t i = 0; i < nNodes; i++) {
        const SurfacePoint& p = curve[i].first;
        if (p.type != SurfacePointType::Vertex) continue;
        Halfedge startHe = p.vertex.halfedge();
        Halfedge currHe = startHe;
        while (true) {
          if (signs[currHe.tipVertex()] != 0 && signs[currHe.tipVertex()] != signs[p.vertex] &&
              signs[currHe.next().tipVertex()] == 0) {
            signs[currHe.next().tipVertex()] = -1;
          }
          currHe = currHe.next().next().twin();
          if (currHe == startHe) break;
        }
      }
    }
  }
  // Initialize distances.
  for (auto& curve : initialDistances) {
    for (auto& x : curve) {
      const SurfacePoint& p = x.first;
      int sourceLabel = nextSourceLabel++;
      switch (p.type) {
      case (SurfacePointType::Vertex): {
        addCandidateDistance(p.vertex, x.second, sourceLabel);
        isSource[p.vertex] = true;
        break;
      }
      case (SurfacePointType::Edge): {
        const Vertex& vA = p.edge.firstVertex();
        const Vertex& vB = p.edge.secondVertex();
        double l = geometry.edgeLengths[p.edge];
        addCandidateDistance(vA, x.second + signs[vA] * p.tEdge * l, sourceLabel);
        addCandidateDistance(vB, x.second + signs[vB] * (1. - p.tEdge) * l, sourceLabel);
        isSource[vA] = true;
        isSource[vB] = true;
        markSourceBaseEdge(p.edge, sourceLabel);
        break;
      }
      case (SurfacePointType::Face): {
        Halfedge he = p.face.halfedge();
        const Vertex& vA = he.vertex();
        const Vertex& vB = he.next().vertex();
        const Vertex& vC = he.next().next().vertex();
        double lAB = geometry.edgeLengths[he.edge()];
        double lBC = geometry.edgeLengths[he.next().edge()];
        double lCA = geometry.edgeLengths[he.next().next().edge()];
        double lAB2 = lAB * lAB;
        double lBC2 = lBC * lBC;
        double lCA2 = lCA * lCA;
        double u = p.faceCoords[0];
        double v = p.faceCoords[1];
        double w = p.faceCoords[2];
        double dist2_A = lAB2 * (v * (1. - u)) + lCA2 * (w * (1. - u)) - lBC2 * v * w; // squared distance from p to vA
        double dist2_B = lAB2 * (u * (1. - v)) + lBC2 * (w * (1. - v)) - lCA2 * u * w; // squared distance from p to vB
        double dist2_C = lCA2 * (u * (1. - w)) + lBC2 * (v * (1. - w)) - lAB2 * u * v; // squared distance from p to vC
        addCandidateDistance(vA, x.second + signs[vA] * std::sqrt(std::max(0., dist2_A)), sourceLabel);
        addCandidateDistance(vB, x.second + signs[vB] * std::sqrt(std::max(0., dist2_B)), sourceLabel);
        addCandidateDistance(vC, x.second + signs[vC] * std::sqrt(std::max(0., dist2_C)), sourceLabel);
        isSource[vA] = true;
        isSource[vB] = true;
        isSource[vC] = true;
        markSourceBaseEdge(he.edge(), sourceLabel);
        markSourceBaseEdge(he.next().edge(), sourceLabel);
        markSourceBaseEdge(he.next().next().edge(), sourceLabel);
        break;
      }
      }
    }
  }

  size_t nFound = 0;
  size_t nVert = mesh.nVertices();

  // Search
  while (nFound < nVert && !frontierPQ.empty()) {

    // Pop the nearest element
    Entry currPair = frontierPQ.top();
    frontierPQ.pop();
    double currDist = std::get<0>(currPair);
    Vertex currV = std::get<1>(currPair);
    int currSourceLabel = std::get<2>(currPair);


    // Accept it if not stale
    if (finalized[currV] || currDist != distances[currV] || currSourceLabel != sourceLabels[currV]) {
      continue;
    }
    distances[currV] = currDist;
    finalized[currV] = true;
    nFound++;

    // Add any eligible neighbors
    for (Halfedge he : currV.incomingHalfedges()) {
      Vertex neighVert = he.vertex();

      // Add with length
      if (!finalized[neighVert] && !isSource[neighVert]) {
        if (signs[neighVert] == 0) signs[neighVert] = signs[currV];
        double newDist = currDist + signs[neighVert] * geometry.edgeLengths[he.edge()];
        addCandidateDistance(neighVert, newDist, currSourceLabel);
        continue;
      }

      bool useSourceBaseEdge = sourceLabels[neighVert] != currSourceLabel && isSource[currV] && isSource[neighVert] &&
                               isSourceBaseEdge[he.edge()];
      if (sourceLabels[neighVert] != currSourceLabel && !useSourceBaseEdge) continue;
      int updateSourceLabel = useSourceBaseEdge ? sourceBaseLabels[he.edge()] : currSourceLabel;

      // Check the third point of the "left" triangle straddling this edge
      if (he.isInterior()) {
        Vertex newVert = he.next().next().vertex();
        if (!finalized[newVert] && !isSource[newVert]) {

          // Compute the distance
          double lenB = geometry.edgeLengths[he.next().next().edge()];
          double distB = currDist;
          double lenA = geometry.edgeLengths[he.next().edge()];
          double distA = distances[neighVert];
          double theta = geometry.cornerAngles[he.next().next().corner()];
          if (signs[newVert] == 0) signs[newVert] = (signs[currV] != 0) ? signs[currV] : signs[he.next().vertex()];
          double newDist = eikonalDistanceSubroutine(lenA, lenB, theta, distA, distB, signs[newVert]);
          addCandidateDistance(newVert, newDist, updateSourceLabel);
        }
      }

      // Check the third point of the "right" triangle straddling this edge
      Halfedge heT = he.twin();
      if (heT.isInterior()) {
        Vertex newVert = heT.next().next().vertex();
        if (!finalized[newVert] && !isSource[newVert]) {

          // Compute the distance
          double lenB = geometry.edgeLengths[heT.next().edge()];
          double distB = currDist;
          double lenA = geometry.edgeLengths[heT.next().next().edge()];
          double distA = distances[neighVert];
          double theta = geometry.cornerAngles[heT.next().next().corner()];
          if (signs[newVert] == 0) signs[newVert] = (signs[currV] != 0) ? signs[currV] : signs[he.next().vertex()];
          double newDist = eikonalDistanceSubroutine(lenA, lenB, theta, distA, distB, signs[newVert]);
          addCandidateDistance(newVert, newDist, updateSourceLabel);
        }
      }
    }
  }

  return distances;
}


VertexData<double> FMMDistance(IntrinsicGeometryInterface& geometry,
                               const std::vector<std::pair<Vertex, double>>& initialDistances, bool sign) {

  std::vector<std::vector<std::pair<SurfacePoint, double>>> initialConditions;
  initialConditions.emplace_back();
  for (const auto& x : initialDistances) {
    initialConditions.back().emplace_back(SurfacePoint(x.first), x.second);
  }
  VertexData<double> distances = FMMDistance(geometry, initialConditions, sign);
  return distances;
}


} // namespace surface
} // namespace geometrycentral
