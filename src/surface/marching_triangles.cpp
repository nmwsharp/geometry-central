#include "geometrycentral/surface/marching_triangles.h"

namespace geometrycentral {
namespace surface {

namespace {

std::vector<std::vector<std::array<size_t, 2>>>
getCurveComponents(SurfaceMesh& mesh, const std::vector<SurfacePoint>& curveNodes,
                   const std::vector<std::array<size_t, 2>>& curveEdges) {

  // Note: This function assumes that `curveNodes` has been de-duplicated.
  std::vector<std::array<size_t, 2>> edgesToAdd = curveEdges;
  std::vector<std::vector<std::array<size_t, 2>>> curves;
  size_t nSegs = curveEdges.size();
  while (edgesToAdd.size() > 0) {
    std::array<size_t, 2> startSeg = edgesToAdd.back();
    edgesToAdd.pop_back();
    curves.emplace_back();
    std::vector<std::array<size_t, 2>>& currCurve = curves.back();
    currCurve.push_back(startSeg);

    // Add segs to the front end until we can't.
    std::array<size_t, 2> currSeg = startSeg;
    while (true) {
      const SurfacePoint& front = curveNodes[currSeg[1]];
      bool didWeFindOne = false;
      for (size_t i = 0; i < edgesToAdd.size(); i++) {
        std::array<size_t, 2> otherSeg = edgesToAdd[i];
        if (curveNodes[otherSeg[0]] == front) {
          currSeg = otherSeg;
          currCurve.push_back(otherSeg);
          edgesToAdd.erase(edgesToAdd.begin() + i);
          didWeFindOne = true;
          break;
        }
      }
      if (!didWeFindOne) break;
    }

    // Add segs to the back end until we can't.
    currSeg = startSeg;
    while (true) {
      const SurfacePoint& back = curveNodes[currSeg[0]];
      bool didWeFindOne = false;
      for (size_t i = 0; i < edgesToAdd.size(); i++) {
        std::array<size_t, 2> otherSeg = edgesToAdd[i];
        if (curveNodes[otherSeg[1]] == back) {
          currSeg = otherSeg;
          currCurve.insert(currCurve.begin(), otherSeg);
          edgesToAdd.erase(edgesToAdd.begin() + i);
          didWeFindOne = true;
          break;
        }
      }
      if (!didWeFindOne) break;
    }
  }
  return curves;
}
} // namespace

std::vector<std::vector<SurfacePoint>> marchingTriangles(IntrinsicGeometryInterface& geom, const VertexData<double>& u,
                                                         double isoval) {

  SurfaceMesh& mesh = *u.getMesh();
  std::vector<SurfacePoint> nodes;
  std::vector<std::array<size_t, 2>> edges;
  for (Face f : mesh.faces()) {
    std::vector<SurfacePoint> hits;
    BarycentricVector gradient(f);
    for (Halfedge he : f.adjacentHalfedges()) {
      // Record edge crossings
      Edge e = he.edge();
      Vertex v0 = e.firstVertex();
      Vertex v1 = e.secondVertex();
      double u0 = u[v0];
      double u1 = u[v1];
      double lB = std::min(u0, u1);
      double uB = std::max(u0, u1);
      if (lB == uB && lB == isoval) {
        hits.clear();
        hits = {SurfacePoint(v0), SurfacePoint(v1)};
        break;
      }
      double t = (isoval - lB) / (uB - lB);
      if (u0 > u1) t = 1. - t;
      if (t <= 1. && t >= 0.) hits.emplace_back(e, t);
      // Compute gradient of the scalar function.
      BarycentricVector heVec(he.next(), f);
      BarycentricVector ePerp = heVec.rotate90(geom);
      gradient += ePerp * u[he.vertex()];
    }
    if (hits.size() != 2) continue;

    // Orient segments so that smaller values are always on the "inside" of the curve.
    std::array<size_t, 2> seg;
    for (int i = 0; i < 2; i++) {
      auto iter = std::find(nodes.begin(), nodes.end(), hits[i]);
      if (iter != nodes.end()) {
        seg[i] = iter - nodes.begin();
      } else {
        nodes.push_back(hits[i]);
        seg[i] = nodes.size() - 1;
      }
    }
    BarycentricVector segTangent(hits[0], hits[1]);
    BarycentricVector segNormal = -segTangent.inFace(f).rotate90(geom);
    if (dot(geom, segNormal, gradient) > 0.) {
      edges.push_back(seg);
    } else {
      edges.push_back({seg[1], seg[0]});
    }
  }
  std::vector<std::vector<std::array<size_t, 2>>> components = getCurveComponents(mesh, nodes, edges);
  size_t nCurves = components.size();
  std::vector<std::vector<SurfacePoint>> curves(nCurves);
  for (size_t i = 0; i < nCurves; i++) {
    size_t nEdges = components[i].size();
    for (size_t j = 0; j < nEdges; j++) {
      curves[i].push_back(nodes[components[i][j][0]]);
    }
    curves[i].push_back(nodes[components[i][nEdges - 1][1]]);
  }
  return curves;
}

} // namespace surface
} // namespace geometrycentral