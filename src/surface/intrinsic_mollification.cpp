#include "geometrycentral/surface/intrinsic_mollification.h"


namespace geometrycentral {
namespace surface {

void mollifyIntrinsic(SurfaceMesh& mesh, EdgeData<double>& edgeLengths, double relativeFactor) {
  // Mean edge length
  double edgeSum = 0.;
  for (Edge e : mesh.edges()) {
    edgeSum += edgeLengths[e];
  }
  double meanEdge = edgeSum / mesh.nEdges();

  double mollifyDelta = meanEdge * relativeFactor;

  mollifyIntrinsicAbsolute(mesh, edgeLengths, mollifyDelta);
}

void mollifyIntrinsicAbsolute(SurfaceMesh& mesh, EdgeData<double>& edgeLengths, double absoluteFactor) {

  // Compute the mollify epsilon
  double mollifyEps = 0.;
  for (Halfedge he : mesh.interiorHalfedges()) {

    double lA = edgeLengths[he.edge()];
    double lB = edgeLengths[he.next().edge()];
    double lC = edgeLengths[he.next().next().edge()];

    double thisEPS = lC - lA - lB + absoluteFactor;
    mollifyEps = std::fmax(mollifyEps, thisEPS);
  }

  // Apply the offset
  for (Edge e : mesh.edges()) {
    edgeLengths[e] += mollifyEps;
  }
}

} // namespace surface
} // namespace geomtrycentral
