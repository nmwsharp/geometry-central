#include "geometrycentral/surface/signed_heat_method.h"

namespace geometrycentral {
namespace surface {

SignedHeatMethodSolver::SignedHeatMethodSolver(IntrinsicGeometryInterface& geom_, double tCoef_)
    : tCoef(tCoef_), mesh(geom_.mesh), geom(geom_)

{
  geom.requireEdgeLengths();
  geom.requireVertexLumpedMassMatrix();

  // Compute mean edge length and set shortTime
  double meanEdgeLength = 0.;
  for (Edge e : mesh.edges()) {
    meanEdgeLength += geom.edgeLengths[e];
  }
  meanEdgeLength /= mesh.nEdges();
  shortTime = tCoef * meanEdgeLength * meanEdgeLength;

  // We always want the mass matrix
  massMat = geom.vertexLumpedMassMatrix;

  geom.unrequireVertexLumpedMassMatrix();
  geom.unrequireEdgeLengths();
}