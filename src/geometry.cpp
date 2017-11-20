#include <limits>
#include <fstream>
#include <geometry.h>

template <>
void Geometry<Euclidean>::normalize() {
  // compute center of mass
  Vector3 cm;
  for (VertexPtr v : mesh.vertices()) {
    cm += position(v);
  }
  cm /= mesh.nVertices();

  // translate to origin and determine radius
  double rMax = 0;
  for (VertexPtr v : mesh.vertices()) {
    Vector3& p = position(v);
    p -= cm;
    rMax = std::max(rMax, norm(p));
  }

  // rescale to unit sphere
  for (VertexPtr v : mesh.vertices()) {
    position(v) /= rMax;
  }
}
