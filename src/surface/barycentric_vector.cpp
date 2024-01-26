#include "geometrycentral/surface/barycentric_vector.h"


namespace geometrycentral {
namespace surface {

BarycentricVector faceVectorRotated90(const BarycentricVector& w, IntrinsicGeometryInterface& geom) {
  double ui = w.faceCoords[0];
  double uj = w.faceCoords[1];
  double uk = w.faceCoords[2];
  geom.requireEdgeLengths();
  double l_ij = geom.edgeLengths[w.face.halfedge().edge()];
  double l_jk = geom.edgeLengths[w.face.halfedge().next().edge()];
  double l_ki = geom.edgeLengths[w.face.halfedge().next().next().edge()];
  geom.unrequireEdgeLengths();
  // coefficients of matrix D taking barycentric coords to 3D local coords, squared and multiplied by 2
  double ai = l_ij * l_ij - l_jk * l_jk + l_ki * l_ki;
  double aj = l_jk * l_jk - l_ki * l_ki + l_ij * l_ij;
  double ak = l_ki * l_ki - l_ij * l_ij + l_jk * l_jk;
  double s = 0.5 * (l_ij + l_jk + l_ki);
  double A = std::sqrt(s * (s - l_ij) * (s - l_jk) * (s - l_ki));
  Vector3 newCoords = {ak * uk - aj * uj, ai * ui - ak * uk, aj * uj - ai * ui};
  newCoords /= (4. * A);
  return BarycentricVector(w.face, newCoords);
}

} // namespace surface
} // namespace geometrycentral
