#include "geometrycentral/surface/barycentric_vector.h"


namespace geometrycentral {
namespace surface {

BarycentricVector faceVectorRotate90(const BarycentricVector& w, IntrinsicGeometryInterface& geom) {
  // Found via (D(fj - fi) x (fk - fi))/(2A) x (Dw), where D is the matrix taking barycentric coords to 3D local coords,
  // where (fj - fi) x (fk - fi)/(2A) is the unit normal vector to the 3D-embedded local triangle. Then we apply D^{-1}
  // to map back to barycentric coords.
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

BarycentricVector faceVectorRotate(const BarycentricVector& w, IntrinsicGeometryInterface& geom, double angle) {
  // Found by applying the 3D axis-angle rotation matrix in 3D local coordinates (then transforming back to barycentric
  // coordinates), which conveniently avoids (almost all) square roots.
  double ui = w.faceCoords[0];
  double uj = w.faceCoords[1];
  double uk = w.faceCoords[2];
  geom.requireEdgeLengths();
  double l_ij = geom.edgeLengths[w.face.halfedge().edge()];
  double l_jk = geom.edgeLengths[w.face.halfedge().next().edge()];
  double l_ki = geom.edgeLengths[w.face.halfedge().next().next().edge()];
  geom.unrequireEdgeLengths();
  double lij2 = l_ij * l_ij;
  double ljk2 = l_jk * l_jk;
  double lki2 = l_ki * l_ki;
  double A2 = -(l_ij - l_jk - l_ki) * (l_ij + l_jk - l_ki) * (l_ij - l_jk + l_ki) * (l_ij + l_jk + l_ki) / 16.;
  double bA = 2. * std::sqrt(A2); // 2A
  double qA2 = 4. * A2;           // 4A^2
  double ai2 = 0.5 * (lij2 - ljk2 + lki2);
  double aj2 = 0.5 * (ljk2 - lki2 + lij2);
  double ak2 = 0.5 * (lki2 - lij2 + ljk2);
  double cosTheta = std::cos(angle);
  double sinTheta = std::sin(angle);
  double oneMinusCos = 1. - cosTheta;
  double aj2ak2 = aj2 * ak2;
  double ai2ak2 = ai2 * ak2;
  double ai2aj2 = ai2 * aj2;
  double xi = (qA2 * cosTheta + oneMinusCos * aj2ak2) * ui;
  double xj = (oneMinusCos * aj2ak2 - bA * sinTheta * aj2) * uj;
  double xk = (oneMinusCos * aj2ak2 + bA * sinTheta * ak2) * uk;
  double yi = (oneMinusCos * ai2ak2 + bA * sinTheta * ai2) * ui;
  double yj = (qA2 * cosTheta + oneMinusCos * ai2ak2) * uj;
  double yk = (oneMinusCos * ai2ak2 - bA * sinTheta * ak2) * uk;
  double zi = (oneMinusCos * ai2aj2 - bA * sinTheta * ai2) * ui;
  double zj = (oneMinusCos * ai2aj2 + bA * sinTheta * aj2) * uj;
  double zk = (qA2 * cosTheta + oneMinusCos * ai2aj2) * uk;
  Vector3 newCoords = {xi + xj + xk, yi + yj + yk, zi + zj + zk};
  newCoords /= qA2;
  return BarycentricVector(w.face, newCoords);
}

} // namespace surface
} // namespace geometrycentral
