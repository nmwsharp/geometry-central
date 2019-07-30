#include "geometrycentral/utilities/unit_vector3.h"

namespace geometrycentral {

UnitVector3 interpolate(UnitVector3& u0, UnitVector3& u1, double t) {
  double theta = angle(u0, u1);
  Vector3 w = cross(u0, u1);

  Vector3 e1 = u0;
  Vector3 e2 = cross(w, e1);

  Vector3 r = cos(t * theta) * e1 + sin(t * theta) * e2;

  return UnitVector3(r);
}

} // namespace geometrycentral
