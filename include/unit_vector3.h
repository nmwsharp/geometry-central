#pragma once

#include <algorithm>
#include <cmath>

#include <utilities.h>
#include <vector3.h>

namespace geometrycentral {

class UnitVector3 : public Vector3 {
 public:
  UnitVector3(void) : Vector3{1., 0., 0.} {}
  UnitVector3(double x, double y, double z) : Vector3{x, y, z} { normalize(); }
  UnitVector3(const Vector3& v) : Vector3{v.x, v.y, v.z} { normalize(); }
  UnitVector3(const UnitVector3& u) : Vector3{u.x, u.y, u.z} { normalize(); }
  const UnitVector3& operator=(const Vector3& v) {
    x = v.x;
    y = v.y;
    z = v.z;
    normalize();
    return *this;
  }
  const UnitVector3& operator=(const UnitVector3& u) {
    x = u.x;
    y = u.y;
    z = u.z;
    return *this;
  }
};

inline double angle(UnitVector3& u0, UnitVector3& u1) {
  return ::acos(clamp(dot(u0, u1), -1., 1.));
}

UnitVector3 interpolate(UnitVector3& u0, UnitVector3& u1, double t);

}  // namespace geometrycentral