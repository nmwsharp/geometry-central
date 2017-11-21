#pragma once

#include "geometrycentral/vector3.h"
#include <cmath>
#include <iostream>

namespace geometrycentral {

// Note: this class avoids any constructors so that it is a POD type
struct Vector2 {
  double x, y;

  // Access-by-index
  double& operator[](int index) { return (&x)[index]; }
  double operator[](int index) const { return (&x)[index]; };

  // Overloaded operators
  Vector2 operator+(const Vector2& v) const;
  Vector2 operator-(const Vector2& v) const;
  Vector2 operator*(double s) const;
  Vector2 operator/(double s) const;
  Vector2& operator+=(const Vector2& other);
  Vector2& operator-=(const Vector2& other);
  Vector2& operator*=(const double& s);
  Vector2& operator/=(const double& s);
  bool operator==(const Vector2& v) const;
  bool operator!=(const Vector2& v) const;
  const Vector2 operator-() const;
  void normalize(void);

  static Vector2 constant(double c) { return Vector2{c, c}; }

  static Vector2 zero(void) { return Vector2{0., 0.}; }

  static Vector2 infinity(void) {
    const double inf = std::numeric_limits<double>::infinity();
    return Vector2{inf, inf};
  }

  static Vector2 undefined(void) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    return Vector2{nan, nan};
  }

  bool isFinite() const;
  bool isDefined() const;
};

Vector2 operator*(const double s, const Vector2& v);

::std::ostream& operator<<(std::ostream& output, const Vector2& v);

double norm(const Vector2& v);
double norm2(const Vector2& v);
Vector2 unit(const Vector2& v);
double angle(const Vector2& u, const Vector2& v);
double dot(const Vector2& u, const Vector2& v);
Vector3 cross(const Vector2& u,
              const Vector2& v);  // assumes arguments are in x-y plane
bool isFinite(const Vector2& u);
Vector2 componentwiseMin(const Vector2& u, const Vector2& v);
Vector2 componentwiseMax(const Vector2& u, const Vector2& v);

}  // namespace geometrycentral

#include "geometrycentral/vector2.ipp"
