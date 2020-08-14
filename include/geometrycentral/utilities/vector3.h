#pragma once

#include "geometrycentral/utilities/utilities.h"

#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>

namespace geometrycentral {

// Note: this class avoids any constructors so that it is a POD type
struct Vector3 {
  // Components
  double x, y, z;


  static Vector3 zero() { return Vector3{0., 0., 0.}; }
  static Vector3 constant(double c) { return Vector3{c, c, c}; }
  static Vector3 infinity() {
    const double inf = ::std::numeric_limits<double>::infinity();
    return Vector3{inf, inf, inf};
  }
  static Vector3 undefined() {
    const double nan = ::std::numeric_limits<double>::quiet_NaN();
    return Vector3{nan, nan, nan};
  }

  // Access-by-index
  double& operator[](int index) { return (&x)[index]; }
  double operator[](int index) const { return (&x)[index]; };

  // Overloaded operators
  Vector3 operator+(const Vector3& v) const;
  Vector3 operator-(const Vector3& v) const;
  Vector3 operator*(double s) const;
  Vector3 operator/(double s) const;
  Vector3& operator+=(const Vector3& v);
  Vector3& operator-=(const Vector3& v);
  Vector3& operator*=(const double& s);
  Vector3& operator/=(const double& s);
  bool operator==(const Vector3& v) const;
  bool operator!=(const Vector3& v) const;
  const Vector3 operator-() const;

  // Other functions
  Vector3 rotateAround(Vector3 axis, double theta) const;
  Vector3 removeComponent(const Vector3& unitDir) const; // removes component in direction D
  std::array<Vector3, 2> buildTangentBasis() const;      // build a basis orthogonal to D (need not be unit already)
  Vector3 normalize() const;
  Vector3 normalizeCutoff(double mag = 0.) const;
  Vector3 unit() const;

  double norm() const;
  double norm2() const;

  bool isFinite() const;
  bool isDefined() const;
};

// Scalar multiplication
template <typename T>
Vector3 operator*(const T s, const Vector3& v);

// Printing
::std::ostream& operator<<(::std::ostream& output, const Vector3& v);
::std::istream& operator>>(::std::istream& intput, Vector3& v);

double norm(const Vector3& v);
double norm2(const Vector3& v);

Vector3 normalize(const Vector3& v);
Vector3 normalizeCutoff(const Vector3& v, double mag = 0.);
Vector3 unit(const Vector3& v);

Vector3 cross(const Vector3& u, const Vector3& v);
double angle(const Vector3& u, const Vector3& v);
double angleInPlane(const Vector3& u, const Vector3& v, const Vector3& normal);
double dot(const Vector3& u, const Vector3& v);
double sum(const Vector3& u);
bool isfinite(const Vector3& u); // break camel case rule to match std
bool isDefined(const Vector3& u);
Vector3 clamp(const Vector3& val, const Vector3& low, const Vector3& high);
Vector3 componentwiseMin(const Vector3& u, const Vector3& v);
Vector3 componentwiseMax(const Vector3& u, const Vector3& v);


} // namespace geometrycentral

namespace std {
template <>
struct hash<geometrycentral::Vector3> {
  std::size_t operator()(const geometrycentral::Vector3& v) const;
};

// overload for std string
std::string to_string(geometrycentral::Vector3 value);

} // namespace std

#include "vector3.ipp"
