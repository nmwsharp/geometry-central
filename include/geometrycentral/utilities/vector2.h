#pragma once

#include "geometrycentral/utilities/vector3.h"

#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>

namespace geometrycentral {

// Note: this class avoids any constructors so that it is a POD type
struct Vector2 {
  double x, y;

  static Vector2 zero() { return Vector2{0., 0.}; }
  static Vector2 constant(double c) { return Vector2{c, c}; }
  static Vector2 fromAngle(double theta) { return Vector2{std::cos(theta), std::sin(theta)}; }
  static Vector2 fromComplex(std::complex<double> c) { return Vector2{c.real(), c.imag()}; }
  static Vector2 infinity() {
    const double inf = std::numeric_limits<double>::infinity();
    return Vector2{inf, inf};
  }

  static Vector2 undefined() {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    return Vector2{nan, nan};
  }


  // Access-by-index
  double& operator[](int index) { return (&x)[index]; }
  double operator[](int index) const { return (&x)[index]; };

  // Overloaded operators
  // Multiplication & division are in the sense of complex numbers
  Vector2 operator+(const Vector2& v) const;
  Vector2 operator-(const Vector2& v) const;
  Vector2 operator*(const Vector2& v) const;
  Vector2 operator/(const Vector2& v) const;
  Vector2 operator*(double s) const;
  Vector2 operator/(double s) const;
  Vector2& operator+=(const Vector2& other);
  Vector2& operator-=(const Vector2& other);
  Vector2& operator*=(const Vector2& other);
  Vector2& operator/=(const Vector2& other);
  Vector2& operator*=(const double& s);
  Vector2& operator/=(const double& s);
  bool operator==(const Vector2& v) const;
  bool operator!=(const Vector2& v) const;
  const Vector2 operator-() const;

  // Conversion to std::complex
  operator std::complex<double>() const;

  // The non-member functions below return a new object; they do not modify in-place.

  Vector2 normalize() const;
  Vector2 normalizeCutoff(double mag = 0.) const;
  Vector2 unit() const; // alias for normalize
  Vector2 rotate(double theta) const;
  Vector2 rotate90() const;

  // Complex functions
  Vector2 pow(double p) const;  // complex power
  Vector2 pow(Vector2 p) const; // complex to complex power
  Vector2 conj() const;
  Vector2 inv() const;

  double arg() const;
  double norm() const;
  double norm2() const;

  bool isFinite() const;
  bool isDefined() const;
};

template <typename T>
Vector2 operator*(const T s, const Vector2& v);

::std::ostream& operator<<(std::ostream& output, const Vector2& v);
::std::istream& operator<<(std::istream& input, Vector2& v);

// Notice that all of these functions return a new vector when applicable.

double arg(const Vector2& v);
double norm(const Vector2& v);
double norm2(const Vector2& v);

double angle(const Vector2& u, const Vector2& v);
double orientedAngle(const Vector2& u, const Vector2& v);
double dot(const Vector2& u, const Vector2& v);
double cross(const Vector2& u, const Vector2& v);
Vector3 cross3(const Vector2& u, const Vector2& v); // assumes arguments are in x-y plane

Vector2 unit(const Vector2& v);
Vector2 normalize(const Vector2& v);
Vector2 normalizeCutoff(const Vector2& v, double mag = 0.);
Vector2 clamp(const Vector2& val, const Vector2& low, const Vector2& high);

bool isfinite(const Vector2& u); // break camel case rule to match std
bool isDefined(const Vector2& u);
Vector2 componentwiseMin(const Vector2& u, const Vector2& v);
Vector2 componentwiseMax(const Vector2& u, const Vector2& v);
double sum(const Vector2& u);

} // namespace geometrycentral

namespace std {
template <>
struct hash<geometrycentral::Vector2> {
  std::size_t operator()(const geometrycentral::Vector2& v) const;
};

// overload for std string
std::string to_string(geometrycentral::Vector2 value);

} // namespace std


#include "geometrycentral/utilities/vector2.ipp"
