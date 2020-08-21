#pragma once

// Quaternion represents an element of the quaternions, along with all the usual
// vectors space operations (addition, multiplication by scalars, etc.).  The
// Hamilton product is expressed using the * operator:
//
//    Quaternion p, q, r;
//    r = q * p;
//
// and conjugation is expressed using the method Quaternion::bar():
//
//    Quaternion q;
//    double normQSquared = -q.bar()*q;
//
// Individual components can be accessed in several ways: the real and imaginary
// parts can be accessed using the methods Quaternion::re() and
// Quaternion::im():
//
//   Quaternion q;
//   double  a = q.re();
//   Vector3 b = q.im();
//
// or by index:
//
//   Quaternion q;
//   double a  = q[0];
//   double bi = q[1];
//   double bj = q[2];
//   double bk = q[3];
//

#include "geometrycentral/utilities/vector3.h"
#include <ostream>

namespace geometrycentral {

class Quaternion {
public:
  Quaternion(void);
  // initializes all components to zero

  // Quaternion(const Quaternion& q); // should be same as implicitly defined
  // initializes from existing quaternion

  explicit Quaternion(double s, double vi, double vj, double vk);
  // initializes with specified real (s) and imaginary (v) components

  explicit Quaternion(double s, const Vector3& v);
  // initializes with specified real (s) and imaginary (v) components

  explicit Quaternion(double s);
  // initializes purely real quaternion with specified real (s) component
  // (imaginary part is zero)

  explicit Quaternion(const Vector3& v);
  // initializes purely imaginary quaternion with specified imaginary (v)
  // component (real part is zero)

  const Quaternion& operator=(double s);
  // assigns a purely real quaternion with real value s

  const Quaternion& operator=(const Vector3& v);
  // assigns a purely real quaternion with imaginary value v

  double& operator[](int index);
  // returns reference to the specified component (0-based indexing: r, i, j, k)

  const double& operator[](int index) const;
  // returns const reference to the specified component (0-based indexing: r, i,
  // j, k)

  void toMatrix(double Q[4][4]) const;
  // builds 4x4 matrix Q representing (left) quaternion multiplication

  double& re(void);
  // returns reference to double part

  const double& re(void) const;
  // returns const reference to double part

  Vector3& im(void);
  // returns reference to imaginary part

  const Vector3& im(void) const;
  // returns const reference to imaginary part

  Quaternion operator+(const Quaternion& q) const;
  // addition

  Quaternion operator-(const Quaternion& q) const;
  // subtraction

  Quaternion operator-(void) const;
  // negation

  Quaternion operator*(double c) const;
  // right scalar multiplication

  Quaternion operator/(double c) const;
  // scalar division

  void operator+=(const Quaternion& q);
  // addition / assignment

  void operator+=(double c);
  // addition / assignment of pure real

  void operator-=(const Quaternion& q);
  // subtraction / assignment

  void operator-=(double c);
  // subtraction / assignment of pure real

  void operator*=(double c);
  // scalar multiplication / assignment

  void operator/=(double c);
  // scalar division / assignment

  Quaternion operator*(const Quaternion& q) const;
  // Hamilton product

  void operator*=(const Quaternion& q);
  // Hamilton product / assignment

  Quaternion bar(void) const;
  // conjugation

  Quaternion inv(void) const;
  // inverse

  double norm(void) const;
  // returns Euclidean length

  double norm2(void) const;
  // returns Euclidean length squared

  Quaternion unit(void) const;
  // returns unit quaternion

  void normalize(void);
  // divides by Euclidean length

protected:
  double s;
  // scalar (double) part

  Vector3 v;
  // vector (imaginary) part
};

Quaternion conjugate(const Quaternion& q);
// conjugation---this method simply invokes Quaternion::bar(),
// and is needed only for providing a uniform interface to templated
// classes that need to handle both Quaternions and atomic types like
// double (specifically, the SparseMatrix and DenseMatrix classes)

Quaternion operator*(double c, const Quaternion& q);
// left scalar multiplication

std::ostream& operator<<(std::ostream& os, const Quaternion& q);
// prints components

} // namespace geometrycentral
