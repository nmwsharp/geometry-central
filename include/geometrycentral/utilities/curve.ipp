#pragma once

#include <cmath>
#include <curve.h>
#include <iostream>

namespace geometrycentral {

// Curve =======================================================================

template <typename T>
T Curve<T>::operator()(double t) const {
  return value(t);
}

template <typename T>
T Curve<T>::operator()(double t, unsigned int k) const {
  return derivative(t, k);
}

// HermiteCurve ================================================================

template <typename T>
HermiteNode<T>::HermiteNode() : value(T::zero()), derivative(T::zero()) {}

template <typename T>
HermiteNode<T>::HermiteNode(const T& value_, const T& derivative_) : value(value_), derivative(derivative_) {}

template <typename T>
T HermiteCurve<T>::value(double t) const {
  size_t n = nodes.size();

  if (n == 0) return T();
  if (n == 1) return nodes[0].value + t * nodes[0].derivative;

  double L = length();
  if (t < 0.) return nodes[0].value + t * nodes[0].derivative;
  if (t > L) return nodes[n - 1].value + (t - static_cast<double>(n - 1)) * nodes[n - 1].derivative;

  double t0 = floor(t);
  size_t i = static_cast<size_t>(t0);

  const T& p0 = nodes[i + 0].value;
  const T& p1 = nodes[i + 1].value;

  const T& m0 = nodes[i + 0].derivative;
  const T& m1 = nodes[i + 1].derivative;

  double u = t - t0;
  double u2 = u * u;
  double u3 = u * u * u;

  return (2. * u3 - 3. * u2 + 1.) * p0 + (u3 - 2. * u2 + u) * m0 + (3. * u2 - 2. * u3) * p1 + (u3 - u2) * m1;
}

template <typename T>
T HermiteCurve<T>::derivative(double t, unsigned int k) const {
  if (k == 0) return HermiteCurve<T>::value(t);

  size_t n = nodes.size();

  if (n == 0) return T::zero();
  if (n == 1) {
    if (k == 1)
      return nodes[0].derivative;
    else
      return T::zero();
  }

  double L = length();
  if (t < 0.) {
    if (k == 1)
      return nodes[0].derivative;
    else
      return T::zero();
  }
  if (t > L) {
    if (k == 1)
      return nodes[n - 1].derivative;
    else
      return T::zero();
  }

  double t0 = floor(t);
  size_t i = static_cast<size_t>(t0);

  const T& p0 = nodes[i + 0].value;
  const T& p1 = nodes[i + 1].value;

  const T& m0 = nodes[i + 0].derivative;
  const T& m1 = nodes[i + 1].derivative;

  double u = t - t0;
  double u2 = u * u;
  double u3 = u * u * u;

  T dk;

  switch (k) {
  case 1:
    dk = (6. * u2 - 6. * u) * p0 + (3. * u2 - 4. * u + 1.) * m0 + (6. * u - 6. * u2) * p1 + (3. * u2 - 2. * u) * m1;
    break;
  case 2:
    dk = (12. * u - 6.) * p0 + (6. * u - 4.) * m0 + (6. - 12. * u) * p1 + (6. * u - 2.) * m1;
    break;
  case 3:
    dk = 12. * p0 + 6. * m0 + -12. * p1 + 6. * m1;
    break;
  default:
    dk = T::zero();
    break;
  }

  return dk;
}

template <typename T>
double HermiteCurve<T>::length() const {
  return static_cast<double>(nodes.size());
}

template <typename T>
void HermiteCurve<T>::addNode(const T& value, const T& derivative) {
  nodes.push_back(HermiteNode<T>(value, derivative));
}

template <typename T>
BezierCurve<T> HermiteCurve<T>::toBezier() const {
  throw std::invalid_argument("Hermite to Bezier conversion not yet implemented!");
}

// BezierCurve =================================================================

// TODO

} // namespace geometrycentral
