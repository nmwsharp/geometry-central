#pragma once

// Curve =======================================================================

// Curve is an abstract base class defining the interface for
// parameterized curves interpolating data of any type T that
// overloads standard vector space operations (+,-,*) and has
// a method T T::zero() returning a zero value (additive identity).
// (Such types include, for instance, Vector2 and Vector3.)
//
// Curves may be evaluated either via named methods or the
// parenthesis operator.  For example,
//
//    Curve<double> curve;
//    double a = curve.value( 1.23 );
//    double b = curve( 1.23 );
//    // Both a and b now hold the value at time t=1.23
//
// Likewise, derivatives can be expressed via
//
//    double u = curve.derivative( 1.23, 2 );
//    double v = curve( 1.23, 2 );
//    // Both u and v now hold the first derivative at time t=1.23
//
// Curve evaluation returns valid results for values in the range
// [0,L], where L is the length returned by Curve::length().  Any
// behavior outside this range is undefined, and may vary depending
// on the curve type.  For instance,
//
//    double L = curve.length();
//    curve.value( 2.*L ); // behavior undefined!
//
// Beyond this basic interface nothing is assumed about the
// representation of the curve, though subclasses may expose
// additional information.  See below for more detail.
//

#include <vector>

namespace geometrycentral {

template <typename T>
class Curve {
public:
  // Evaluate the curve at time t
  virtual T value(double t) const = 0;
  T operator()(double t) const;

  // Evaluate the kth derivative of the curve at time t
  // (NOTE: may return NaN if curve is not k times differentiable)
  virtual T derivative(double t, unsigned int k) const = 0;
  T operator()(double t, unsigned int k) const;

  // The curve is defined between [0,length]
  virtual double length() const = 0;
};

template <typename T>
class BezierCurve;
template <typename T>
class HermiteCurve;

// Hermite =====================================================================

// A HermiteCurve encodes a curve as a composite cubic in Hermite
// form, i.e., a collection of points and tangents interpolated at
// segment endpoints.  (For the correpsonding Bézier form, see
// BezierCurve.) Each node is encoded as an instance of HermiteNode.
// A HermiteCurve with no nodes defined evaluates like a constant
// curve with value equal to the default value of type T.  A curve
// with only one node likewise behaves like an affine function
// interpolating the given value/derivative at t=0.  Two or more nodes
// yield a standard composite cubic Hermite curve, where consecutive
// nodes are separated by unit length (i.e., a curve with N nodes has
// length N-1).

template <typename T>
struct HermiteNode {
public:
  HermiteNode();
  HermiteNode(const T& value, const T& derivative);

  T value;
  T derivative;
};

template <typename T>
class HermiteCurve : public Curve<T> {
public:
  virtual T value(double t) const override;
  virtual T derivative(double t, unsigned int k) const override;
  virtual double length() const override;

  void addNode(const T& value, const T& derivative);

  BezierCurve<T> toBezier() const;

  std::vector<HermiteNode<T>> nodes;
};

// Bézier ======================================================================

template <typename T>
class BezierCurve : public Curve<T> {
public:
  virtual T value(double t) const override;
  virtual T derivative(double t, unsigned int k) const override;
  virtual double length() const override;

  HermiteCurve<T> toHermite() const;

  std::vector<T> nodes;
};

} // namespace geometrycentral

// Implementation
#include <curve.ipp>
