#include <spline_surface.h>

// BezierPatch =================================================================

template <>
Vector3 BezierPatch<Vector3>::position(const PatchParam& point) {
  return evaluate(point);
}

template <>
Vector3 BezierPatch<Vector3>::tangentU(const PatchParam& point) {
  double u = point.x;
  double v = point.y;
  double a = 1. - u;
  double b = 1. - v;

  // Evaluate cubic Bernstein polynomials in v
  double Bv0 = b * b * b;
  double Bv1 = 3. * v * b * b;
  double Bv2 = 3. * v * v * b;
  double Bv3 = v * v * v;

  // Evaluate first derivatives of Bernstein polynomials in u
  double dBu0 = -3. * a * a;
  double dBu1 = 3. * a * a - 6. * a * u;
  double dBu2 = 6. * a * u - 3. * u * u;
  double dBu3 = 3. * u * u;

  // Evaluate first u-derivative of bicubic interpolant
  array4x4<Vector3>& c = controlPoints;
  return dBu0 *
             (Bv0 * c[0][0] + Bv1 * c[0][1] + Bv2 * c[0][2] + Bv3 * c[0][3]) +
         dBu1 *
             (Bv0 * c[1][0] + Bv1 * c[1][1] + Bv2 * c[1][2] + Bv3 * c[1][3]) +
         dBu2 *
             (Bv0 * c[2][0] + Bv1 * c[2][1] + Bv2 * c[2][2] + Bv3 * c[2][3]) +
         dBu3 * (Bv0 * c[3][0] + Bv1 * c[3][1] + Bv2 * c[3][2] + Bv3 * c[3][3]);
}

template <>
Vector3 BezierPatch<Vector3>::tangentV(const PatchParam& point) {
  double u = point.x;
  double v = point.y;
  double a = 1. - u;
  double b = 1. - v;

  // Evaluate cubic Bernstein polynomials in u
  double Bu0 = a * a * a;
  double Bu1 = 3. * u * a * a;
  double Bu2 = 3. * u * u * a;
  double Bu3 = u * u * u;

  // Evaluate first derivatives of Bernstein polynomials in v
  double dBv0 = -3. * b * b;
  double dBv1 = 3. * b * b - 6. * b * v;
  double dBv2 = 6. * b * v - 3. * v * v;
  double dBv3 = 3. * v * v;

  // Evaluate first v-derivative of bicubic interpolant
  array4x4<Vector3>& c = controlPoints;
  return dBv0 *
             (Bu0 * c[0][0] + Bu1 * c[1][0] + Bu2 * c[2][0] + Bu3 * c[3][0]) +
         dBv1 *
             (Bu0 * c[0][1] + Bu1 * c[1][1] + Bu2 * c[2][1] + Bu3 * c[3][1]) +
         dBv2 *
             (Bu0 * c[0][2] + Bu1 * c[1][2] + Bu2 * c[2][2] + Bu3 * c[3][2]) +
         dBv3 * (Bu0 * c[0][3] + Bu1 * c[1][3] + Bu2 * c[2][3] + Bu3 * c[3][3]);
}

template <>
Vector3 BezierPatch<Vector3>::normal(const PatchParam& point) {
  return unit(cross(tangentU(point), tangentV(point)));
}

// GregoryPatchTri =============================================================

template <>
Vector3 GregoryPatchTri<Vector3>::position(const PatchParam& point) {
  return evaluate(point);
}

template <>
Vector3 GregoryPatchTri<Vector3>::tangentU(const PatchParam& point) {
  double u = point.x;
  double v = point.y;
  double w = point.z;

  // Evaluate face points, handling degeneracy at corners
  Vector3 F0 = (v + w == 0.) ? f0p : (w * f0m + v * f0p) / (v + w);
  Vector3 F1 = (w + u == 0.) ? f1p : (u * f1m + w * f1p) / (w + u);
  Vector3 F2 = (u + v == 0.) ? f2p : (v * f2m + u * f2p) / (u + v);

  // Evaluate u-derivative of triangular Gregory interpolant,
  // making the approximation that the face points are fixed
  // (This expression was automatically generated using Mathematica)
  return 3. *
         (u * (2. * e0m - 3. * e0m * u + p0 * u) - 3. * e2m * v + 4. * F2 * v +
          u * (-4. * e0m + 6. * e2m + 8. * F0 - 16. * F2 +
               3. * (e0m + e0p - e2m - 4. * F0 + 4. * F2) * u) *
              v +
          2. *
              (-e1p + 2. * e2m + 2. * F1 - 4. * F2 +
               (e0m + e0p + e1m + e1p - 2. * (e2m + 2. * (F0 + F1 - 2. * F2))) *
                   u) *
              ((v) * (v)) +
          (e1m + e1p - e2m - 4. * F1 + 4. * F2) * (v * v * v) -
          p2 * ((-1. + u + v) * (-1. + u + v)) -
          e2p * (-1. + v) * (-1. + u + v) * (-1. + 3. * u + v));
}

template <>
Vector3 GregoryPatchTri<Vector3>::tangentV(const PatchParam& point) {
  double u = point.x;
  double v = point.y;
  double w = point.z;

  // Evaluate face points, handling degeneracy at corners
  Vector3 F0 = (v + w == 0.) ? f0p : (w * f0m + v * f0p) / (v + w);
  Vector3 F1 = (w + u == 0.) ? f1p : (u * f1m + w * f1p) / (w + u);
  Vector3 F2 = (u + v == 0.) ? f2p : (v * f2m + u * f2p) / (u + v);

  // Evaluate v-derivative of triangular Gregory interpolant,
  // making the approximation that the face points are fixed
  // (This expression was automatically generated using Mathematica)
  return 3 *
         (u * (-3 * e2p + 4 * F2 - 2 * e0m * u + 4 * (e2p + F0 - 2 * F2) * u +
               (e0m + e0p - e2p - 4 * F0 + 4 * F2) * ((u) * (u))) +
          2 * e1p * v +
          2 * u *
              (-2 * e1p + 3 * e2p + 4 * F1 - 8 * F2 +
               (e0m + e0p + e1m + e1p - 2 * (e2p + 2 * (F0 + F1 - 2 * F2))) *
                   u) *
              v +
          (-3 * e1p + p1 + 3 * (e1m + e1p - e2p - 4 * F1 + 4 * F2) * u) *
              ((v) * (v)) -
          p2 * ((-1 + u + v) * (-1 + u + v)) -
          e2m * (-1 + u) * (-1 + u + v) * (-1 + u + 3 * v));
}

template <>
Vector3 GregoryPatchTri<Vector3>::normal(const PatchParam& point) {
  return unit(cross(tangentU(point), tangentV(point)));
}

// GregoryPatchQuad ============================================================

template <>
Vector3 GregoryPatchQuad<Vector3>::position(const PatchParam& point) {
  return evaluate(point);
}

template <>
Vector3 GregoryPatchQuad<Vector3>::tangentU(const PatchParam& point) {
  double u = point.x;
  double v = point.y;
  double a = 1. - u;
  double b = 1. - v;

  // Evaluate cubic Bernstein polynomials in u
  double Bu0 = a * a * a;
  double Bu1 = 3. * u * a * a;
  double Bu2 = 3. * u * u * a;
  double Bu3 = u * u * u;

  // Evaluate first derivatives of Bernstein polynomials in v
  float dBv0 = -3. * b * b;
  float dBv1 = 3. * b * b - 6. * b * v;
  float dBv2 = 6. * b * v - 3. * v * v;
  float dBv3 = 3. * v * v;

  // Evaluate face points
  Vector3 F0 = (u + v == 0.) ? f0p : (u * f0p + v * f0m) / (u + v);
  Vector3 F1 = (a + v == 0.) ? f1p : (a * f1m + v * f1p) / (a + v);
  Vector3 F2 = (a + b == 0.) ? f2p : (a * f2p + b * f2m) / (a + b);
  Vector3 F3 = (u + b == 0.) ? f3p : (u * f3m + b * f3p) / (u + b);

  // Evaluate derivative of bicubic Bézier in u direction (which only
  // approximates the true derivative of the rational Gregory interpolant)
  return dBv0 * (Bu0 * p0 + Bu1 * e0p + Bu2 * e1m + Bu3 * p1) +
         dBv1 * (Bu0 * e0m + Bu1 * F0 + Bu2 * F1 + Bu3 * e1p) +
         dBv2 * (Bu0 * e3p + Bu1 * F3 + Bu2 * F2 + Bu3 * e2m) +
         dBv3 * (Bu0 * p3 + Bu1 * e3m + Bu2 * e2p + Bu3 * p2);

  // // Alt: exact derivative of rational expression (computed in Mathematica)
  // return 9*(1 - u)*sqr(u)*(-(e1m*sqr(-1 + v)) + e2p*sqr(v) + (f1m*sqr(-1 +
  // v)*v)/(1 - u + v) + ((f2p*(-1 + u) + f2m*(-1 + v))*(-1 + v)*sqr(v))/sqr(-2
  // + u + v) - (2*(f2p*(-1 + u) + f2m*(-1 + v))*(-1 + v)*v)/(-2 + u + v) -
  // ((f2p*(-1 + u) + f2m*(-1 + v))*sqr(v))/(-2 + u + v) - (f2m*(-1 +
  // v)*sqr(v))/(-2 + u + v) - (sqr(-1 + v)*v*(f1p - f1p*u + f1m*v))/sqr(1 - u +
  // v) + (sqr(-1 + v)*(f1p - f1p*u + f1m*v))/(1 - u + v) + (2*(-1 + v)*v*(f1p -
  // f1p*u + f1m*v))/(1 - u + v)) - 3*cube(-1 + u)*(-(p0*sqr(-1 + v)) + v*(2*e3p
  // - 3*e3p*v + p3*v) + e0m*(1 - 4*v + 3*sqr(v))) + 3*cube(u)*(-(p1*sqr(-1 +
  // v)) + v*(2*e2m - 3*e2m*v + p2*v) + e1p*(1 - 4*v + 3*sqr(v))) + 9*sqr(-1 +
  // u)*u*(-(e0p*sqr(-1 + v)) + (f0p*u*sqr(1 + u - v)*(-1 + v)*(2*sqr(v) + u*(-1
  // + 3*v)) + v*(f0m*sqr(1 + u - v)*(-1 + v)*(v*(-1 + 3*v) + u*(-2 + 4*v)) +
  // sqr(u + v)*(f3p*u*(u*(2 - 3*v) + 2*sqr(-1 + v)) + e3m*sqr(1 + u - v)*v +
  // f3m*(-1 + v)*(-2 + 5*v - 3*sqr(v) + u*(-2 + 4*v)))))/(sqr(1 + u - v)*sqr(u
  // + v)));
}

template <>
Vector3 GregoryPatchQuad<Vector3>::tangentV(const PatchParam& point) {
  double u = point.x;
  double v = point.y;
  double a = 1. - u;
  double b = 1. - v;

  // Evaluate cubic Bernstein polynomials in v
  double Bv0 = b * b * b;
  double Bv1 = 3. * v * b * b;
  double Bv2 = 3. * v * v * b;
  double Bv3 = v * v * v;

  // Evaluate first derivatives of Bernstein polynomials in u
  float dBu0 = -3. * a * a;
  float dBu1 = 3. * a * a - 6. * a * u;
  float dBu2 = 6. * a * u - 3. * u * u;
  float dBu3 = 3. * u * u;

  // Evaluate face points
  Vector3 F0 = (u + v == 0.) ? f0p : (u * f0p + v * f0m) / (u + v);
  Vector3 F1 = (a + v == 0.) ? f1p : (a * f1m + v * f1p) / (a + v);
  Vector3 F2 = (a + b == 0.) ? f2p : (a * f2p + b * f2m) / (a + b);
  Vector3 F3 = (u + b == 0.) ? f3p : (u * f3m + b * f3p) / (u + b);

  // Evaluate derivative of bicubic Bézier in u direction (which only
  // approximates the true derivative of the rational Gregory interpolant)
  return dBu0 * (Bv0 * p0 + Bv1 * e0m + Bv2 * e3p + Bv3 * p3) +
         dBu1 * (Bv0 * e0p + Bv1 * F0 + Bv2 * F3 + Bv3 * e3m) +
         dBu2 * (Bv0 * e1m + Bv1 * F1 + Bv2 * F2 + Bv3 * e2p) +
         dBu3 * (Bv0 * p1 + Bv1 * e1p + Bv2 * e2m + Bv3 * p2);

  // // Alt: exact derivative of rational expression (computed in Mathematica)
  // return 3.*((3.*(-1. + u)*sqr(u)*sqr(-1. + v)*sqr(v)*(-((f2m - f2p)*sqr(1. -
  // u + v)) - f1m*sqr(-2. + u + v) + f1p*sqr(-2. + u + v)))/(sqr(1. - u +
  // v)*sqr(-2. + u + v)) - (3.*sqr(-1. + u)*u*sqr(-1. + v)*sqr(v)*(f0m*sqr(1. +
  // u - v) - f0p*sqr(1. + u - v) + (f3m - f3p)*sqr(u + v)))/(sqr(1. + u -
  // v)*sqr(u + v)) + 2.*(1. - u)*u*(-(e1m*cube(-1. + v)) + e2p*cube(v) -
  // (3.*(f2p*(-1. + u) + f2m*(-1. + v))*(-1. + v)*sqr(v))/(-2. + u + v) +
  // (3.*sqr(-1. + v)*v*(f1p - f1p*u + f1m*v))/(1. - u + v)) -
  // sqr(u)*(-(e1m*cube(-1. + v)) + e2p*cube(v) - (3.*(f2p*(-1. + u) + f2m*(-1.
  // + v))*(-1. + v)*sqr(v))/(-2. + u + v) + (3.*sqr(-1. + v)*v*(f1p - f1p*u +
  // f1m*v))/(1. - u + v)) + sqr(-1. + u)*(-(e0p*cube(-1. + v)) + e3m*cube(v) +
  // (3.*sqr(-1. + v)*v*(f0p*u + f0m*v))/(u + v) - (3.*(-1. + v)*sqr(v)*(f3m +
  // f3p*u - f3m*v))/(1. + u - v)) - 2.*(1. - u)*u*(-(e0p*cube(-1. + v)) +
  // e3m*cube(v) + (3.*sqr(-1. + v)*v*(f0p*u + f0m*v))/(u + v) - (3.*(-1. +
  // v)*sqr(v)*(f3m + f3p*u - f3m*v))/(1. + u - v)) + sqr(u)*(-(p1*cube(-1. +
  // v)) + v*(3.*e1p*sqr(-1. + v) + v*(3.*e2m - 3.*e2m*v + p2*v))) - sqr(-1. +
  // u)*(-(p0*cube(-1. + v)) + v*(3.*e0m*sqr(-1. + v) + v*(3.*e3p - 3.*e3p*v +
  // p3*v))));
}

template <>
Vector3 GregoryPatchQuad<Vector3>::normal(const PatchParam& point) {
  return unit(cross(tangentU(point), tangentV(point)));
}

// SplineSurface ===============================================================

// Differential geometry (for T=Vector3 only) --- these methods just find the
// patch containing the query point and then call the correpsonding method from
// the patch.
template <>
Vector3 SplineSurface<Vector3>::position(const PatchPoint& p) {
  return evaluate(p);
}
template <>
Vector3 SplineSurface<Vector3>::normal(const PatchPoint& p) {
  SplinePatchPtr<Vector3> patch = getPatch(p.face);
  return patch->normal(p.coords);
}
template <>
Vector3 SplineSurface<Vector3>::tangentU(const PatchPoint& p) {
  SplinePatchPtr<Vector3> patch = getPatch(p.face);
  return patch->tangentU(p.coords);
}
template <>
Vector3 SplineSurface<Vector3>::tangentV(const PatchPoint& p) {
  SplinePatchPtr<Vector3> patch = getPatch(p.face);
  return patch->tangentV(p.coords);
}
// template<> Vector3 SplineSurface<Vector3>::pushForward(const SplineTangent&
// tangent) { SplinePatchPtr<Vector3> patch = getPatch(p.face); return
// patch->tangentV(p.coords); }
