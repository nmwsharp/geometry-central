#pragma once

#include "geometrycentral/utilities/elementary_geometry.h"

namespace geometrycentral {

inline double triangleArea(double lAB, double lBC, double lCA) {
  double s = (lAB + lBC + lCA) / 2.0;
  double arg = std::max(0., s * (s - lAB) * (s - lBC) * (s - lCA));
  double area = std::sqrt(arg);
  return area;
}

inline bool triangleIsValid(double a, double b, double c, double tolerance) {
  return a + b >= c - tolerance && b + c >= a - tolerance && c + a >= b - tolerance;
}

inline Vector2 layoutTriangleVertex(const Vector2& pA, const Vector2& pB, const double& lBC, const double& lCA) {

  const double lAB = norm(pB - pA);
  double tArea = triangleArea(lAB, lBC, lCA);

  // Compute (signed) width and height of right triangle formed via altitude from C
  double h = 2. * tArea / lAB;
  double w = (lAB * lAB - lBC * lBC + lCA * lCA) / (2. * lAB);

  // Project some vectors to get the actual position
  Vector2 vABn = (pB - pA) / lAB;
  Vector2 vABnPerp{-vABn.y, vABn.x};
  Vector2 pC = pA + w * vABn + h * vABnPerp;

  return pC;
}

inline double pointLineSegmentDistance(Vector2 p, Vector2 lineA, Vector2 lineB) {
  double len2 = (lineA - lineB).norm2();
  if (len2 == 0.0) return (lineA - p).norm();                     // degenerate line case
  double t = clamp(dot(p - lineA, lineB - lineA) / len2, 0., 1.); // nearest point on segment
  Vector2 proj = lineA + t * (lineB - lineA);
  return (p - proj).norm();
}


inline SegmentSegmentIntersectionResult2D segmentSegmentIntersection(Vector2 pAStart, Vector2 pAEnd, Vector2 pBStart,
                                                                     Vector2 pBEnd) {

  // Lean on the ray version for now (can do better numerically)
  double lenA = norm(pAEnd - pAStart);
  RaySegmentIntersectionResult2D rayResult = raySegmentIntersection(pAStart, (pAEnd - pAStart) / lenA, pBStart, pBEnd);


  SegmentSegmentIntersectionResult2D result{rayResult.tRay / lenA, rayResult.tLine, false};
  result.hit = result.tA >= 0 && result.tA <= 1 && result.tB >= 0 && result.tB <= 1;
  return result;
}

inline RaySegmentIntersectionResult2D raySegmentIntersection(Vector2 rayStart, Vector2 rayDir, Vector2 lineStart,
                                                             Vector2 lineEnd) {

  Vector2 v1 = rayStart - lineStart;
  Vector2 v2 = lineEnd - lineStart;
  Vector2 v3{-rayDir.y, rayDir.x};

  double cross21 = cross(v2, v1);
  double tRay = cross21 / dot(v2, v3);
  double tLine = dot(v1, v3) / dot(v2, v3);

  return RaySegmentIntersectionResult2D{tRay, tLine};
}

inline RayRayIntersectionResult2D rayRayIntersection(Vector2 ray1Start, Vector2 ray1Dir, Vector2 ray2Start,
                                                     Vector2 ray2Dir) {

  Vector2 v1 = ray1Start - ray2Start;
  Vector2 v2 = ray2Dir;
  Vector2 v3{-ray1Dir.y, ray1Dir.x};

  double cross21 = cross(v2, v1);
  double tRay1 = cross21 / dot(v2, v3);
  double tRay2 = dot(v1, v3) / dot(v2, v3);

  return RayRayIntersectionResult2D{tRay1, tRay2};
}


inline double pointLineSegmentNeaestLocation(Vector3 p, Vector3 lineA, Vector3 lineB) {
  double len2 = (lineA - lineB).norm2();
  if (len2 == 0.0) return 0;                                      // degenerate line case
  double t = clamp(dot(p - lineA, lineB - lineA) / len2, 0., 1.); // nearest point on segment
  return t;
}
inline double pointLineSegmentDistance(Vector3 p, Vector3 lineA, Vector3 lineB) {
  double t = pointLineSegmentNeaestLocation(p, lineA, lineB);
  Vector3 proj = lineA + t * (lineB - lineA);
  return (p - proj).norm();
}

inline Vector3 tetFourthPoint(const Vector3& a, const Vector3& b, const Vector3& c, double Lad, double Lbd,
                              double Lcd) {
  // construct local orthonormal coordinate
  // system for the base triangle
  Vector3 e1, e2, e3;
  e1 = (b - a) / norm(b - a);
  e2 = (c - a);
  e2 -= dot(e2, e1) * e1;
  e2 = e2 / norm(e2);
  e3 = cross(e1, e2);

  // assuming this coordinate system is based
  // at the point a, get the nontrivial components
  // of the three points a,b,c, in this coordinate
  // system (i.e., we now have points (0,0), (x2,0),
  // and (x3,y3)
  double x2 = dot(b - a, e1);
  double x3 = dot(c - a, e1);
  double y3 = dot(c - a, e2);

  // compute the coordinates of the fourth point by
  // solving a system of three quadratic equations
  double x4 = (Lad * Lad - Lbd * Lbd + x2 * x2) / (2. * x2);
  double y4 = (Lad * Lad - Lcd * Lcd + x3 * x3 + y3 * y3 - 2. * x4 * x3) / (2. * y3);
  double z4 = sqrt(std::max(0., Lad * Lad - x4 * x4 - y4 * y4));

  // transform these coordinates back into the original coordinate system,
  // considering the two possible solutions (reflected across the base triangle)
  Vector3 d1 = a + x4 * e1 + y4 * e2 + z4 * e3;
  Vector3 d2 = a + x4 * e1 + y4 * e2 - z4 * e3;

  // compute the signed volume for the first solution
  double vol1 = dot(d1 - a, cross(b - a, c - a));

  // return the positively-oriented solution
  if (vol1 < 0.) {
    return d1;
  }
  return d2;
}

} // namespace geometrycentral
