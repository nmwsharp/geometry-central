
namespace geometrycentral {

inline double triangleArea(double lAB, double lBC, double lCA) {
  double s = (lAB + lBC + lCA) / 2.0;
  double arg = std::max(0., s * (s - lAB) * (s - lBC) * (s - lCA));
  double area = std::sqrt(arg);
  return area;
}

inline Vector2 layoutTriangleVertex(const Vector2& pA, const Vector2& pB, const double& lBC, const double& lCA) {

  const double lAB = norm(pB - pA);
  double tArea = triangleArea(lAB, lBC, lCA);

  // Compute width and height of right triangle formed via altitude from C
  double h = 2. * tArea / lAB;
  double w = std::sqrt(std::max(0., lCA * lCA - h * h));

  // Take the closer of the positive and negative solutions
  if (lBC * lBC > (lAB * lAB + lCA * lCA)) w *= -1.0;

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



} // namespace geometrycentral
