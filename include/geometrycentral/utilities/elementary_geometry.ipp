
namespace geometrycentral {

inline double triangleArea(double lAB, double lBC, double lCA) {
  double s = (lAB + lBC + lCA) / 2.0;
  double arg = std::max(0., s * (s - lAB) * (s - lBC) * (s - lCA));
  double area = std::sqrt(arg);
  return area;
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


} // namespace geometrycentral
