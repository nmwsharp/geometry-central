#pragma once

#include "geometrycentral/utilities/vector2.h"
#include "geometrycentral/utilities/vector3.h"


namespace geometrycentral {

// Computes the area of a triangle using Heron's rule. Will fail silently with 0 if inputs to not satisfy triangle
// inequality
double triangleArea(double lA, double lB, double lC);


// For triangle A,B,C, given vertex positions pA and pB, compute pC, such that lengths lBC and lCA are realized.
// pC will be placed on the side which gives CCW winding of ABC.
Vector2 layoutTriangleVertex(const Vector2& pA, const Vector2& pB, const double& lBC, const double& lCA);


double pointLineSegmentDistance(Vector2 p, Vector2 lineA, Vector2 lineB);

// === Line-line intersections

double pointLineSegmentDistance(Vector2 p, Vector2 lineA, Vector2 lineB);

// TODO for now, makes no special attempt at numerical robustness. In particular, may behave badly for colinear lines.
struct SegmentSegmentIntersectionResult2D {
  double tA;
  double tB;
  bool hit;
};
SegmentSegmentIntersectionResult2D segmentSegmentIntersection(Vector2 pAStart, Vector2 pAEnd, Vector2 pBStart,
                                                              Vector2 pBEnd);


// TODO for now, makes no special attempt at numerical robustness. In particular, may behave badly for colinear lines.
struct RaySegmentIntersectionResult2D {
  double tRay; // inf if no intersection
  double tLine;
};
RaySegmentIntersectionResult2D raySegmentIntersection(Vector2 rayStart, Vector2 rayDir, Vector2 lineStart,
                                                      Vector2 lineEnd);


} // namespace geometrycentral

#include "geometrycentral/utilities/elementary_geometry.ipp"
