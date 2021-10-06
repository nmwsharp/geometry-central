#pragma once

#include "geometrycentral/utilities/vector2.h"
#include "geometrycentral/utilities/vector3.h"


namespace geometrycentral {


// ======================================================
// ======== 2 dimensional
// ======================================================

// Computes the area of a triangle using Heron's rule. Will fail silently with 0 if inputs to not satisfy triangle
// inequality
double triangleArea(double lA, double lB, double lC);


// For triangle A,B,C, given vertex positions pA and pB, compute pC, such that lengths lBC and lCA are realized.
// pC will be placed on the side which gives CCW winding of ABC.
Vector2 layoutTriangleVertex(const Vector2& pA, const Vector2& pB, const double& lBC, const double& lCA);


double pointLineSegmentDistance(Vector2 p, Vector2 lineA, Vector2 lineB);


// === Line-line intersections


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
struct RayRayIntersectionResult2D {
  double tRay1; // inf if no intersection
  double tRay2;
};
RayRayIntersectionResult2D rayRayIntersection(Vector2 ray1Start, Vector2 ray1Dir, Vector2 ray2Start, Vector2 ray2Dir);



// Returns true if pTest is inside the incircle of triangle pA-pB-pC. No snazzy exact predicates or numerical epsilons are used; this is just a simple det(A) < 0 test.
bool inCircleTest(Vector2 pA, Vector2 pB, Vector2 pC, Vector2 pTest);

// ======================================================
// ======== 3 dimensional
// ======================================================

// Compute the distance to a line segment, returning either distance, location of the nearest point along the line, or
// both.
double pointLineSegmentDistance(Vector3 p, Vector3 lineA, Vector3 lineB);
double pointLineSegmentNeaestLocation(Vector3 p, Vector3 lineA, Vector3 lineB);

} // namespace geometrycentral

#include "geometrycentral/utilities/elementary_geometry.ipp"
