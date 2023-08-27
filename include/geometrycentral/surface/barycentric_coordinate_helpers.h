#pragma once

#include "geometrycentral/surface/surface_mesh.h"

#include "geometrycentral/utilities/vector2.h"
#include "geometrycentral/utilities/vector3.h"

// Various utility arithmetic for working with barycentric coordinates
//
// All barycentric coordinates are assumed to be normalized (sum to 1) unless otherwise stated.

// Useful references:
//   - "Barycentric Coordinates in Olympiad Geometry"
//   https://s3.amazonaws.com/aops-cdn.artofproblemsolving.com/resources/articles/bary.pdf
//   - https://math.stackexchange.com/a/2385307 (handy barycentric answer to "Determine if a line segment passes
//   “through” a triangle...")

namespace geometrycentral {
namespace surface {

// Length of a vector
double displacementLength2(Vector3 displacement, Vector3 triangleLengths);
double displacementLength(Vector3 displacement, Vector3 triangleLengths);

// Convert between cartesian and barycentric coordinates
Vector3 cartesianVectorToBarycentric(const std::array<Vector2, 3>& vertCoords, Vector2 faceVec);
Vector2 barycentricDisplacementToCartesian(const std::array<Vector2, 3>& vertCoords, Vector3 baryVec);

// Normalize to sum to 1 (does nothing else)
Vector3 normalizeBarycentric(Vector3 baryCoords); // Allows not-normalized input

// Shift to sum to 0 (does nothing else)
// (note that displacements between normalized barycentric coords should sum to 0)
Vector3 normalizeBarycentricDisplacement(Vector3 baryVec);

// Ensure all coordinates are in [0,1] and sum to 1
Vector3 projectInsideTriangle(Vector3 baryCoords); // Allows not-normalized input

// Is the point inside the triangle?
bool isInsideTriangle(Vector3 baryCoords);

// The index of the halfedge in a triangular face,
int halfedgeIndexInTriangle(Halfedge he);

// Given barycentric coordinates defined by treating refHe.vertex() as the `x` coordinate, permute to the canonical
// coordinate ordering for the face, which has face.halfedge().vertex() as the `x` coordinate.
Vector3 permuteBarycentricToCanonical(Vector3 baryCoords, Halfedge refHe);

// Inverse of above
Vector3 permuteBarycentricFromCanonical(Vector3 baryCoords, Halfedge refHe);


} // namespace surface
} // namespace geometrycentral

#include "geometrycentral/surface/barycentric_coordinate_helpers.ipp"
