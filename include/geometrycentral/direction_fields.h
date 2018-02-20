#pragma once

#include <vector>

#include "geometrycentral/geometry.h"

// Compute useful geometric quantities relating to optimal (Levi-Civita)
// transport

// Note: All of these functions implicitly follow the convention that angles in
// tangent space
// are measured against vertex.halfedge.

namespace geometrycentral {

// === Completely compute direction fields
//     If the mesh has boundary, imposes dirichlet boundary conditions to
//     conform to the boundary.
//     Otherwise, computes the unit-norm solution
//     t \in [0,1] controls the strength of alignment with principal directions

VertexData<Complex> computeSmoothestVertexDirectionField(Geometry<Euclidean>* geometry, int nSym = 1,
                                                         bool alignCurvature = false);

FaceData<Complex> computeSmoothestFaceDirectionField(Geometry<Euclidean>* geometry, int nSym = 1,
                                                     bool alignCurvature = false);

// Find singularities in direction fields
FaceData<int> computeFaceIndex(Geometry<Euclidean>* geometry, VertexData<Complex> directionField, int nSym = 1);
VertexData<int> computeVertexIndex(Geometry<Euclidean>* geometry, FaceData<Complex> directionField, int nSym = 1);

} // namespace geometrycentral
