#pragma once

#include <vector>

#include <geometry.h>

// Compute useful geometric quantities relating to optimal (Levi-Civita)
// transport

// Note: All of these functions implicitly follow the convention that angles in
// tangent space
// are measured against vertex.halfedge.

// === Completely compute direction fields
//     If the mesh has boundary, imposes dirichlet boundary conditions to
//     conform to the boundary.
//     Otherwise, computes the unit-norm solution
//     t \in [0,1] controls the strength of alignment with principal directions
VertexData<Complex> computeSmoothestDirectionField(
    Geometry<Euclidean>* geometry, int nSym = 1, bool alignCurvature = false);

// Find singularities in direction fields
FaceData<int> computeFaceIndex(Geometry<Euclidean>* geometry,
                               VertexData<Complex> directionField,
                               int nSym = 1);

// === Compute intermediate vaues useful for transport and direction fields

VertexData<double> computeAngleDefects(Geometry<Euclidean>* geometry);

HalfedgeData<double> computeRescaledHalfedgeAngles(
    Geometry<Euclidean>* geometry);
HalfedgeData<double> computeRescaledHalfedgeAngles(
    Geometry<Euclidean>* geometry, const VertexData<double>& angleDefects);

HalfedgeData<double> computeTransportAngles(Geometry<Euclidean>* geometry);
HalfedgeData<double> computeTransportAngles(
    Geometry<Euclidean>* geometry,
    const HalfedgeData<double>& rescaledHalfedgeAngles);

// === Conversion functions

VertexData<Vector3> convertTangentAnglesToR3Vectors(
    Geometry<Euclidean>* geometry, const VertexData<double>& tangentAngles);
VertexData<Vector3> convertComplexDirectionsToR3Vectors(
    Geometry<Euclidean>* geometry, const VertexData<Complex>& directionField,
    int nSym = 1);
