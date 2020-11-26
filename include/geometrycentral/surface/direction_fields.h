#pragma once

#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/extrinsic_geometry_interface.h"
#include "geometrycentral/surface/embedded_geometry_interface.h"


namespace geometrycentral {
namespace surface {

// Compute useful geometric quantities relating to optimal (Levi-Civita) transport

// === Compute smoothest direction fields

// Smoothest unit-norm direction field
VertexData<Vector2> computeSmoothestVertexDirectionField(IntrinsicGeometryInterface& geometry, int nSym = 1);

// Like above, but with Dirichlet boundary conditions to align to hte boundary
VertexData<Vector2> computeSmoothestBoundaryAlignedVertexDirectionField(IntrinsicGeometryInterface& geometry, int nSym = 1);

FaceData<Vector2> computeSmoothestFaceDirectionField(IntrinsicGeometryInterface& geometry, int nSym = 1);

FaceData<Vector2> computeSmoothestBoundaryAlignedFaceDirectionField(IntrinsicGeometryInterface& geometry, int nSym = 1);

// Same as above, but aligned to curvatures
VertexData<Vector2> computeCurvatureAlignedVertexDirectionField(ExtrinsicGeometryInterface& geometry, int nSym = 2);

FaceData<Vector2> computeCurvatureAlignedFaceDirectionField(EmbeddedGeometryInterface& geometry, int nSym = 2);


// Find singularities in direction fields
FaceData<int> computeFaceIndex(IntrinsicGeometryInterface& geometry, const VertexData<Vector2>& directionField,
                               int nSym = 1);
VertexData<int> computeVertexIndex(IntrinsicGeometryInterface& geometry, const FaceData<Vector2>& directionField,
                                   int nSym = 1);

} // namespace surface
} // namespace geometrycentral
