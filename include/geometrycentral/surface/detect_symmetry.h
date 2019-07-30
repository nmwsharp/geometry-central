#pragma once

#include "geometrycentral/surface/geometry.h"

#include <vector>

namespace geometrycentral {
namespace surface {

struct SymmetryResult {
  bool symmetryFound = false;                  // was a symmetry found? is this data valid?
  std::vector<Vertex> canonicalVertices;       // a representative entry from each
                                               // set of symmetry pairs
  VertexData<std::vector<Vertex>> symmetrySet; // for each unique vertex,
                                               // all others vertices that
                                               // are symmetry pairs
};

// Look for a symmetry about a mirror plane
SymmetryResult detectSymmetryMirror(Geometry<Euclidean>* geom, Vector3 planeNormal, Vector3 planePoint);

// Look for a rotational symmetry
SymmetryResult detectSymmetryRotation(Geometry<Euclidean>* geom, Vector3 rotAxis, Vector3 rotPoint, int nSym);

// Automatically search for the typical mirror and rotation symmetries about the
// shape center
// Returns any symmetry which is found.
SymmetryResult detectSymmetryAuto(Geometry<Euclidean>* geom);
SymmetryResult detectSymmetryAutoRotation(Geometry<Euclidean>* geom);
SymmetryResult detectSymmetryAutoMirror(Geometry<Euclidean>* geom); // Look for a symmetry about a mirror plane

// Look for symmetry which is mirrored over the y and z planes
SymmetryResult detectSymmetryDoubleMirror(Geometry<Euclidean>* geom);

} // namespace surface
} // namespace geometrycentral
