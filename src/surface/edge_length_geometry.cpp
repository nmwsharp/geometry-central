#include "geometrycentral/surface/edge_length_geometry.h"

#include <fstream>
#include <limits>

namespace geometrycentral {
namespace surface {

EdgeLengthGeometry::EdgeLengthGeometry(SurfaceMesh& mesh_)
    : IntrinsicGeometryInterface(mesh_), inputEdgeLengths(edgeLengths) {

  edgeLengths = EdgeData<double>(mesh_, 0.);

  // The input edge lengths share storage with edgeLengths, increment the required counter and make sure they
  // never get cleared
  requireEdgeLengths();
  edgeLengthsQ.clearable = false;
}

EdgeLengthGeometry::EdgeLengthGeometry(SurfaceMesh& mesh_, const EdgeData<double>& inputEdgeLengths_)
    : IntrinsicGeometryInterface(mesh_), inputEdgeLengths(edgeLengths) {

  edgeLengths = inputEdgeLengths_;

  // The input edge lengths share storage with edgeLengths, increment the required counter and make sure they
  // never get cleared
  requireEdgeLengths();
  edgeLengthsQ.clearable = false;
}

std::unique_ptr<EdgeLengthGeometry> EdgeLengthGeometry::copy() { return reinterpretTo(mesh); }

std::unique_ptr<EdgeLengthGeometry> EdgeLengthGeometry::reinterpretTo(SurfaceMesh& targetMesh) {
  std::unique_ptr<EdgeLengthGeometry> newGeom(new EdgeLengthGeometry(targetMesh));
  newGeom->inputEdgeLengths = inputEdgeLengths.reinterpretTo(targetMesh);
  return newGeom;
}


void EdgeLengthGeometry::computeEdgeLengths() {
  // The input edge lengthss share storage with edgeLengths, so this is a no-op
}

} // namespace surface
} // namespace geometrycentral
