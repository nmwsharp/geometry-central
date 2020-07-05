#include "geometrycentral/surface/edge_length_geometry.h"

#include <fstream>
#include <limits>

namespace geometrycentral {
namespace surface {

EdgeLengthGeometry::EdgeLengthGeometry(SurfaceMesh& mesh_)
    : IntrinsicGeometryInterface(mesh_), inputEdgeLengths(mesh_, 0.)
{}

EdgeLengthGeometry::EdgeLengthGeometry(SurfaceMesh& mesh_, EdgeData<double>& inputEdgeLengths_)
    : IntrinsicGeometryInterface(mesh_), inputEdgeLengths(inputEdgeLengths_) {}

std::unique_ptr<EdgeLengthGeometry> EdgeLengthGeometry::copy() { return reinterpretTo(mesh); }

std::unique_ptr<EdgeLengthGeometry> EdgeLengthGeometry::reinterpretTo(SurfaceMesh& targetMesh) {
  std::unique_ptr<EdgeLengthGeometry> newGeom(new EdgeLengthGeometry(targetMesh));
  newGeom->inputEdgeLengths = inputEdgeLengths.reinterpretTo(targetMesh);
  return newGeom;
}


void EdgeLengthGeometry::computeEdgeLengths() { edgeLengths = inputEdgeLengths; }

} // namespace surface
} // namespace geometrycentral
