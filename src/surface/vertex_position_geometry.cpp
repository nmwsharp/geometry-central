#include "geometrycentral/surface/vertex_position_geometry.h"

#include <fstream>
#include <limits>

namespace geometrycentral {
namespace surface {

VertexPositionGeometry::VertexPositionGeometry(HalfedgeMesh& mesh_)
    // clang-format off
// formatter is having an outburst here...
    : EmbeddedGeometryInterface(mesh_), inputVertexPositions(mesh_, Vector3{ 0., 0., 0, })
// clang-format on

{}

VertexPositionGeometry::VertexPositionGeometry(HalfedgeMesh& mesh_, VertexData<Vector3>& inputVertexPositions_)
    : EmbeddedGeometryInterface(mesh_), inputVertexPositions(inputVertexPositions_) {}


std::unique_ptr<VertexPositionGeometry> VertexPositionGeometry::copy() { return reinterpretTo(mesh); }

std::unique_ptr<VertexPositionGeometry> VertexPositionGeometry::reinterpretTo(HalfedgeMesh& targetMesh) {
  std::unique_ptr<VertexPositionGeometry> newGeom(new VertexPositionGeometry(targetMesh));
  newGeom->inputVertexPositions = inputVertexPositions.reinterpretTo(targetMesh);
  return newGeom;
}

void VertexPositionGeometry::computeVertexPositions() { vertexPositions = inputVertexPositions; }


} // namespace surface
} // namespace geometrycentral
