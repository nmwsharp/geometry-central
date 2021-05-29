#include "geometrycentral/surface/vertex_position_geometry.h"

#include <fstream>
#include <limits>

namespace geometrycentral {
namespace surface {

VertexPositionGeometry::VertexPositionGeometry(SurfaceMesh& mesh_)
    : EmbeddedGeometryInterface(mesh_), inputVertexPositions(vertexPositions) {

  vertexPositions = VertexData<Vector3>(mesh_, Vector3{0., 0., 0.});

  // The input vertex positions share storage with vertexPositions, incremented the required counter and make sure they
  // never get cleared
  requireVertexPositions();
  vertexPositionsQ.clearable = false;
}

VertexPositionGeometry::VertexPositionGeometry(SurfaceMesh& mesh_, const VertexData<Vector3>& inputVertexPositions_)
    : EmbeddedGeometryInterface(mesh_), inputVertexPositions(vertexPositions) {

  vertexPositions = inputVertexPositions_;

  // The input vertex positions share storage with vertexPositions, incremented the required counter and make sure they
  // never get cleared
  requireVertexPositions();
  vertexPositionsQ.clearable = false;
}


std::unique_ptr<VertexPositionGeometry> VertexPositionGeometry::copy() { return reinterpretTo(mesh); }

std::unique_ptr<VertexPositionGeometry> VertexPositionGeometry::reinterpretTo(SurfaceMesh& targetMesh) {
  std::unique_ptr<VertexPositionGeometry> newGeom(new VertexPositionGeometry(targetMesh));
  newGeom->inputVertexPositions = inputVertexPositions.reinterpretTo(targetMesh);
  return newGeom;
}

void VertexPositionGeometry::computeVertexPositions() {
  // The input vertex positions share storage with vertexPositions, so this is a no-op
}


} // namespace surface
} // namespace geometrycentral
