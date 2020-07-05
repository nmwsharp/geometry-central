#include "geometrycentral/surface/extrinsic_geometry_interface.h"

#include <limits>

namespace geometrycentral {
namespace surface {

// clang-format off
ExtrinsicGeometryInterface::ExtrinsicGeometryInterface(SurfaceMesh& mesh_) : 
  IntrinsicGeometryInterface(mesh_),

  edgeDihedralAnglesQ                    (&edgeDihedralAngles,                      std::bind(&ExtrinsicGeometryInterface::computeEdgeDihedralAngles, this),            quantities),
  vertexPrincipalCurvatureDirectionsQ    (&vertexPrincipalCurvatureDirections,      std::bind(&ExtrinsicGeometryInterface::computeVertexPrincipalCurvatureDirections, this),     quantities)
  
  {
  }
// clang-format on

// Edge dihedral angle
void ExtrinsicGeometryInterface::requireEdgeDihedralAngles() { edgeDihedralAnglesQ.require(); }
void ExtrinsicGeometryInterface::unrequireEdgeDihedralAngles() { edgeDihedralAnglesQ.unrequire(); }


void ExtrinsicGeometryInterface::computeVertexPrincipalCurvatureDirections() {
  edgeLengthsQ.ensureHave();
  halfedgeVectorsInVertexQ.ensureHave();
  edgeDihedralAnglesQ.ensureHave();

  vertexPrincipalCurvatureDirections = VertexData<Vector2>(mesh);

  for (Vertex v : mesh.vertices()) {
    Vector2 principalDir{0.0, 0.0};
    for (Halfedge he : v.outgoingHalfedges()) {
      double len = edgeLengths[he.edge()];
      double alpha = edgeDihedralAngles[he.edge()];
      Vector2 vec = halfedgeVectorsInVertex[he];
      principalDir += -vec * vec / len * std::abs(alpha);
    }

    vertexPrincipalCurvatureDirections[v] = principalDir / 4.0;
  }
}
void ExtrinsicGeometryInterface::requireVertexPrincipalCurvatureDirections() {
  vertexPrincipalCurvatureDirectionsQ.require();
}
void ExtrinsicGeometryInterface::unrequireVertexPrincipalCurvatureDirections() {
  vertexPrincipalCurvatureDirectionsQ.unrequire();
}


} // namespace surface
} // namespace geometrycentral
