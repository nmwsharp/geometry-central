#include "geometrycentral/edge_length_geometry.h"

#include "geometrycentral/discrete_operators.h"

#include <fstream>
#include <limits>

namespace geometrycentral {

EdgeLengthGeometry::EdgeLengthGeometry(HalfedgeMesh* mesh_, EdgeData<double>& edgeLengths_)
    : IntrinsicGeometry(mesh_), geodesicEdgeLengths(edgeLengths_)

{
  buildDependencies();
}

EdgeLengthGeometry::EdgeLengthGeometry(HalfedgeMesh* mesh_, VertexData<Vector3>& vertexPositions)
    : IntrinsicGeometry(mesh_) {

  edgeLengths = EdgeData<double>(mesh);
  for (EdgePtr e : mesh->edges()) {
    edgeLengths[e] = norm(vertexPositions[e.halfedge().vertex()] - vertexPositions[e.halfedge().twin().vertex()]);
  }

  buildDependencies();
}


// === Quantity implementations

void EdgeLengthGeometry::computeEdgeLengths() { edgeLengths = geodesicEdgeLengths; }

} // namespace geometrycentral
