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

  geodesicEdgeLengths = EdgeData<double>(mesh);
  for (EdgePtr e : mesh->edges()) {
    geodesicEdgeLengths[e] = norm(vertexPositions[e.halfedge().vertex()] - vertexPositions[e.halfedge().twin().vertex()]);
  }

  buildDependencies();
}
  
EdgeLengthGeometry::~EdgeLengthGeometry() {};

// Set new edgelengths to define the geometry, immediately recalculating any quantities that have been required.
void EdgeLengthGeometry::update(EdgeData<double> edgeLengths) {
  geodesicEdgeLengths = edgeLengths;
  recomputeQuantities();
}

// === Quantity implementations

void EdgeLengthGeometry::computeEdgeLengths() { 
  edgeLengths = geodesicEdgeLengths; 
}

} // namespace geometrycentral
