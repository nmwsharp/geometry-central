#include "geometrycentral/halfedge_mesh.h"

namespace geometrycentral {

unsigned int DualVertexPtr::degree() {
  unsigned int k = 0;

  for (DualEdgePtr e : adjacentEdges()) {
    k++;
  }

  return k;
}

}  // namespace geometrycentral