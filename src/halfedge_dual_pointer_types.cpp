#include <halfedge_mesh.h>

unsigned int DualVertexPtr::degree() {
  unsigned int k = 0;

  for (DualEdgePtr e : adjacentEdges()) {
    k++;
  }

  return k;
}
