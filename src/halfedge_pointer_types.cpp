#include <halfedge_mesh.h>

namespace geometrycentral {

unsigned int VertexPtr::degree() {
  unsigned int k = 0;

  for (EdgePtr e : adjacentEdges()) {
    k++;
  }

  return k;
}

std::vector<Triangle> FacePtr::triangulation() {
  // Get list of vertices
  std::vector<VertexPtr> vertices;
  unsigned int k = 0;
  for (VertexPtr v : adjacentVertices()) {
    vertices.push_back(v);
    k++;
  }

  // Construct alternating triangulation
  std::vector<Triangle> triangles(k - 2);
  unsigned int a = 0;
  unsigned int b = k - 1;
  unsigned int i = 0;
  while (b - a > 1) {
    if (i % 2 == 0) {
      triangles[i][0] = vertices[a];
      triangles[i][1] = vertices[a + 1];
      triangles[i][2] = vertices[b];
      a++;
    } else {
      triangles[i][0] = vertices[b - 1];
      triangles[i][1] = vertices[b];
      triangles[i][2] = vertices[a];
      b--;
    }
    i++;
  }

  if (k % 2 == 1) {
    Triangle tmp = triangles[k - 3];
    triangles[k - 3][0] = tmp[1];
    triangles[k - 3][1] = tmp[2];
    triangles[k - 3][2] = tmp[0];
  }

  return triangles;
}

}  // namespace geometrycentral