#include "geometrycentral/surface/subdivide.h"

namespace geometrycentral {
namespace surface {

void linearSubdivide(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geo) {

  VertexData<Vector3>& pos = geo.inputVertexPositions;

  // Compute new positions for original vertices
  VertexData<Vector3> newPositions = geo.inputVertexPositions;

  FaceData<Vector3> splitFacePositions(mesh);
  for (Face f : mesh.faces()) {
    double D = (double)f.degree();
    splitFacePositions[f] = Vector3::zero();
    for (Vertex v : f.adjacentVertices()) {
      splitFacePositions[f] += geo.inputVertexPositions[v] / D;
    }
  }

  EdgeData<Vector3> splitEdgePositions(mesh);
  for (Edge e : mesh.edges()) {
    splitEdgePositions[e] = (pos[e.halfedge().tailVertex()] + pos[e.halfedge().tipVertex()]) / 2.;
  }

  // Subdivide edges
  VertexData<bool> isOrigVert(mesh, true);
  EdgeData<bool> isOrigEdge(mesh, true);
  for (Edge e : mesh.edges()) {
    if (!isOrigEdge[e]) continue;

    Vector3 newPos = splitEdgePositions[e];

    // split the edge
    Halfedge newHe = mesh.insertVertexAlongEdge(e);
    Vertex newV = newHe.vertex();

    isOrigVert[newV] = false;
    isOrigEdge[newHe.edge()] = false;
    isOrigEdge[newHe.twin().next().edge()] = false;

    GC_SAFETY_ASSERT(isOrigVert[newHe.twin().vertex()] && isOrigVert[newHe.twin().next().twin().vertex()], "???");

    newPositions[newV] = newPos;
  }

  // Subdivide faces
  FaceData<bool> isOrigFace(mesh, true);
  for (Face f : mesh.faces()) {
    if (!isOrigFace[f]) continue;

    Vector3 newPos = splitFacePositions[f];

    // split the face
    Vertex newV = mesh.insertVertex(f);
    isOrigVert[newV] = false;
    newPositions[newV] = newPos;

    for (Face f : newV.adjacentFaces()) {
      isOrigFace[f] = false;
    }

    // get incoming halfedges before we start deleting edges
    std::vector<Halfedge> incomingHalfedges;
    for (Halfedge he : newV.incomingHalfedges()) {
      incomingHalfedges.push_back(he);
    }

    for (Halfedge he : incomingHalfedges) {
      if (isOrigVert[he.vertex()]) {
        Face mergedF = mesh.removeEdge(he.edge());
        if (mergedF == Face()) {
          std::cout << "merge was impossible?" << std::endl;
        } else {
          isOrigFace[mergedF] = false;
        }
      }
    }
  }

  mesh.compress();
  geo.inputVertexPositions = newPositions;
  geo.refreshQuantities();
}

void catmullClarkSubdivide(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geo) {
  VertexData<Vector3>& pos = geo.inputVertexPositions;

  // Compute new positions for original vertices
  VertexData<Vector3> newPositions(mesh);

  FaceData<Vector3> splitFacePositions(mesh);
  for (Face f : mesh.faces()) {
    double D = (double)f.degree();
    splitFacePositions[f] = Vector3::zero();
    for (Vertex v : f.adjacentVertices()) {
      splitFacePositions[f] += geo.inputVertexPositions[v] / D;
    }
  }

  EdgeData<Vector3> splitEdgePositions(mesh);
  for (Edge e : mesh.edges()) {
    std::array<Face, 2> neigh{e.halfedge().face(), e.halfedge().twin().face()};
    splitEdgePositions[e] = (splitFacePositions[neigh[0]] + splitFacePositions[neigh[1]]) / 2.;
  }

  for (Vertex v : mesh.vertices()) {
    double D = (double)v.degree();
    Vector3 S = geo.inputVertexPositions[v];
    Vector3 R = Vector3::zero();
    for (Edge e : v.adjacentEdges()) {
      R += splitEdgePositions[e] / D;
    }
    Vector3 Q = Vector3::zero();
    for (Face f : v.adjacentFaces()) {
      Q += splitFacePositions[f] / D;
    }
    newPositions[v] = (Q + 2 * R + (D - 3) * S) / D;
  }

  // Subdivide edges
  VertexData<bool> isOrigVert(mesh, true);
  EdgeData<bool> isOrigEdge(mesh, true);
  for (Edge e : mesh.edges()) {
    if (!isOrigEdge[e]) continue;

    Vector3 newPos = splitEdgePositions[e];

    // split the edge
    Halfedge newHe = mesh.insertVertexAlongEdge(e);
    Vertex newV = newHe.vertex();

    isOrigVert[newV] = false;
    isOrigEdge[newHe.edge()] = false;
    isOrigEdge[newHe.twin().next().edge()] = false;

    GC_SAFETY_ASSERT(isOrigVert[newHe.twin().vertex()] && isOrigVert[newHe.twin().next().twin().vertex()], "???");

    newPositions[newV] = newPos;
  }

  // Subdivide faces
  FaceData<bool> isOrigFace(mesh, true);
  for (Face f : mesh.faces()) {
    if (!isOrigFace[f]) continue;

    Vector3 newPos = splitFacePositions[f];

    // split the face
    Vertex newV = mesh.insertVertex(f);
    isOrigVert[newV] = false;
    newPositions[newV] = newPos;

    for (Face f : newV.adjacentFaces()) {
      isOrigFace[f] = false;
    }

    // get incoming halfedges before we start deleting edges
    std::vector<Halfedge> incomingHalfedges;
    for (Halfedge he : newV.incomingHalfedges()) {
      incomingHalfedges.push_back(he);
    }

    for (Halfedge he : incomingHalfedges) {
      if (isOrigVert[he.vertex()]) {
        Face mergedF = mesh.removeEdge(he.edge());
        if (mergedF == Face()) {
          std::cout << "merge was impossible?" << std::endl;
        } else {
          isOrigFace[mergedF] = false;
        }
      }
    }
  }

  mesh.compress();
  geo.inputVertexPositions = newPositions;
  geo.refreshQuantities();
}

void loopSubdivide(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geo) {
  MutationManager mm(mesh, geo);
  return loopSubdivide(mesh, geo, mm);
}

void loopSubdivide(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geo, MutationManager& mm) {
  GC_SAFETY_ASSERT(mesh.isTriangular(), "Cannot run loop subdivision on a mesh with non-triangular faces");

  VertexData<Vector3>& pos = geo.inputVertexPositions;
  VertexData<bool> isOrigVert(mesh, true);
  EdgeData<bool> isOrigEdge(mesh, true);

  // Compute new positions for original vertices
  VertexData<Vector3> newPositions(mesh);
  for (Vertex v : mesh.vertices()) {
    double n = (double)v.degree();
    double u = (n == 3) ? 3. / 16. : 3. / (8. * n);
    newPositions[v] = (1 - n * u) * pos[v];
    for (Vertex w : v.adjacentVertices()) {
      newPositions[v] += u * pos[w];
    }
  }

  EdgeData<Vector3> splitVertexPositions(mesh);
  for (Edge e : mesh.edges()) {
    if (e.isBoundary()) {
      splitVertexPositions[e] = (pos[e.halfedge().tailVertex()] + pos[e.halfedge().tipVertex()]) / 2.;
    } else {
      std::array<Vertex, 2> ends{e.halfedge().tailVertex(), e.halfedge().tipVertex()};
      std::array<Vertex, 2> sides{e.halfedge().next().next().vertex(), e.halfedge().twin().next().next().vertex()};
      splitVertexPositions[e] = 3. / 8. * (pos[ends[0]] + pos[ends[1]]) + 1. / 8. * (pos[sides[0]] + pos[sides[1]]);
    }
  }

  // Subdivide edges
  std::vector<Edge> toFlip;
  for (Edge e : mesh.edges()) {
    if (!isOrigEdge[e]) continue;

    // gather both vertices incident on the edge
    Vertex oldA = e.halfedge().tipVertex();
    Vertex oldB = e.halfedge().tailVertex();

    Vector3 newPos = splitVertexPositions[e];

    // split the edge
    Vertex newV = mm.splitEdge(e, newPos).vertex();
    isOrigVert[newV] = false;
    newPositions[newV] = newPos;

    // iterate through the edges incident on the new vertex
    for (Edge e : newV.adjacentEdges()) {
      isOrigEdge[e] = false;               // mark the new edges
      Vertex otherV = e.otherVertex(newV); // other side of edge

      // if this is a new edge between an old an new vertex, save for
      // flipping
      if (isOrigVert[otherV] && otherV != oldA && otherV != oldB) {
        toFlip.push_back(e);
      }
    }
  }

  for (Edge e : toFlip) {
    mesh.flip(e);
  }

  mesh.compress();
  geo.inputVertexPositions = newPositions;
  geo.refreshQuantities();
}

} // namespace surface
} // namespace geometrycentral
