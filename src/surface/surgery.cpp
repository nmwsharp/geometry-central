#include "geometrycentral/surface/surgery.h"

#include <queue>

namespace geometrycentral {
namespace surface {

std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, HalfedgeData<Halfedge>> cutAlongEdges(ManifoldSurfaceMesh& origMesh,
                                                                                const EdgeData<char>& origCut) {

  // Create a copy of the input mesh
  std::unique_ptr<ManifoldSurfaceMesh> mesh = origMesh.copy();

  // Initialize parent references
  // The built-in dynamic container updates will automatically keep this container in sync as we modify the mesh
  HalfedgeData<Halfedge> parentHalfedges(*mesh, Halfedge());
  for (size_t i = 0; i < mesh->nHalfedges(); i++) {
    parentHalfedges[i] = origMesh.halfedge(i);
  }

  // Transfer the cut data to the new mesh.
  EdgeData<char> cut = origCut.reinterpretTo(*mesh);
  cut.setDefault(false);

  // Cut along all cut edges
  // TODO use this instead of tree trick below
  // TODO modifying while iterating
  /*
  for(Edge e : mesh->edges()) {
    if(cut[e]) {
      mesh->separateEdge(e); // all the hard works happens here
    }
  }
  */


  // TODO Right now separateEdge() can only handle cutting along a tree, so walk along tree
  std::queue<Edge> queue;
  EdgeData<char> considered(*mesh, false);

  // find any leaf edge
  Edge firstEdge;
  for (Edge e : mesh->edges()) {
    if (cut[e]) {
      int topCount = 0;
      for (Edge en : e.halfedge().vertex().adjacentEdges()) {
        if (cut[en]) topCount++;
      }
      int botCount = 0;
      for (Edge en : e.halfedge().twin().vertex().adjacentEdges()) {
        if (cut[en]) botCount++;
      }
      if (topCount == 1 || botCount == 1) {
        firstEdge = e;
        break;
      }
    }
  }
  if (firstEdge == Edge())
    throw std::runtime_error("could not find leaf edge. must cut along simple disk tree. see note");
  queue.emplace(firstEdge);
  considered[firstEdge] = true;

  // process until queue is empty
  while (!queue.empty()) {
    Edge e = queue.front();
    queue.pop();

    // Cache to restore after separate
    Halfedge oldParent = parentHalfedges[e.halfedge()];
    Halfedge oldParentT = parentHalfedges[e.halfedge().twin()];

    Halfedge newHe, newHeOpp;
    std::tie(newHe, newHeOpp) = mesh->separateEdge(e); // all the hard works happens here

    // Keep the parent map accurate
    parentHalfedges[newHe] = oldParent;
    parentHalfedges[newHeOpp] = oldParentT;

    // Add new neighbors for processing
    for (Edge en : e.halfedge().vertex().adjacentEdges()) {
      if (cut[en] && !considered[en]) {
        considered[en] = true;
        queue.emplace(en);
      }
    }
    for (Edge en : e.halfedge().twin().vertex().adjacentEdges()) {
      if (cut[en] && !considered[en]) {
        considered[en] = true;
        queue.emplace(en);
      }
    }
  }

  // Clear out any stale boundary parents
  for (Halfedge he : mesh->exteriorHalfedges()) {
    parentHalfedges[he] = Halfedge();
  }

  return std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, HalfedgeData<Halfedge>>{std::move(mesh), parentHalfedges};
}


} // namespace surface
} // namespace geometrycentral

