#include "geometrycentral/mesh_graph_algorithms.h"

#include "geometrycentral/disjoint_sets.h"

#include <algorithm>
#include <utility>
#include <vector>

namespace geometrycentral {

// Note: Assumes mesh is a single connected component
EdgeData<char> minimalSpanningTree(Geometry<Euclidean>* geometry) {

  // Preliminaries
  HalfedgeMesh* mesh = geometry->getMesh();
  GeometryCache<Euclidean>& gc = geometry->cache;
  gc.edgeLengthsQ.require();
  VertexData<size_t> vInd = mesh->getVertexIndices();

  // Store result here
  EdgeData<char> spanningTree(mesh, false);

  // Track which vertices have been connected
  DisjointSets dj(mesh->nVertices());
  size_t nConnected = 1;

  // Process the edges in order of length
  std::vector<std::pair<double, EdgePtr>> edgesByLength;
  for (EdgePtr e : mesh->edges()) {
    edgesByLength.push_back(std::make_pair(gc.edgeLengths[e], e));
  }
  std::sort(edgesByLength.begin(), edgesByLength.end());

  for (auto& edgePair : edgesByLength) {

    double len = edgePair.first;
    EdgePtr e = edgePair.second;
    VertexPtr v1 = e.halfedge().vertex();
    VertexPtr v2 = e.halfedge().twin().vertex();

    // Pass if already connected
    if (dj.find(vInd[v1]) == dj.find(vInd[v2])) {
      continue;
    }

    // Otherwise, accept the edge and union
    spanningTree[e] = true;
    dj.merge(vInd[v1], vInd[v2]);
    nConnected++;

    // Can early-out once we have connected all vertices
    if (nConnected == mesh->nVertices()) {
      break;
    }
  }

  return spanningTree;
}


// Note: Assumes mesh is a single connected component
EdgeData<char> spanningTreeBetweenVertices(Geometry<Euclidean>* geometry,
                                           const std::vector<VertexPtr>& requiredVertices) {

  // Preliminaries
  HalfedgeMesh* mesh = geometry->getMesh();
  GeometryCache<Euclidean>& gc = geometry->cache;
  gc.edgeLengthsQ.require();
  VertexData<size_t> vInd = mesh->getVertexIndices();

  // Handle special case of no required vertices
  if (requiredVertices.size() == 0) {
    return EdgeData<char>(mesh, false);
  }

  // Find a spanning tree amongst all vertices in the graph
  EdgeData<char> spanningTree = minimalSpanningTree(geometry);

  // == Trim all uneeded edges from the spanning tree
  // Process inward from leaves of tree, trimming until we hit required vertex

  // = Initialize
  // Mark needed vertices for O(1) lookup
  VertexData<char> vertexNeeded(mesh, false);
  for (VertexPtr v : requiredVertices) {
    vertexNeeded[v] = true;
  }

  // Initialize a count of vertex degrees in the spanning tree
  std::vector<VertexPtr> degree1Verts;
  VertexData<int> vDegree(mesh);
  for (VertexPtr v : mesh->vertices()) {
    int treeDegree = 0;
    for (EdgePtr e : v.adjacentEdges()) {
      if (spanningTree[e]) {
        treeDegree++;
      }
    }
    vDegree[v] = treeDegree;
    if (treeDegree == 1) {
      degree1Verts.push_back(v);
    }
  }

  // Trim the tree until there is nothing more to trim
  while (degree1Verts.size() > 0) {

    VertexPtr currV = degree1Verts.back();
    degree1Verts.pop_back();

    // Keep needed vertices
    if(vertexNeeded[currV]) {
      continue;
    }

    // Find that one edge
    HalfedgePtr treeHe;
    for (HalfedgePtr he : currV.incomingHalfedges()) {
      if (spanningTree[he.edge()]) {
        treeHe = he;
        break;
      }
    }

    // This can happen if we're removing the last edge of a tree component, such as if there is a connected component of
    // the graph with no required vertices on it
    if (treeHe == HalfedgePtr()) {
      continue;
    }

    // Un-mark the edge
    spanningTree[treeHe.edge()] = false;

    // Reduce the degree of the other vertex, and add it for processing if appropriate
    VertexPtr oppV = treeHe.vertex();
    vDegree[oppV]--;
    if (vDegree[oppV] == 1) {
      degree1Verts.push_back(oppV);
    }
  }

  return spanningTree;
}


} // namespace geometrycentral