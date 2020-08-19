#include "geometrycentral/surface/mesh_graph_algorithms.h"

#include "geometrycentral/utilities/disjoint_sets.h"

#include <algorithm>
#include <queue>
#include <unordered_map>
#include <utility>
#include <vector>

namespace geometrycentral {
namespace surface {

std::vector<Halfedge> shortestEdgePath(IntrinsicGeometryInterface& geom, Vertex startVert, Vertex endVert) {

  // Early out for empty case
  if (startVert == endVert) {
    return std::vector<Halfedge>();
  }

  // Gather values
  SurfaceMesh& mesh = geom.mesh;
  geom.requireEdgeLengths();

  // Search state: incoming halfedges to each vertex, once discovered
  std::unordered_map<Vertex, Halfedge> incomingHalfedge;

  // Search state: visible neighbors eligible to expand to
  using WeightedHalfedge = std::tuple<double, Halfedge>;
  std::priority_queue<WeightedHalfedge, std::vector<WeightedHalfedge>, std::greater<WeightedHalfedge>> pq;

  // Helper to add a vertex's
  auto vertexDiscovered = [&](Vertex v) {
    return v == startVert || incomingHalfedge.find(v) != incomingHalfedge.end();
  };
  auto enqueueVertexNeighbors = [&](Vertex v, double dist) {
    for (Halfedge he : v.outgoingHalfedges()) {
      if (!vertexDiscovered(he.twin().vertex())) {
        double len = geom.edgeLengths[he.edge()];
        double targetDist = dist + len;
        pq.emplace(targetDist, he);
      }
    }
  };

  // Add initial halfedges
  enqueueVertexNeighbors(startVert, 0.);

  while (!pq.empty()) {

    // Get the next closest neighbor off the queue
    double currDist = std::get<0>(pq.top());
    Halfedge currIncomingHalfedge = std::get<1>(pq.top());
    pq.pop();

    Vertex currVert = currIncomingHalfedge.twin().vertex();
    if (vertexDiscovered(currVert)) continue;

    // Accept the neighbor
    incomingHalfedge[currVert] = currIncomingHalfedge;

    // Found path! Walk backwards to reconstruct it and return
    if (currVert == endVert) {
      std::vector<Halfedge> path;
      Vertex walkV = currVert;
      while (walkV != startVert) {
        Halfedge prevHe = incomingHalfedge[walkV];
        path.push_back(prevHe);
        walkV = prevHe.vertex();
      }

      std::reverse(std::begin(path), std::end(path));

      geom.unrequireEdgeLengths();
      return path;
    }

    // Enqueue neighbors
    enqueueVertexNeighbors(currVert, currDist);
  }

  // Didn't find path
  geom.unrequireEdgeLengths();
  return std::vector<Halfedge>();
}


std::unordered_map<Vertex, double> vertexDijkstraDistanceWithinRadius(IntrinsicGeometryInterface& geom,
                                                                      Vertex startVert, double ballRad) {

  // Gather values
  SurfaceMesh& mesh = geom.mesh;
  geom.requireEdgeLengths();

  // Search state: incoming halfedges to each vertex, once discovered
  std::unordered_map<Vertex, double> shortestDist;

  // Search state: visible neighbors eligible to expand to
  using WeightedVertex = std::tuple<double, Vertex>;
  std::priority_queue<WeightedVertex, std::vector<WeightedVertex>, std::greater<WeightedVertex>> toProcess;

  // Add initial source
  toProcess.emplace(0., startVert);

  while (!toProcess.empty()) {

    // Get the next closest neighbor off the queue
    double currDist = std::get<0>(toProcess.top());
    Vertex currVert = std::get<1>(toProcess.top());
    toProcess.pop();

    // skips stale entries
    if (shortestDist.find(currVert) != shortestDist.end()) continue;

    shortestDist[currVert] = currDist;

    for (Edge e : currVert.adjacentEdges()) {
      double len = geom.edgeLengths[e];
      Vertex targetVert = e.otherVertex(currVert);
      double targetDist = currDist + len;

      if (targetDist <= ballRad && shortestDist.find(targetVert) == shortestDist.end()) {
        toProcess.emplace(targetDist, targetVert);
      }
    }
  }

  return shortestDist;
}


/*

// Note: Assumes mesh is a single connected component
EdgeData<char> minimalSpanningTree(Geometry<Euclidean>* geometry) {

  // Preliminaries
  SurfaceMesh* mesh = geometry->getMesh();
  GeometryCache<Euclidean>& gc = geometry->cache;
  gc.requireEdgeLengths();
  VertexData<size_t> vInd = mesh->getVertexIndices();

  // Store result here
  EdgeData<char> spanningTree(*mesh, false);

  // Track which vertices have been connected
  DisjointSets dj(mesh->nVertices());
  size_t nConnected = 1;

  // Process the edges in order of length
  std::vector<std::pair<double, Edge>> edgesByLength;
  for (Edge e : mesh->edges()) {
    edgesByLength.push_back(std::make_pair(gc.edgeLengths[e], e));
  }
  std::sort(edgesByLength.begin(), edgesByLength.end());

  for (auto& edgePair : edgesByLength) {

    double len = edgePair.first;
    Edge e = edgePair.second;
    Vertex v1 = e.halfedge().vertex();
    Vertex v2 = e.halfedge().twin().vertex();

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
EdgeData<char> minimalSpanningTree(EdgeLengthGeometry* geometry) {

  // Preliminaries
  SurfaceMesh* mesh = geometry->mesh;
  geometry->requireEdgeLengths();
  VertexData<size_t> vInd = mesh->getVertexIndices();

  // Store result here
  EdgeData<char> spanningTree(*mesh, false);

  // Track which vertices have been connected
  DisjointSets dj(mesh->nVertices());
  size_t nConnected = 1;

  // Process the edges in order of length
  std::vector<std::pair<double, Edge>> edgesByLength;
  for (Edge e : mesh->edges()) {
    edgesByLength.push_back(std::make_pair(geometry->edgeLengths[e], e));
  }
  std::sort(edgesByLength.begin(), edgesByLength.end());

  for (auto& edgePair : edgesByLength) {

    double len = edgePair.first;
    Edge e = edgePair.second;
    Vertex v1 = e.halfedge().vertex();
    Vertex v2 = e.halfedge().twin().vertex();

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
EdgeData<char> spanningTreeBetweenVertices(Geometry<Euclidean>* geometry, const std::vector<Vertex>& requiredVertices)
{

  // Preliminaries
  SurfaceMesh* mesh = geometry->getMesh();
  GeometryCache<Euclidean>& gc = geometry->cache;
  gc.requireEdgeLengths();
  VertexData<size_t> vInd = mesh->getVertexIndices();

  // Handle special case of no required vertices
  if (requiredVertices.size() == 0) {
    return EdgeData<char>(*mesh, false);
  }

  // Find a spanning tree amongst all vertices in the graph
  EdgeData<char> spanningTree = minimalSpanningTree(geometry);

  // == Trim all uneeded edges from the spanning tree
  // Process inward from leaves of tree, trimming until we hit required vertex

  // = Initialize
  // Mark needed vertices for O(1) lookup
  VertexData<char> vertexNeeded(*mesh, false);
  for (Vertex v : requiredVertices) {
    vertexNeeded[v] = true;
  }

  // Initialize a count of vertex degrees in the spanning tree
  std::vector<Vertex> degree1Verts;
  VertexData<int> vDegree(*mesh);
  for (Vertex v : mesh->vertices()) {
    int treeDegree = 0;
    for (Edge e : v.adjacentEdges()) {
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

    Vertex currV = degree1Verts.back();
    degree1Verts.pop_back();

    // Keep needed vertices
    if (vertexNeeded[currV]) {
      continue;
    }

    // Find that one edge
    Halfedge treeHe;
    for (Halfedge he : currV.incomingHalfedges()) {
      if (spanningTree[he.edge()]) {
        treeHe = he;
        break;
      }
    }

    // This can happen if we're removing the last edge of a tree component, such as if there is a connected component
of
    // the graph with no required vertices on it
    if (treeHe == Halfedge()) {
      continue;
    }

    // Un-mark the edge
    spanningTree[treeHe.edge()] = false;

    // Reduce the degree of the other vertex, and add it for processing if appropriate
    Vertex oppV = treeHe.vertex();
    vDegree[oppV]--;
    if (vDegree[oppV] == 1) {
      degree1Verts.push_back(oppV);
    }
  }

  return spanningTree;
}
*/

} // namespace surface
} // namespace geometrycentral
