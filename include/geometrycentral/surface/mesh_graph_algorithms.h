#pragma once

#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"

#include <unordered_map>

namespace geometrycentral {
namespace surface {

// Returns the shortest path between two vertices via Dijkstra's algorithm
// Uses sparse data structures, so complexity is output-sensitive (but still O(n logn) worst time of course)
// Returns an empty vector if the target is unreachable
std::vector<Halfedge> shortestEdgePath(IntrinsicGeometryInterface& geom, Vertex startVert, Vertex endVert);

// Return the Dijstra distance to all vertices within the ball radius
std::unordered_map<Vertex, double> vertexDijkstraDistanceWithinRadius(IntrinsicGeometryInterface& geom, Vertex startVert, double ballRad);

// Find a subset of edges which connects all vertices
// Return value holds 'true' for an edge if it is in the tree
EdgeData<char> minimalSpanningTree(IntrinsicGeometryInterface& geom);

// Returns a set of edges which connect all vertices
// Note: Uses an MST+pruning approach to find short trees in O(N logN), but not guaranteed to be minimal; that's an
// NP-hard Steiner tree problem
EdgeData<char> spanningTreeBetweenVertices(IntrinsicGeometryInterface& geom, const std::vector<Vertex>& requiredVertices);
    


} // namespace surface
} // namespace geometrycentral
