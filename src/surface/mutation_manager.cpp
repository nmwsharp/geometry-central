#include "geometrycentral/surface/mutation_manager.h"

#include "geometrycentral/utilities/elementary_geometry.h"

namespace geometrycentral {
namespace surface {

// ======================================================
// ======== Construtors
// ======================================================

MutationManager::MutationManager(ManifoldSurfaceMesh& mesh_, VertexPositionGeometry& geometry_)
    : mesh(mesh_), geometry(geometry_), pos(geometry.inputVertexPositions) {}

// ======================================================
// ======== Low-level mutations
// ======================================================

void MutationManager::repositionVertex(Vertex vert, Vector3 offset) {
  pos[vert] += offset;

  // Invoke any callbacks
  for (auto& fn : repositionVertexCallbackList) {
    fn(vert, offset);
  }
}

// Flip an edge.
bool MutationManager::flipEdge(Edge e) {

  bool flipped = mesh.flip(e);

  // Might not have been flippable for connectivity reasons
  if (!flipped) {
    return false;
  }

  // Invoke any edge flip callbacks
  for (auto& fn : edgeFlipCallbackList) {
    fn(e);
  }

  return true;
}

void MutationManager::splitEdge(Edge e, double tSplit) {
  Vector3 newPos = (1. - tSplit) * pos[e.halfedge().tailVertex()] + tSplit * pos[e.halfedge().tipVertex()];
  splitEdge(e, tSplit, newPos);
}

void MutationManager::splitEdge(Edge e, Vector3 newVertexPosition) {
  // Find the nearest tCoord
  Vector3 posTail = pos[e.halfedge().tailVertex()];
  Vector3 posTip = pos[e.halfedge().tipVertex()];
  double tSplit = pointLineSegmentNeaestLocation(newVertexPosition, posTail, posTip);

  splitEdge(e, tSplit, newVertexPosition);
}

void MutationManager::splitEdge(Edge e, double tSplit, Vector3 newVertexPosition) {

  Halfedge newHeFront = mesh.splitEdgeTriangular(e);
  Vertex newV = newHeFront.vertex();
  pos[newV] = newVertexPosition;
  Halfedge newHeBack = newHeFront.prevOrbitFace().twin().prevOrbitFace().twin();

  // Invoke any callbacks
  for (auto& fn : edgeSplitCallbackList) {
    fn(newHeFront, newHeBack, tSplit);
  }
}

// Collapse an edge.
// Returns true if the edge could actually be collapsed.
bool MutationManager::collapseEdge(Edge e, double tCollapse) {
  Vector3 newPos = (1. - tCollapse) * pos[e.halfedge().tailVertex()] + tCollapse * pos[e.halfedge().tipVertex()];
  return collapseEdge(e, tCollapse, newPos);
}

bool MutationManager::collapseEdge(Edge e, Vector3 newVertexPosition) {
  // Find the nearest tCoord
  Vector3 posTail = pos[e.halfedge().tailVertex()];
  Vector3 posTip = pos[e.halfedge().tipVertex()];
  double tCollapse = pointLineSegmentNeaestLocation(newVertexPosition, posTail, posTip);

  return collapseEdge(e, tCollapse, newVertexPosition);
}

bool MutationManager::collapseEdge(Edge e, double tCollapse, Vector3 newVertexPosition) {
  Vertex newV = mesh.collapseEdgeTriangular(e);
  if (newV == Vertex()) return false;
  pos[newV] = newVertexPosition;

  // Invoke any callbacks
  for (auto& fn : edgeCollapseCallbackList) {
    fn(e, newV, tCollapse);
  }

  return true;
}


// ======================================================
// ======== High-level mutations
// ======================================================


// ======================================================
// ======== Callbacks
// ======================================================

void MutationManager::managePointwiseScalarData(VertexData<double>& data) {

  { // Reposition callback
    auto updateFromGradient = [&](Vertex v, Vector3 offset) {
      // TODO estimate gradient in neighborhood then update scalar value?
    };
    // repositionVertexCallbackList.push_back(updateFromGradient);
  }

  {
      // Edge flip callback
      // (nothing needed)
  }

  { // Edge split callback
    auto updateOnEdgeSplit = [&](Halfedge newHe1, Halfedge newHe2, double tSplit) {
      Vertex newV = newHe1.vertex();
      Vertex oldVA = newHe2.tipVertex();
      Vertex oldVB = newHe1.tipVertex();
      data[newV] = (1. - tSplit) * data[oldVA] + tSplit * data[oldVB];
    };
    edgeSplitCallbackList.push_back(updateOnEdgeSplit);
  }

  { // Edge collapse callback
    auto updateOnEdgeCollapse = [&](Edge oldE, Vertex newV, double tSplit) {
      Vertex oldVA = oldE.halfedge().tailVertex();
      Vertex oldVB = oldE.halfedge().tipVertex();
      data[newV] = (1. - tSplit) * data[oldVA] + tSplit * data[oldVB];
    };
    edgeCollapseCallbackList.push_back(updateOnEdgeCollapse);
  }
}

} // namespace surface
} // namespace geometrycentral
