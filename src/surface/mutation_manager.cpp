#include "geometrycentral/surface/mutation_manager.h"

#include "geometrycentral/utilities/elementary_geometry.h"

namespace geometrycentral {
namespace surface {

// ======================================================
// ======== Construtors
// ======================================================

MutationManager::MutationManager(ManifoldSurfaceMesh& mesh_, VertexPositionGeometry& geometry_)
    : mesh(mesh_), geometry(&geometry_) {
  setConstraintCallbacks();
}

MutationManager::MutationManager(ManifoldSurfaceMesh& mesh_) : mesh(mesh_) { setConstraintCallbacks(); }

// ======================================================
// ======== Constraints
// ======================================================

void MutationManager::setConstraintCallbacks() {
  // == Register a callback to maintain edge constraint data after edge splits
  auto cacheConstraintsBeforeEdgeSplit = [&](Edge oldE, double tSplit) -> std::array<bool, 3> {
    bool wasSplittable = maySplitEdge(oldE);
    bool wasCollapsible = mayCollapseEdge(oldE);
    bool wasFlippable = mayFlipEdge(oldE);
    return {wasSplittable, wasCollapsible, wasFlippable};
  };
  auto updateConstraintsAfterEdgeSplit = [&](Halfedge newHe1, Halfedge newHe2, double tSplit,
                                             std::array<bool, 3> oldConstraints) {
    if (splittableEdges.size() > 0) {
      splittableEdges[newHe1.edge()] = oldConstraints[0];
      splittableEdges[newHe2.edge()] = oldConstraints[0];
    }
    if (collapsibleEdges.size() > 0) {
      collapsibleEdges[newHe1.edge()] = oldConstraints[1];
      collapsibleEdges[newHe2.edge()] = oldConstraints[1];
    }
    if (flippableEdges.size() > 0) {
      flippableEdges[newHe1.edge()] = oldConstraints[2];
      flippableEdges[newHe2.edge()] = oldConstraints[2];
    }
  };
  registerEdgeSplitHandlers(cacheConstraintsBeforeEdgeSplit, updateConstraintsAfterEdgeSplit);
}

void MutationManager::setRepositionableVertices(const VertexData<bool>& repositionableVertices_, bool defaultValue) {
  repositionableVertices = repositionableVertices_;
  repositionableVertices.setDefault(defaultValue);
}
void MutationManager::clearRepositionableVertices() { repositionableVertices = VertexData<bool>(); }
bool MutationManager::mayRepositionVertex(Vertex v) const {
  return repositionableVertices.size() == 0 || repositionableVertices[v];
}

void MutationManager::setSplittableEdges(const EdgeData<bool>& splittableEdges_, bool defaultValue) {
  splittableEdges = splittableEdges_;
  splittableEdges.setDefault(defaultValue);
}
void MutationManager::clearSplittableEdges() { splittableEdges = EdgeData<bool>(); }
bool MutationManager::maySplitEdge(Edge e) const { return splittableEdges.size() == 0 || splittableEdges[e]; }

void MutationManager::setCollapsibleEdges(const EdgeData<bool>& collapsibleEdges_, bool defaultValue) {
  collapsibleEdges = collapsibleEdges_;
  collapsibleEdges.setDefault(defaultValue);
}
void MutationManager::clearCollapsibleEdges() { collapsibleEdges = EdgeData<bool>(); }
bool MutationManager::mayCollapseEdge(Edge e) const {
  if (collapsibleEdges.size() > 0 && !collapsibleEdges[e]) return false;
  if (repositionableVertices.size() > 0) {
    for (Vertex v : e.adjacentVertices()) {
      if (!repositionableVertices[v]) return false;
    }
  }
  return collapsibleEdges.size() == 0 || collapsibleEdges[e];
}

void MutationManager::setFlippableEdges(const EdgeData<bool>& flippableEdges_, bool defaultValue) {
  flippableEdges = flippableEdges_;
  flippableEdges.setDefault(defaultValue);
}
void MutationManager::clearFlippableEdges() { flippableEdges = EdgeData<bool>(); }
bool MutationManager::mayFlipEdge(Edge e) const { return flippableEdges.size() == 0 || flippableEdges[e]; }

void MutationManager::setSplittableFaces(const FaceData<bool>& splittableFaces_, bool defaultValue) {
  splittableFaces = splittableFaces_;
  splittableFaces.setDefault(defaultValue);
}
void MutationManager::clearSplittableFaces() { splittableFaces = FaceData<bool>(); }
bool MutationManager::maySplitFace(Face f) const { return splittableFaces.size() == 0 || splittableFaces[f]; }

// ======================================================
// ======== Low-level mutations
// ======================================================

bool MutationManager::repositionVertex(Vertex vert, Vector3 offset) {
  if (!mayRepositionVertex(vert)) return false;

  // Invoke callbacks
  for (VertexRepositionPolicy* policy : vertexRepositionPolicies) {
    policy->beforeVertexReposition(vert, offset);
  }

  geometry->vertexPositions[vert] += offset;
  return true;
}

// Flip an edge.
bool MutationManager::flipEdge(Edge e) {
  if (!mayFlipEdge(e)) return false;

  // First do a test flip to see if its possible
  // TODO implement canFlip() to avoid this
  bool canFlip = mesh.flip(e);

  // Might not have been flippable for connectivity reasons
  if (!canFlip) {
    return false;
  }

  // Undo the test flip
  mesh.flip(e);

  // Invoke before callbacks
  for (EdgeFlipPolicy* policy : edgeFlipPolicies) {
    policy->beforeEdgeFlip(e);
  }

  // Do the actual flip
  mesh.flip(e);

  // Invoke after callbacks
  for (EdgeFlipPolicy* policy : edgeFlipPolicies) {
    policy->afterEdgeFlip(e);
  }

  return true;
}

Halfedge MutationManager::insertVertexAlongEdge(Edge e, double tSplit) {

  Vector3 newPos = Vector3::zero();
  if (geometry) {
    VertexData<Vector3>& pos = geometry->vertexPositions;
    newPos = (1. - tSplit) * pos[e.firstVertex()] + tSplit * pos[e.secondVertex()];
  }
  Halfedge he = mesh.insertVertexAlongEdge(e);
  if (geometry) {
    geometry->vertexPositions[he.vertex()] = newPos;
  }
  return he;
}

Halfedge MutationManager::splitEdge(Edge e, double tSplit) {
  if (!maySplitEdge(e)) return Halfedge();

  Vector3 newPos{0., 0., 0.};
  if (geometry) {
    VertexData<Vector3>& pos = geometry->vertexPositions;
    newPos = (1. - tSplit) * pos[e.halfedge().tailVertex()] + tSplit * pos[e.halfedge().tipVertex()];
  }
  return splitEdge(e, tSplit, newPos);
}

Halfedge MutationManager::splitEdge(Edge e, Vector3 newVertexPosition) {
  if (!maySplitEdge(e)) return Halfedge();

  double tSplit = -1;
  GC_SAFETY_ASSERT(geometry, "must have geometry to split by position");
  if (geometry) {
    // Find the nearest tCoord
    VertexData<Vector3>& pos = geometry->vertexPositions;
    Vector3 posTail = pos[e.halfedge().tailVertex()];
    Vector3 posTip = pos[e.halfedge().tipVertex()];
    tSplit = pointLineSegmentNeaestLocation(newVertexPosition, posTail, posTip);
  }

  return splitEdge(e, tSplit, newVertexPosition);
}

Halfedge MutationManager::splitEdge(Edge e, double tSplit, Vector3 newVertexPosition) {
  if (!maySplitEdge(e)) return Halfedge();

  // Invoke before callbacks
  for (EdgeSplitPolicy* policy : edgeSplitPolicies) {
    policy->beforeEdgeSplit(e, tSplit);
  }

  Halfedge newHeFront = mesh.splitEdgeTriangular(e);
  Vertex newV = newHeFront.vertex();
  if (geometry) {
    VertexData<Vector3>& pos = geometry->vertexPositions;
    pos[newV] = newVertexPosition;
  }
  Halfedge newHeBack = newHeFront.prevOrbitFace().twin().prevOrbitFace().twin();

  // Invoke after callbacks
  for (EdgeSplitPolicy* policy : edgeSplitPolicies) {
    policy->afterEdgeSplit(newHeFront, newHeBack, tSplit);
  }

  return newHeFront;
}

// Collapse an edge.
// Returns the new vertex if the edge could be collapsed, and Vertex() otherwise
Vertex MutationManager::collapseEdge(Edge e, double tCollapse) {
  if (!mayCollapseEdge(e)) return Vertex();

  Vector3 newPos{0., 0., 0.};
  if (geometry) {
    // Find the nearest tCoord
    VertexData<Vector3>& pos = geometry->vertexPositions;
    newPos = (1. - tCollapse) * pos[e.halfedge().tailVertex()] + tCollapse * pos[e.halfedge().tipVertex()];
  }
  return collapseEdge(e, tCollapse, newPos);
}

Vertex MutationManager::collapseEdge(Edge e, Vector3 newVertexPosition) {
  if (!mayCollapseEdge(e)) return Vertex();

  double tCollapse = -1;
  GC_SAFETY_ASSERT(geometry, "must have geometry to split by position");
  if (geometry) {
    // Find the nearest tCoord
    VertexData<Vector3>& pos = geometry->vertexPositions;
    Vector3 posTail = pos[e.halfedge().tailVertex()];
    Vector3 posTip = pos[e.halfedge().tipVertex()];
    double tCollapse = pointLineSegmentNeaestLocation(newVertexPosition, posTail, posTip);
  }

  return collapseEdge(e, tCollapse, newVertexPosition);
}
Vertex MutationManager::collapseEdge(Edge e, double tCollapse, Vector3 newVertexPosition) {
  if (!mayCollapseEdge(e)) return Vertex();

  // Invoke before callbacks
  // TODO need to handle possiblity that collapse fails -- check before calling
  for (EdgeCollapsePolicy* policy : edgeCollapsePolicies) {
    policy->beforeEdgeCollapse(e, tCollapse);
  }

  Vertex newV = mesh.collapseEdgeTriangular(e);
  if (newV == Vertex()) return Vertex();

  if (geometry) {
    VertexData<Vector3>& pos = geometry->vertexPositions;
    pos[newV] = newVertexPosition;
  }

  // Invoke after callbacks
  for (EdgeCollapsePolicy* policy : edgeCollapsePolicies) {
    policy->afterEdgeCollapse(newV, tCollapse);
  }

  return newV;
}

// Split a face (i.e. insert a vertex into the face)
Vertex MutationManager::splitFace(Face f, const std::vector<double>& bSplit) {
  if (!maySplitFace(f)) return Vertex();

  Vector3 newPos = Vector3::zero();
  if (geometry) {
    size_t iV = 0;
    VertexData<Vector3>& pos = geometry->vertexPositions;
    for (Vertex v : f.adjacentVertices()) {
      newPos += bSplit[iV] * pos[v];
      iV++;
    }
  }

  return splitFace(f, bSplit, newPos);
}

Vertex MutationManager::splitFace(Face f, Vector3 newVertexPosition) {
  if (!maySplitFace(f)) return Vertex();

  // TODO
  throw std::runtime_error("Face split based on vertex position not implemented yet");
  return Vertex();
}


Vertex MutationManager::splitFace(Face f, const std::vector<double>& bSplit, Vector3 newVertexPosition) {
  if (!maySplitFace(f)) return Vertex();

  // Invoke before callbacks
  for (FaceSplitPolicy* policy : faceSplitPolicies) {
    policy->beforeFaceSplit(f, bSplit);
  }

  Vertex newV = mesh.insertVertex(f);
  if (geometry) {
    VertexData<Vector3>& pos = geometry->vertexPositions;
    pos[newV] = newVertexPosition;
  }

  // Invoke after callbacks
  for (FaceSplitPolicy* policy : faceSplitPolicies) {
    policy->afterFaceSplit(newV, bSplit);
  }

  return newV;
}

Halfedge MutationManager::cutFace(Vertex vA, Vertex vB) {

  Face commonFace = sharedFace(SurfacePoint(vA), SurfacePoint(vB));
  if (commonFace == Face())
    throw std::logic_error("Can only cut mesh if the vertices at the endpoints "
                           "of the cut share a face.");

  // If segment already lies along an edge, return an existing halfedge.
  for (Halfedge he : vA.outgoingHalfedges()) {
    if (he.tipVertex() == vB) return he;
  }

  // Find out which halfedges correspond to these vertices, then call the
  // already-written connectVertices() function.
  Halfedge heA, heB;
  for (Halfedge he : commonFace.adjacentHalfedges()) {
    if (he.vertex() == vA) heA = he;
    if (he.vertex() == vB) heB = he;
  }
  return mesh.connectVertices(heA, heB);
}

Halfedge MutationManager::cutFace(SurfacePoint pA, SurfacePoint pB) {

  Face commonFace = sharedFace(pA, pB);
  if (commonFace == Face()) throw std::logic_error("Can only cut mesh if endpoints of cut share a face.");

  // TODO: seems buggy
  switch (pA.type) {
  case (SurfacePointType::Vertex): {
    switch (pB.type) {
    case (SurfacePointType::Vertex):
      return cutFace(pA.vertex, pB.vertex);
      break;
    case (SurfacePointType::Edge): {
      Halfedge heB = insertVertexAlongEdge(pB.edge, pB.tEdge);
      return cutFace(pA.vertex, heB.vertex());
      break;
    }
    case (SurfacePointType::Face):
      throw std::logic_error("Cutting into faces, i.e. from a point on the "
                             "boundary of the face to an interior point, "
                             "is not yet supported.");
      break;
    }
    break;
  }
  case (SurfacePointType::Edge): {
    switch (pB.type) {
    case (SurfacePointType::Vertex): {
      Halfedge heA = insertVertexAlongEdge(pA.edge, pA.tEdge);
      Vertex vA = heA.vertex();
      return cutFace(vA, pB.vertex);
      break;
    }
    case (SurfacePointType::Edge): {
      Halfedge heA = insertVertexAlongEdge(pA.edge, pA.tEdge);
      Halfedge heB;
      if (pA.edge == pB.edge) {
        // If points lie on same edge, then make sure we add both points in
        // correct places. ManifoldSurfaceMesh::insertVertexAlongEdge()
        // repurposes the input edge as one of the two new edges; I think the
        // one corresponding to the returned halfedge, but I'm not sure (and
        // this might change if insertVertexAlongEdge() ever changes.)
        double tB = pB.tEdge;
        double tA = pA.tEdge;
        if (tB > tA) {
          // pB gets inserted along same edge that pA is now the firstVertex of.
          heB = insertVertexAlongEdge(heA.edge(), (tB - tA) / (1. - tA));
        } else {
          // pB gets inserted along the *previous* edge along the old edge
          heB = insertVertexAlongEdge(heA.twin().next().edge(), tB / tA);
        }
      } else {
        // If pA and pB don't lie on the same edge, it's business as usual.
        Halfedge heB = insertVertexAlongEdge(pB.edge, pB.tEdge);
      }
      // instead of ManifoldSurfaceMesh::connectVertices(), which doesn't handle
      // endpoints on the same edge.
      return cutFace(heA.vertex(), heB.vertex());
      break;
    }
    case (SurfacePointType::Face): {
      throw std::logic_error("Cutting into faces, i.e. from a point on the "
                             "boundary of the face to an interior point, "
                             "is not yet supported.");
      break;
    }
    }
    break;
  }
  case (SurfacePointType::Face): {
    switch (pB.type) {
    case (SurfacePointType::Vertex): {
      throw std::logic_error("Cutting into faces, i.e. from a point on the "
                             "boundary of the face to an interior point, "
                             "is not yet supported.");
      break;
    }
    case (SurfacePointType::Edge): {
      throw std::logic_error("Cutting into faces, i.e. from a point on the "
                             "boundary of the face to an interior point, "
                             "is not yet supported.");
      break;
    }
    case (SurfacePointType::Face):
      throw std::logic_error("Cannot cut mesh along a cut within the interior "
                             "of a face; geometry-central does not "
                             "support faces with holes.");
      break;
    }
    break;
  }
  }
  return Halfedge(); // shouldn't get here
}

std::vector<Halfedge> MutationManager::cutAlongPath(const std::vector<std::vector<SurfacePoint>>& pathCurves) {

  std::vector<SurfacePoint> pathNodes;
  std::vector<std::array<size_t, 2>> pathEdges;

  for (auto curve : pathCurves) {
    size_t N = pathNodes.size();
    pathNodes.insert(pathNodes.end(), curve.begin(), curve.end());
    size_t M = curve.size();
    for (size_t i = 0; i < M - 1; i++) {
      pathEdges.push_back({N + i, N + i + 1});
    }
  }

  // Call general version
  return cutAlongPath(pathNodes, pathEdges);
}

std::vector<Halfedge> MutationManager::cutAlongPath(const std::vector<SurfacePoint>& pathNodes,
                                                    const std::vector<std::array<size_t, 2>>& pathEdges) {

  // Safety check before we mutate the mesh.
  for (auto seg : pathEdges) {
    SurfacePoint pA = pathNodes[seg[0]];
    SurfacePoint pB = pathNodes[seg[1]];

    Face commonFace = sharedFace(pA, pB);
    if (commonFace == Face()) throw std::logic_error("Each segment of cut path must lie within a single face.");
  }

  // Make sure we don't duplicate inserted vertices.
  std::vector<SurfacePoint> dedupNodes;
  std::vector<std::array<size_t, 2>> newEdges;
  for (auto seg : pathEdges) {
    SurfacePoint pA = pathNodes[seg[0]];
    SurfacePoint pB = pathNodes[seg[1]];

    if (pA == pB) continue;

    auto itA = find(dedupNodes.begin(), dedupNodes.end(), pA);
    auto itB = find(dedupNodes.begin(), dedupNodes.end(), pB);

    size_t idxA, idxB;
    std::vector<SurfacePoint> ptsToAdd; // because idk if push_back() invalidates begin() iterator
    if (itA == dedupNodes.end()) {
      ptsToAdd.push_back(pA);
      idxA = dedupNodes.size();
    } else {
      idxA = itA - dedupNodes.begin();
    }

    if (itB == dedupNodes.end()) {
      idxB = dedupNodes.size() + ptsToAdd.size();
      ptsToAdd.push_back(pB);
    } else {
      idxB = itB - dedupNodes.begin();
    }

    dedupNodes.insert(dedupNodes.end(), ptsToAdd.begin(), ptsToAdd.end());
    newEdges.push_back({idxA, idxB});
  }

  // Insert new vertices. TODO: multiple linear cuts within face not yet
  // supported.
  std::vector<Vertex> newVertices;
  for (SurfacePoint& pt : dedupNodes) {
    Halfedge he;
    Vertex newVert;
    switch (pt.type) {
    case (SurfacePointType::Vertex):
      newVert = pt.vertex;
      break;
    case (SurfacePointType::Edge):
      he = insertVertexAlongEdge(pt.edge, pt.tEdge);
      assert(he.tipVertex() == pt.edge.secondVertex());
      newVert = he.vertex();
      break;
    default:
      throw std::logic_error("Each SurfacePoint in cut path must be either "
                             "Type::Vertex or Type::Edge.");
      break;
    }
    // mesh.compress();
    newVertices.push_back(newVert);
  }

  // Then perform the face-splitting. TODO: This function doesn't handle cuts
  // which intersect other cuts.
  std::vector<Halfedge> cutHalfedges;
  for (auto& seg : newEdges) {
    Halfedge newHe = cutFace(newVertices[seg[0]], newVertices[seg[1]]);
    cutHalfedges.push_back(newHe);
  }

  return cutHalfedges;
}


// ======================================================
// ======== High-level mutations
// ======================================================


// ======================================================
// ======== Callbacks
// ======================================================

/*
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
*/

namespace {

// Helper function for registerPolicy()
template <typename T>
void pushIfSubclass(std::vector<T*>& vec, MutationPolicy* pol) {
  T* sub = dynamic_cast<T*>(pol);
  if (sub != nullptr) {
    vec.push_back(sub);
  }
}

template <typename T, typename O>
void removeUniquePtrFromVector(std::vector<std::unique_ptr<T>>& vec, O& obj) {
  for (auto it = vec.begin(); it != vec.end(); it++) {
    if (it->get() == obj) {
      vec.erase(it);
      return;
    }
  }
}

}; // namespace


MutationPolicyHandle MutationManager::registerPolicy(MutationPolicy* policyObject) {

  // Add policy to the global (memory-managed) list
  allPolicies.emplace_back(policyObject);

  // == Conditionally add policy to the various lists

  // Notice that we circumvent the type hierarchy with dynamic casts, which is generally bad practice.
  // Here, it's not too bad since it's mainly just a performance optimization---a more idiomatic approach would be to
  // have default virtual implementations for all policies in the base MutationPolicy, but this would cause lots of
  // wasted dynamic function dereferences to call no-op functions. Instead, we optimize here to pay the price once at
  // registration-time by sorting the callbacks in to categories with a dynamic_cast(). These yields an identical design
  // from the user's point of view, and (hopefully) better performance.

  pushIfSubclass(vertexRepositionPolicies, policyObject);
  pushIfSubclass(edgeFlipPolicies, policyObject);
  pushIfSubclass(edgeSplitPolicies, policyObject);
  pushIfSubclass(edgeCollapsePolicies, policyObject);
  pushIfSubclass(faceSplitPolicies, policyObject);


  // Return a handle so the user can (optionally) remove it.
  return MutationPolicyHandle(*this, policyObject);
}

void MutationManager::removePolicy(const MutationPolicyHandle& toRemove) {
  // Remove from all lists
  removeFromVector(vertexRepositionPolicies, toRemove.policy);
  removeFromVector(edgeFlipPolicies, toRemove.policy);
  removeFromVector(edgeSplitPolicies, toRemove.policy);
  removeFromVector(edgeCollapsePolicies, toRemove.policy);

  // deletion happens here
  removeUniquePtrFromVector(allPolicies, toRemove.policy);
}


MutationPolicyHandle::MutationPolicyHandle(MutationManager& manager_, MutationPolicy* policy_)
    : manager(manager_), policy(policy_) {}
void MutationPolicyHandle::remove() { manager.removePolicy(*this); }

} // namespace surface
} // namespace geometrycentral
