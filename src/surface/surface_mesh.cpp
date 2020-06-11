#include "geometrycentral/surface/surface_mesh.h"

#include "geometrycentral/utilities/combining_hash_functions.h"
#include "geometrycentral/utilities/disjoint_sets.h"
#include "geometrycentral/utilities/timing.h"

#include <algorithm>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>


using std::cout;
using std::endl;

namespace geometrycentral {
namespace surface {

// ==========================================================
// ================      Construction      ==================
// ==========================================================

// Helpers for below
namespace {

// Find an element in a sorted list
size_t halfedgeLookup(const std::vector<size_t>& compressedList, size_t target, size_t start, size_t end) {
  // Linear search is fast for small searches
  if (end - start < 20) {
    for (size_t i = start; i < end; i++) {
      if (compressedList[i] == target) {
        return i;
      }
    }
    return std::numeric_limits<size_t>::max();
  }
  // ...but we don't want to degrade to O(N^2) for really high valence vertices,
  // so fall back to a binary search
  else {
    auto loc = std::lower_bound(compressedList.begin() + start, compressedList.begin() + end, target);

    if (loc != (compressedList.begin() + end) && (target == *loc)) {
      return loc - compressedList.begin();
    } else {
      return std::numeric_limits<size_t>::max();
    }
  }
}
} // namespace


SurfaceMesh::SurfaceMesh(bool useImplicitTwin) : useImplicitTwinFlag(useImplicitTwin) {}

SurfaceMesh::SurfaceMesh(const std::vector<std::vector<size_t>>& polygons) : SurfaceMesh(polygons, {}) {}


SurfaceMesh::SurfaceMesh(const std::vector<std::vector<size_t>>& polygons,
                         const std::vector<std::vector<std::tuple<size_t, size_t>>>& twins)

    : useImplicitTwinFlag(false) {

  // Assumes that the input index set is dense. This sometimes isn't true of (eg) obj files floating around the
  // internet, so consider removing unused vertices first when reading from foreign sources.

  // START_TIMING(construction)

  // Check input list and measure some element counts
  nFacesCount = polygons.size();
  nVerticesCount = 0;
  for (const std::vector<size_t>& poly : polygons) {
    GC_SAFETY_ASSERT(poly.size() >= 3, "faces must have degree >= 3");
    for (auto i : poly) {
      nVerticesCount = std::max(nVerticesCount, i);
    }
  }
  nVerticesCount++; // 0-based means count is max+1

  // Pre-allocate face and vertex arrays
  vHalfedgeArr = std::vector<size_t>(nVerticesCount, INVALID_IND);
  fHalfedgeArr = std::vector<size_t>(nFacesCount, INVALID_IND);
  nVerticesCapacityCount = nVerticesCount;
  nVerticesFillCount = nVerticesCount;
  nFacesCapacityCount = nFacesCount;
  nFacesFillCount = nFacesCount;

  // === Walk the faces, creating halfedges. For now, don't hook up any twin or edge pointers.
  for (size_t iFace = 0; iFace < nFacesCount; iFace++) {
    const std::vector<size_t>& poly = polygons[iFace];

    // Walk around this face
    size_t faceDegree = poly.size();
    size_t prevHeInd = INVALID_IND;
    size_t firstHeInd = INVALID_IND;
    for (size_t iFaceHe = 0; iFaceHe < faceDegree; iFaceHe++) {

      size_t indTail = poly[iFaceHe];
      size_t indTip = poly[(iFaceHe + 1) % faceDegree];

      // Get an index for this halfedge
      std::tuple<size_t, size_t> heKey{indTail, indTip};
      std::tuple<size_t, size_t> heTwinKey{indTip, indTail};

      // Create a halfedge
      size_t halfedgeInd = getNewHalfedge(true).getIndex();

      // Fill arrays with nknown values and placeholders
      heNextArr[halfedgeInd] = INVALID_IND;
      heVertexArr[halfedgeInd] = indTail;
      heFaceArr[halfedgeInd] = iFace;

      // Hook up a bunch of pointers
      vHalfedgeArr[indTail] = halfedgeInd;
      if (iFaceHe == 0) {
        fHalfedgeArr[iFace] = halfedgeInd;
        firstHeInd = halfedgeInd;
      } else {
        heNextArr[prevHeInd] = halfedgeInd;
      }
      prevHeInd = halfedgeInd;
    }
    heNextArr[prevHeInd] = firstHeInd; // hook up the first next() pointer, which we missed in the loop above
  }

  // === Create edges and hook up twins
  std::vector<size_t> boundaryHalfedges;
  if (twins.empty()) {
    // Any halfedges between a pair of vertices are considered to be incident on the same edge

    // If we've already seen this vertex pair at least once, this holds index of the most recent halfedge encountered
    // incident on that edge.
    std::unordered_map<std::tuple<size_t, size_t>, size_t> edgeHistory;
    size_t iHe = 0;
    for (size_t iFace = 0; iFace < nFacesCount; iFace++) {
      const std::vector<size_t>& poly = polygons[iFace];

      // Walk around this face
      size_t faceDegree = poly.size();
      for (size_t iFaceHe = 0; iFaceHe < faceDegree; iFaceHe++) {

        size_t indTail = poly[iFaceHe];
        size_t indTip = poly[(iFaceHe + 1) % faceDegree];

        // Get a key for this edge
        std::tuple<size_t, size_t> eKey{std::min(indTail, indTip), std::max(indTail, indTip)};

        if (edgeHistory.find(eKey) == edgeHistory.end()) {
          // This is the first time we've ever seen this edge, create a new edge object
          size_t newEdgeInd = getNewEdge().getIndex();
          heEdgeArr[iHe] = newEdgeInd;
          heSiblingArr[iHe] = INVALID_IND;
          eHalfedgeArr[newEdgeInd] = iHe;
          edgeHistory[eKey] = iHe;
        } else {
          // We're already seen this edge, connect to the previous halfedge incident on the edge
          size_t iPrevHe = edgeHistory[eKey];
          heSiblingArr[iHe] = iPrevHe;
          heEdgeArr[iHe] = heEdgeArr[iPrevHe];
          edgeHistory[eKey] = iHe;
        }
        iHe++;
      }
    }

    // Complete the sibling cycle by follwing backwards each edge until we reach the first sibling-less entry
    for (auto& entry : edgeHistory) {
      size_t lastHe = entry.second;

      if (heSiblingArr[lastHe] == INVALID_IND) {
        // Any edges which never got any sibling entries at all are boundary halfedges
        boundaryHalfedges.push_back(lastHe);
        continue;
      }

      // Get the index of the first halfedge in the sibling cycle to complete the cycle
      size_t currHe = lastHe;
      while (heSiblingArr[currHe] != INVALID_IND) {
        currHe = heSiblingArr[currHe];
      }
      heSiblingArr[currHe] = lastHe; // connect the first to the last
    }

  } else {
    // DisjointSets djSet
  }


  // == Process boundaries
  // Note: nonmanfiold meshes do not necessarily have coherent boundary loops, so there's not necessarily too much we
  // can do here.

  // Create exterior halfedges
  for (size_t bHe : boundaryHalfedges) {
    size_t tHe = getNewHalfedge(false).getIndex();
    heSiblingArr[bHe] = tHe;
    heSiblingArr[tHe] = bHe;
    heNextArr[tHe] = INVALID_IND;
    heVertexArr[tHe] = heVertexArr[heNextArr[bHe]];
    heFaceArr[tHe] = INVALID_IND;
    heEdgeArr[tHe] = heEdgeArr[bHe];

    // ensure that vertex.halfedge() always starts a half-disk for boundary vertices
    vHalfedgeArr[heVertexArr[bHe]] = bHe;
  }

  // Some vertices will have nice boundary structures, like a half-disk, or collection of half-disks, where there are
  // reasonable choices of he.next().
  //
  // However, there can be other nonmanifold mesh structures which simply do not admit any good choice of he.next().
  // (imagine a sundial structure, where we take a manifold vertex and tack on a single triangle).
  //
  // As such, our basic strategy is to form maximal chains of reasonable next() references wherever vertices have
  // incoming and outgoing boundary halfedges that can be matched up. Then, we will close any unclosed chains by
  // connecting the first to the last entry. Finally, each chain will correspond to a boundary loop.

  // We will need to loop around vertice below (on a possibly nonmanifold mesh)
  populateVertexIterationCache(false);

  // Hook up next() references point to a next exterior halfedge at the same vertex
  std::vector<size_t> incomingExteriorHalfedges;
  std::vector<size_t> outgoingExteriorHalfedges;
  for (Vertex v : vertices()) {
    std::cout << "Testing in out at " << v << std::endl;
    incomingExteriorHalfedges.clear();
    outgoingExteriorHalfedges.clear();

    // Manually traverse the halfegdes incident on this vertex
    size_t rangeStart = vertexIterationCacheVertexStart[v.getIndex()];
    size_t rangeEnd = vertexIterationCacheVertexStart[v.getIndex() + 1];
    for (size_t traverseInd = rangeStart; traverseInd < rangeEnd; traverseInd++) {
      size_t iHe = vertexIterationCacheHeIndex[traverseInd];

      std::cout << "  he out " << iHe << std::endl;
      if (heFaceArr[iHe] == INVALID_IND) { // test if its exterior
        std::cout << "    exterior outgoing! " << iHe << std::endl;
        outgoingExteriorHalfedges.push_back(iHe);
      }
      if (heFaceArr[heTwin(iHe)] == INVALID_IND) {
        std::cout << "    twin interior outgoing! " << heTwin(iHe) << std::endl;
        incomingExteriorHalfedges.push_back(heTwin(iHe));
      }

      // Match up as many next()s as we can (some will still be unpaired)
      while (!outgoingExteriorHalfedges.empty() && !incomingExteriorHalfedges.empty()) {
        size_t iHeOut = outgoingExteriorHalfedges.back();
        size_t iHeIn = incomingExteriorHalfedges.back();
        outgoingExteriorHalfedges.pop_back();
        incomingExteriorHalfedges.pop_back();

        std::cout << "    at " << v << " matching next of incoming " << iHeIn << " as " << iHeOut << std::endl;
        heNextArr[iHeIn] = iHeOut;
      }
    }
  }

  // Check which halfedges are the next() of some halfedge
  std::vector<char> hasIncomingNext(nHalfedgesCount, false);
  for (size_t iHe = 0; iHe < nHalfedgesCount; iHe++) {
    if (heNextArr[iHe] != INVALID_IND) {
      if (hasIncomingNext[heNextArr[iHe]]) {
        // TODO FIXME sanity check
        throw std::runtime_error("Halfedge has multiple incoming next");
      }
      hasIncomingNext[heNextArr[iHe]] = true;
    }
  }

  // Close chains of next()
  for (size_t iHe = 0; iHe < nHalfedgesCount; iHe++) {
    if (hasIncomingNext[iHe]) continue;
    // should only proceed for exterior halfedges
    // std::cout << "halfedge " << iHe << " has no incoming next" << std::endl;

    size_t firstHe = iHe;
    size_t currHe = firstHe;
    while (heNextArr[currHe] != INVALID_IND) {
      currHe = heNextArr[currHe];
    }

    // connect last to first
    std::cout << "connecting heNext[" << currHe << "] = " << firstHe << std::endl;
    heNextArr[currHe] = firstHe;
  }


  // Create a boundary loop for each cycle of next()
  for (size_t iHe = 0; iHe < nHalfedgesCount; iHe++) {
    if (heFaceArr[iHe] != INVALID_IND) continue;

    // Create a new boundary loop
    size_t iBl = getNewBoundaryLoop().getIndex();
    fHalfedgeArr[iBl] = iHe;
    size_t currHe = iHe;
    do {
      heFaceArr[currHe] = iBl;
      currHe = heNextArr[currHe];
    } while (currHe != iHe);
  }


  isCompressedFlag = true;
  // TODO compress here?
  validateConnectivity(); // TODO FIXME
}


SurfaceMesh::~SurfaceMesh() {
  for (auto& f : meshDeleteCallbackList) {
    f();
  }
}


// ==========================================================
// ================       Utilities        ==================
// ==========================================================

void SurfaceMesh::printStatistics() const {
  std::cout << "Halfedge mesh with: " << std::endl;
  std::cout << "    # verts =  " << nVertices() << std::endl;
  std::cout << "    # edges =  " << nEdges() << std::endl;
  std::cout << "    # faces =  " << nFaces() << std::endl;
  std::cout << "    # halfedges =  " << nHalfedges() << "  (" << nInteriorHalfedges() << " interior, "
            << nExteriorHalfedges() << " exterior)" << std::endl;
  std::cout << "      and " << nBoundaryLoops() << " boundary components. " << std::endl;
}

bool SurfaceMesh::isTriangular() {
  for (Face f : faces()) {
    if (!f.isTriangle()) {
      return false;
    }
  }
  return true;
}


bool SurfaceMesh::isManifold() {
  for (Edge e : edges()) {
    if (!e.isManifold()) return false;
  }
  for (Vertex v : vertices()) {
    if (!v.isManifold()) return false;
  }
  return true;
}

size_t SurfaceMesh::nConnectedComponents() {
  VertexData<size_t> vertInd = getVertexIndices();
  DisjointSets dj(nVertices());
  for (Edge e : edges()) {
    dj.merge(vertInd[e.halfedge().vertex()], vertInd[e.halfedge().twin().vertex()]);
  }
  std::unordered_set<size_t> distinctComponents;
  for (size_t i = 0; i < nVertices(); i++) {
    distinctComponents.insert(dj.find(i));
  }
  return distinctComponents.size();
}

size_t SurfaceMesh::nInteriorVertices() {
  size_t nInteriorVertices = 0;
  for (const Vertex v : vertices()) {
    if (!v.isBoundary()) {
      nInteriorVertices++;
    }
  }
  return nInteriorVertices;
}


VertexData<size_t> SurfaceMesh::getVertexIndices() {
  VertexData<size_t> indices(*this);
  size_t i = 0;
  for (Vertex v : vertices()) {
    indices[v] = i;
    i++;
  }
  return indices;
}

VertexData<size_t> SurfaceMesh::getInteriorVertexIndices() {
  VertexData<size_t> indices(*this);
  size_t i = 0;
  for (Vertex v : vertices()) {
    if (v.isBoundary()) {
      indices[v] = INVALID_IND;
    } else {
      indices[v] = i;
      i++;
    }
  }
  return indices;
}

HalfedgeData<size_t> SurfaceMesh::getHalfedgeIndices() {
  HalfedgeData<size_t> indices(*this);
  size_t i = 0;
  for (Halfedge he : halfedges()) {
    indices[he] = i;
    i++;
  }
  return indices;
}

CornerData<size_t> SurfaceMesh::getCornerIndices() {
  CornerData<size_t> indices(*this);
  size_t i = 0;
  for (Corner c : corners()) {
    indices[c] = i;
    i++;
  }
  return indices;
}

EdgeData<size_t> SurfaceMesh::getEdgeIndices() {
  EdgeData<size_t> indices(*this);
  size_t i = 0;
  for (Edge e : edges()) {
    indices[e] = i;
    i++;
  }
  return indices;
}


FaceData<size_t> SurfaceMesh::getFaceIndices() {
  FaceData<size_t> indices(*this);
  size_t i = 0;
  for (Face f : faces()) {
    indices[f] = i;
    i++;
  }
  return indices;
}

BoundaryLoopData<size_t> SurfaceMesh::getBoundaryLoopIndices() {
  BoundaryLoopData<size_t> indices(*this);
  size_t i = 0;
  for (BoundaryLoop bl : boundaryLoops()) {
    indices[bl] = i;
    i++;
  }
  return indices;
}


std::unique_ptr<SurfaceMesh> SurfaceMesh::copy() const {
  return copyToSurfaceMesh();
}

std::unique_ptr<SurfaceMesh> SurfaceMesh::copyToSurfaceMesh() const {
  SurfaceMesh* newMesh = new SurfaceMesh(false);
  copyInternalFields(*newMesh);
  return std::unique_ptr<SurfaceMesh>(newMesh);
}

void SurfaceMesh::copyInternalFields(SurfaceMesh& target) const {
  // == Copy _all_ the fields!

  // Raw data buffers (underlying std::vectors duplicate storage automatically)
  target.heNextArr = heNextArr;
  target.heVertexArr = heVertexArr;
  target.heFaceArr = heFaceArr;
  target.vHalfedgeArr = vHalfedgeArr;
  target.fHalfedgeArr = fHalfedgeArr;
  target.heSiblingArr = heSiblingArr;
  target.heEdgeArr = heEdgeArr;
  target.eHalfedgeArr = eHalfedgeArr;

  // counts and flags
  target.nHalfedgesCount = nHalfedgesCount;
  target.nInteriorHalfedgesCount = nInteriorHalfedgesCount;
  target.nEdgesCount = nEdgesCount;
  target.nVerticesCount = nVerticesCount;
  target.nFacesCount = nFacesCount;
  target.nBoundaryLoopsCount = nBoundaryLoopsCount;
  target.nVerticesCapacityCount = nVerticesCapacityCount;
  target.nHalfedgesCapacityCount = nHalfedgesCapacityCount;
  target.nEdgesCapacityCount = nEdgesCapacityCount;
  target.nFacesCapacityCount = nFacesCapacityCount;
  target.nVerticesFillCount = nVerticesFillCount;
  target.nHalfedgesFillCount = nHalfedgesFillCount;
  target.nEdgesFillCount = nEdgesFillCount;
  target.nFacesFillCount = nFacesFillCount;
  target.nBoundaryLoopsFillCount = nBoundaryLoopsFillCount;

  target.isCompressedFlag = isCompressedFlag;

  // Note: _don't_ copy callbacks lists! New mesh has new callbacks
}

std::vector<std::vector<size_t>> SurfaceMesh::getFaceVertexList() {

  std::vector<std::vector<size_t>> result;

  VertexData<size_t> vInd = getVertexIndices();
  for (Face f : faces()) {
    std::vector<size_t> faceList;
    for (Vertex v : f.adjacentVertices()) {
      faceList.push_back(vInd[v]);
    }
    result.push_back(faceList);
  }

  return result;
}


void SurfaceMesh::populateVertexIterationCache(bool skipDead) {
  // TODO will have problems with constness here, need mutable

  std::cout << "========= POPULATING VERTEX ITERATION CACHE ==============" << std::endl;

  // First, count the degree of every vertex
  std::vector<size_t> vDegree(nVerticesFillCount, 0);
  for (size_t iHe = 0; iHe < nHalfedgesFillCount; iHe++) {
    if (skipDead && halfedgeIsDead(iHe)) continue;
    vDegree[heVertex(iHe)]++;
  }

  // Build a sum-array of the number of vertices up to that point
  vertexIterationCacheVertexStart.resize(nVerticesFillCount + 1);
  size_t heSum = 0;
  for (size_t iV = 0; iV < nVerticesFillCount; iV++) {
    vertexIterationCacheVertexStart[iV] = heSum;
    heSum += vDegree[iV];
  }
  // add one extra element, so we can check end bound by indexing+1
  vertexIterationCacheVertexStart[nVerticesFillCount] = heSum;

  // Build a compressed array of the halfedges at each vertex
  std::vector<size_t> currVertexCacheEntry = vertexIterationCacheVertexStart;
  vertexIterationCacheHeIndex.resize(nHalfedgesFillCount);
  for (size_t iHe = 0; iHe < nHalfedgesFillCount; iHe++) {
    if (skipDead && halfedgeIsDead(iHe)) continue;
    size_t iV = heVertex(iHe);
    size_t entryInd = currVertexCacheEntry[iV];
    vertexIterationCacheHeIndex[entryInd] = iHe;
    currVertexCacheEntry[iV]++;
  }

  vertexIterationCacheTick = modificationTick;
}

void SurfaceMesh::validateConnectivity() {

  // Sanity check sizes and counts
  if (nInteriorHalfedges() + nExteriorHalfedges() != nHalfedges())
    throw std::logic_error("nInterior + nImaginary != nTotal halfedges");
  if (nHalfedgesCount > nHalfedgesFillCount) throw std::logic_error("halfedge count > halfedge fill");
  if (nHalfedgesFillCount > nHalfedgesCapacityCount) throw std::logic_error("halfedge fill > halfedge capacity");

  if (nVerticesCount > nVerticesFillCount) throw std::logic_error("vertex count > vertex fill");
  if (nVerticesFillCount > nVerticesCapacityCount) throw std::logic_error("vertex fill > vertex capacity");

  if (nFacesCount > nFacesFillCount) throw std::logic_error("face count > face fill");
  if (nFacesFillCount + nBoundaryLoopsCount > nFacesCapacityCount)
    throw std::logic_error("face + bl fill > face capacity");

  // Helpers to check the validity of references
  auto validateVertex = [&](size_t iV, std::string msg) {
    if (iV >= nVerticesFillCount || vertexIsDead(iV)) throw std::logic_error(msg + " - bad vertex reference");
  };
  auto validateHalfedge = [&](size_t iHe, std::string msg) {
    if (iHe >= nHalfedgesFillCount || halfedgeIsDead(iHe)) throw std::logic_error(msg + " - bad halfedge reference");
  };
  auto validateEdge = [&](size_t iE, std::string msg) {
    if (iE >= nEdgesFillCount || edgeIsDead(iE)) throw std::logic_error(msg + " - bad edge reference");
  };
  auto validateFace = [&](size_t iF, std::string msg) {
    if (iF >= nFacesCapacityCount || faceIsDead(iF) ||
        (iF >= nFacesFillCount && iF < nFacesCapacityCount - nBoundaryLoopsCount)) {
      // third case checks the dead zone between faces and boundary loop indices
      throw std::logic_error(msg + " - bad face reference");
    }
  };

  // == Halfedges

  // Check valid pointers
  // Note: we intentionally avoid using iterators here, because they can be hard to debug when things are broken.
  for (size_t iHe = 0; iHe < nHalfedgesFillCount; iHe++) {
    if (halfedgeIsDead(iHe)) continue;
    validateHalfedge(heTwin(iHe), "he.twin()");
    validateHalfedge(heNextArr[iHe], "he.next()");
    validateVertex(heVertexArr[iHe], "he.vertex()");
    validateEdge(heEdge(iHe), "he.edge()");
    validateFace(heFaceArr[iHe], "he.face()");
  }
  for (size_t iV = 0; iV < nVerticesFillCount; iV++) {
    if (vertexIsDead(iV)) continue;
    validateHalfedge(vHalfedgeArr[iV], "v.halfedge()");
  }
  for (size_t iE = 0; iE < nEdgesFillCount; iE++) {
    if (edgeIsDead(iE)) continue;
    validateHalfedge(eHalfedge(iE), "e.halfedge()");
  }
  for (size_t iF = 0; iF < nFacesFillCount; iF++) {
    if (faceIsDead(iF)) continue;
    validateHalfedge(fHalfedgeArr[iF], "f.halfedge()");
  }

  // Check edge and twin sanity
  for (Halfedge he : halfedges()) {
    if (he != he.twin().twin()) throw std::logic_error("twins not reflective");
    if (he == he.twin()) throw std::logic_error("self-twin");
  }
  for (Edge e : edges()) {
    for (Halfedge he : e.adjacentHalfedges()) {
      if (e != he.edge()) throw std::logic_error("edge.halfedge doesn't match halfedge.edge");
    }
  }

  // Check face & next sanity
  for (Face f : faces()) {
    if (f.halfedge().face() != f) throw std::logic_error("f.halfedge().face() is not f");

    Halfedge currHe = f.halfedge();
    Halfedge firstHe = f.halfedge();
    size_t count = 0;
    do {
      if (currHe.face() != f) throw std::logic_error("face.halfedge doesn't match halfedge.face");
      currHe = currHe.next();
      count++;
      if (count > nHalfedgesCount) throw std::logic_error("next forms non-face loop");
    } while (currHe != firstHe);

    if (count < 2) throw std::logic_error("face of degree < 2");
  }


  // Check face & next sanity
  for (BoundaryLoop b : boundaryLoops()) {
    if (b.asFace().asBoundaryLoop() != b) throw std::logic_error("b.asFace().asBoundaryLoop() is not fixed point");

    if (b.halfedge().face() != b.asFace()) throw std::logic_error("bl.halfedge().face() is not bl");

    if (!b.halfedge().face().isBoundaryLoop()) throw std::logic_error("bl.halfedge().face() is not a boundary loop");

    Halfedge currHe = b.halfedge();
    Halfedge firstHe = b.halfedge();
    size_t count = 0;
    do {
      if (!currHe.face().isBoundaryLoop())
        throw std::logic_error("walking around boundary loop yielded he.face() which is not a boundary loop");
      if (currHe.face().asBoundaryLoop() != b)
        throw std::logic_error("(boundary loop) face.halfedge doesn't match halfedge.face");
      currHe = currHe.next();
      count++;
      if (count > nHalfedgesCount) throw std::logic_error("(boundary loop) next forms non-face loop");
    } while (currHe != firstHe);

    if (count < 2) throw std::logic_error("(boundary loop) face of degree < 2");
  }


  // Check halfedge neighborhood sanity
  std::vector<char> halfedgeSeen(nHalfedgesCapacityCount, false);
  for (Halfedge he : halfedges()) {

    // Check that he.twin().next() locally orbis a vertex
    if (he.vertex() != he.twin().next().vertex()) throw std::logic_error("halfedge vertices don't match");

    // Check that boundary rules are observed
    if (!he.isInterior()) {
      if (he == he.edge().halfedge()) throw std::logic_error("exterior halfedge is e.halfedge()");
      if (!he.twin().isInterior()) throw std::logic_error("he and he.twin() are both exterior");
    }

    // This can happen in irregular triangulations
    // if (he.vertex == he.next->twin->vertex) throw std::logic_error("halfedge face spur");

    // Check halfedge orbit sanity (useful if halfedge doesn't appear in face)
    Halfedge currHe = he;
    if (halfedgeSeen[currHe.getIndex()]) continue;
    Halfedge firstHe = he;
    size_t count = 0;
    do {
      if (currHe.face() != he.face()) throw std::logic_error("he.next.**.face doesn't match he.face");
      halfedgeSeen[currHe.getIndex()] = true;
      currHe = currHe.next();
      count++;
      if (count > nHalfedgesCount) throw std::logic_error("next forms non-face loop");
    } while (currHe != firstHe);
  }

  // Check vertex orbit sanity
  for (Vertex v : vertices()) {
    Halfedge currHe = v.halfedge();
    Halfedge firstHe = v.halfedge();
    size_t count = 0;
    do {
      if (currHe.vertex() != v) throw std::logic_error("vertex.halfedge doesn't match halfedge.vertex");
      currHe = currHe.twin().next();
      count++;
      if (count > nHalfedgesCount) throw std::logic_error("twin->next forms non-vertex loop");
    } while (currHe != firstHe);
  }

  // Verify boundary rules are correct
  for (Vertex v : vertices()) {

    // Manually check if this is a boundary vertex
    size_t boundaryHeCount = 0;
    for (Halfedge he : v.outgoingHalfedges()) {
      if (!he.isInterior()) boundaryHeCount++;
    }

    if (boundaryHeCount > 1) {
      throw std::logic_error("multiple boundaries incident on vertex");
    }
    bool hasBoundaryHe = boundaryHeCount > 0;

    if (hasBoundaryHe) {
      if (!v.halfedge().isInterior()) {
        throw std::logic_error("v.halfedge() is exterior");
      }
      if (v.halfedge().twin().isInterior()) {
        throw std::logic_error("v.halfedge() does not border boundary on a boundary vertex");
      }
      if (!v.isBoundary()) {
        throw std::logic_error("computed v.isBoundary is wrong");
      }
    }
  }

  // Check manifoldness, if the mesh should be manifold
  if (usesImplictTwin()) {
    // Right now (and for the forseeable future) the usesImplictTwin is really equivalent to being of subclass
    // ManifoldSurfaceMesh, so this error prints that, which is easier to understand to the user.
    if (!isManifold()) {
      throw std::logic_error("Mesh with underlying type ManifoldSurfaceMesh is not actually manifold, likely due to a "
                             "nonmanifold vertex with disconnected neighborhoods.");
    }
  }
}


Vertex SurfaceMesh::getNewVertex() {

  // The boring case, when no resize is needed
  if (nVerticesFillCount < nVerticesCapacityCount) {
    // No work needed
  }
  // The intesting case, where vectors resize
  else {
    size_t newCapacity = nVerticesCapacityCount * 2;

    // Resize internal arrays
    vHalfedgeArr.resize(newCapacity);

    nVerticesCapacityCount = newCapacity;

    // Invoke relevant callback functions
    for (auto& f : vertexExpandCallbackList) {
      f(newCapacity);
    }
  }

  nVerticesFillCount++;
  nVerticesCount++;

  modificationTick++;
  return Vertex(this, nVerticesFillCount - 1);
}

Halfedge SurfaceMesh::getNewHalfedge(bool isInterior) {

  if (usesImplictTwin()) {
    throw std::logic_error("cannot construct a single new halfedge with implicit twin convention");
  }

  // The boring case, when no resize is needed
  if (nHalfedgesFillCount < nHalfedgesCapacityCount) {
    // No work needed
  }
  // The intesting case, where vectors resize
  else {
    size_t newHalfedgeCapacity = std::max(nHalfedgesCapacityCount * 2, (size_t)1);

    // Resize internal arrays
    heNextArr.resize(newHalfedgeCapacity);
    heVertexArr.resize(newHalfedgeCapacity);
    heFaceArr.resize(newHalfedgeCapacity);
    if (!usesImplictTwin()) { // must enter this case, see test above
      heSiblingArr.resize(newHalfedgeCapacity);
      heEdgeArr.resize(newHalfedgeCapacity);
    }

    nHalfedgesCapacityCount = newHalfedgeCapacity;

    // Invoke relevant callback functions
    for (auto& f : halfedgeExpandCallbackList) {
      f(newHalfedgeCapacity);
    }
  }

  nHalfedgesFillCount++;
  nHalfedgesCount++;
  if (isInterior) {
    nInteriorHalfedgesCount++;
  }

  modificationTick++;
  return Halfedge(this, nHalfedgesFillCount - 1);
}

Edge SurfaceMesh::getNewEdge() {

  if (usesImplictTwin()) {
    throw std::logic_error("cannot construct a single new edge with implicit twin convention");
  }

  // The boring case, when no resize is needed
  if (nEdgesFillCount < nEdgesCapacityCount) {
    // No work needed
  }
  // The intesting case, where vectors resize
  else {
    size_t newEdgeCapacity = std::max(nEdgesCapacityCount * 2, (size_t)1);

    nEdgesCapacityCount = newEdgeCapacity;

    if (!usesImplictTwin()) { // must enter this case, see test above
      eHalfedgeArr.resize(newEdgeCapacity);
    }

    // Invoke relevant callback functions
    for (auto& f : edgeExpandCallbackList) {
      f(newEdgeCapacity);
    }
  }

  nEdgesFillCount++;
  nEdgesCount++;

  modificationTick++;
  return Edge(this, nEdgesFillCount - 1);
}

Halfedge SurfaceMesh::getNewEdgeTriple(bool onBoundary) {

  // == Get two halfedges and one edge
  // recall that these capacities should always be in sync, so we resize and expand for either both edges and
  // halfedges, or neither

  if (nHalfedgesFillCount + 1 < nHalfedgesCapacityCount) {
    GC_SAFETY_ASSERT(nEdgesFillCount < nEdgesCapacityCount, "edge capacity is out of sync with halfedge capacity");

    // No work needed
  } else {

    size_t initHalfedgeCapacity = nHalfedgesCapacityCount; // keep track before we start modifying for clarity
    size_t initEdgeCapacity = nEdgesCapacityCount;
    size_t newHalfedgeCapacity = std::max(initHalfedgeCapacity * 2, (size_t)2); // double the capacity
    size_t newEdgeCapacity = std::max(initEdgeCapacity * 2, (size_t)1);

    { // expand halfedge list

      // Resize internal arrays
      heNextArr.resize(newHalfedgeCapacity);
      heVertexArr.resize(newHalfedgeCapacity);
      heFaceArr.resize(newHalfedgeCapacity);
      if (!usesImplictTwin()) {
        heSiblingArr.resize(newHalfedgeCapacity);
        heEdgeArr.resize(newHalfedgeCapacity);
      }

      nHalfedgesCapacityCount = newHalfedgeCapacity;

      // Invoke relevant callback functions
      for (auto& f : halfedgeExpandCallbackList) {
        f(newHalfedgeCapacity);
      }
    }

    { // expand edges
      nEdgesCapacityCount = newEdgeCapacity;

      if (!usesImplictTwin()) {
        eHalfedgeArr.resize(newEdgeCapacity);
      }

      // Invoke relevant callback functions
      for (auto& f : edgeExpandCallbackList) {
        f(newEdgeCapacity);
      }
    }
  }


  // == Get one

  // Fill connectivity buffers if needed
  if (!usesImplictTwin()) {
    heSiblingArr[nHalfedgesFillCount] = nHalfedgesFillCount + 1;
    heSiblingArr[nHalfedgesFillCount + 1] = nHalfedgesFillCount;
    heEdgeArr[nHalfedgesFillCount] = nEdgesFillCount;
    heEdgeArr[nHalfedgesFillCount + 1] = nEdgesFillCount;
    eHalfedgeArr[nEdgesFillCount] = nHalfedgesFillCount;
  }

  nHalfedgesFillCount += 2;
  nHalfedgesCount += 2;
  if (onBoundary) {
    nInteriorHalfedgesCount += 1;
  } else {
    nInteriorHalfedgesCount += 2;
  }
  nEdgesFillCount++;
  nEdgesCount++;

  modificationTick++;
  return Halfedge(this, nHalfedgesFillCount - 2);
}


Face SurfaceMesh::getNewFace() {

  // The boring case, when no resize is needed
  if (nFacesFillCount + nBoundaryLoopsCount < nFacesCapacityCount) {
    // No work needed
  }
  // The intesting case, where vectors resize
  else {
    expandFaceStorage();
  }

  nFacesCount++;
  nFacesFillCount++;

  modificationTick++;
  return Face(this, nFacesFillCount - 1);
}

BoundaryLoop SurfaceMesh::getNewBoundaryLoop() {

  // The boring case, when no resize is needed
  if (nFacesFillCount + nBoundaryLoopsCount < nFacesCapacityCount) {
    // No work needed
  }
  // The intesting case, where vectors resize
  else {
    expandFaceStorage();
  }

  nBoundaryLoopsCount++;
  nBoundaryLoopsFillCount++;

  modificationTick++;
  return BoundaryLoop(this, nFacesCapacityCount - nBoundaryLoopsFillCount);
}

void SurfaceMesh::expandFaceStorage() {
  size_t newCapacity = nFacesCapacityCount * 2;

  // Resize internal arrays
  fHalfedgeArr.resize(newCapacity);

  // Scooch boundary data back
  for (size_t iBack = 0; iBack < nBoundaryLoopsFillCount; iBack++) {
    size_t iOld = nFacesCapacityCount - iBack - 1;
    size_t iNew = fHalfedgeArr.size() - iBack - 1;
    fHalfedgeArr[iNew] = fHalfedgeArr[iOld];
    fHalfedgeArr[iOld] = INVALID_IND; // will help catch bugs
  }

  // Scooch back he.face() indices that point to boundary loops
  for (size_t iHe = 0; iHe < nHalfedgesFillCount; iHe++) {
    if (halfedgeIsDead(iHe)) {
      continue;
    }
    if (heFaceArr[iHe] >= nFacesFillCount) {
      heFaceArr[iHe] += (newCapacity - nFacesCapacityCount);
    }
  }

  nFacesCapacityCount = newCapacity;

  // Invoke relevant callback functions
  for (auto& f : faceExpandCallbackList) {
    f(newCapacity);
  }

  modificationTick++;
}


void SurfaceMesh::deleteEdgeTriple(Halfedge he) {
  // Be sure we have the canonical halfedge
  he = he.edge().halfedge();
  bool isBoundary = he.twin().isInterior();
  size_t iHe = he.getIndex();
  size_t iHeT = he.twin().getIndex();

  heNextArr[iHe] = INVALID_IND;
  heVertexArr[iHe] = INVALID_IND;
  heFaceArr[iHe] = INVALID_IND;

  heNextArr[iHeT] = INVALID_IND;
  heVertexArr[iHeT] = INVALID_IND;
  heFaceArr[iHeT] = INVALID_IND;

  nHalfedgesCount -= 2;
  if (isBoundary) {
    nInteriorHalfedgesCount--;
  } else {
    nInteriorHalfedgesCount -= 2;
  }

  modificationTick++;
  isCompressedFlag = false;
}

void SurfaceMesh::deleteElement(Vertex v) {
  size_t iV = v.getIndex();

  vHalfedgeArr[iV] = INVALID_IND;
  nVerticesCount--;

  modificationTick++;
  isCompressedFlag = false;
}

void SurfaceMesh::deleteElement(Face f) {
  size_t iF = f.getIndex();

  fHalfedgeArr[iF] = INVALID_IND;
  nFacesCount--;

  modificationTick++;
  isCompressedFlag = false;
}

void SurfaceMesh::deleteElement(BoundaryLoop bl) {
  size_t iF = boundaryLoopIndToFaceInd(bl.getIndex());

  fHalfedgeArr[iF] = INVALID_IND;
  nBoundaryLoopsCount--;

  modificationTick++;
  isCompressedFlag = false;
}


/*

void SurfaceMesh::deleteElement(Halfedge he) {
  he->markDead();
  isCompressedFlag = false;

  if (he.isReal()) {
    nRealHalfedgesCount--;
  } else {
    nImaginaryHalfedgesCount--;
  }
}

void SurfaceMesh::deleteElement(Edge e) {
  e->markDead();
  isCompressedFlag = false;
  nEdgesCount--;
}

void SurfaceMesh::deleteElement(Vertex v) {
  v->markDead();
  isCompressedFlag = false;
  nVerticesCount--;
}

void SurfaceMesh::deleteElement(Face f) {
  f->markDead();
  isCompressedFlag = false;
  nFacesCount--;
}

void SurfaceMesh::compressHalfedges() {

  // Build the compressing shift
  std::vector<size_t> newIndMap; // maps new ind -> old ind
  std::vector<size_t> oldIndMap; // maps new old -> new ind
  oldIndMap.resize(rawHalfedges.size());
  for (size_t i = 0; i < rawHalfedges.size(); i++) {
    if (!rawHalfedges[i].isDead()) {
      oldIndMap[i] = newIndMap.size();
      newIndMap.push_back(i);
    }
  }
  // === Prep the "before" lists
  Halfedge* oldStart = &rawHalfedges.front();
  std::vector<size_t> offsetsTwin(rawHalfedges.size());
  std::vector<size_t> offsetsNext(rawHalfedges.size());
  std::vector<size_t> offsetsV(rawVertices.size());
  std::vector<size_t> offsetsE(rawEdges.size());
  std::vector<size_t> offsetsF(rawFaces.size());
  std::vector<size_t> offsetsB(rawBoundaryLoops.size());

  for (size_t iHe = 0; iHe < rawHalfedges.size(); iHe++) {
    if (!rawHalfedges[iHe].isDead()) {
      offsetsTwin[iHe] = rawHalfedges[iHe].twin - oldStart;
      offsetsNext[iHe] = rawHalfedges[iHe].next - oldStart;
    }
  }
  for (size_t iV = 0; iV < rawVertices.size(); iV++) {
    if (!rawVertices[iV].isDead()) {
      offsetsV[iV] = rawVertices[iV].halfedge - oldStart;
    }
  }
  for (size_t iE = 0; iE < rawEdges.size(); iE++) {
    if (!rawEdges[iE].isDead()) {
      offsetsE[iE] = rawEdges[iE].halfedge - oldStart;
    }
  }
  for (size_t iF = 0; iF < rawFaces.size(); iF++) {
    if (!rawFaces[iF].isDead()) {
      offsetsF[iF] = rawFaces[iF].halfedge - oldStart;
    }
  }
  for (size_t iB = 0; iB < rawBoundaryLoops.size(); iB++) {
    if (!rawBoundaryLoops[iB].isDead()) {
      offsetsB[iB] = rawBoundaryLoops[iB].halfedge - oldStart;
    }
  }

  // Apply the permutation
  rawHalfedges = applyPermutation(rawHalfedges, newIndMap);
  Halfedge* newStart = &rawHalfedges.front();

  // === Loop back through, shifting all pointers
  // TODO since this is in compress(), should never be dead, right?
  for (size_t iHe = 0; iHe < rawHalfedges.size(); iHe++) {
    if (!rawHalfedges[iHe].isDead()) {
      rawHalfedges[iHe].twin = newStart + oldIndMap[offsetsTwin[newIndMap[iHe]]];
      rawHalfedges[iHe].next = newStart + oldIndMap[offsetsNext[newIndMap[iHe]]];
    }
  }
  for (size_t iV = 0; iV < rawVertices.size(); iV++) {
    if (!rawVertices[iV].isDead()) {
      rawVertices[iV].halfedge = newStart + oldIndMap[offsetsV[iV]];
    }
  }
  for (size_t iE = 0; iE < rawEdges.size(); iE++) {
    if (!rawEdges[iE].isDead()) {
      rawEdges[iE].halfedge = newStart + oldIndMap[offsetsE[iE]];
    }
  }
  for (size_t iF = 0; iF < rawFaces.size(); iF++) {
    if (!rawFaces[iF].isDead()) {
      rawFaces[iF].halfedge = newStart + oldIndMap[offsetsF[iF]];
    }
  }
  for (size_t iB = 0; iB < rawBoundaryLoops.size(); iB++) {
    if (!rawBoundaryLoops[iB].isDead()) {
      rawBoundaryLoops[iB].halfedge = newStart + oldIndMap[offsetsB[iB]];
    }
  }

  // Invoke callbacks
  for (auto& f : halfedgePermuteCallbackList) {
    f(newIndMap);
  }
}

void SurfaceMesh::compressEdges() {

  // Build the compressing shift
  std::vector<size_t> newIndMap; // maps new ind -> old ind
  std::vector<size_t> oldIndMap; // maps new old -> new ind
  oldIndMap.resize(rawEdges.size());
  for (size_t i = 0; i < rawEdges.size(); i++) {
    if (!rawEdges[i].isDead()) {
      oldIndMap[i] = newIndMap.size();
      newIndMap.push_back(i);
    }
  }


  Edge* oldStart = &rawEdges.front();
  std::vector<size_t> offsets(rawHalfedges.size());
  for (size_t iHe = 0; iHe < rawHalfedges.size(); iHe++) {
    if (!rawHalfedges[iHe].isDead()) {
      offsets[iHe] = rawHalfedges[iHe].edge - oldStart;
    }
  }

  // Apply the permutation
  rawEdges = applyPermutation(rawEdges, newIndMap);

  Edge* newStart = &rawEdges.front();
  for (size_t iHe = 0; iHe < rawHalfedges.size(); iHe++) {
    if (!rawHalfedges[iHe].isDead()) {
      rawHalfedges[iHe].edge = newStart + oldIndMap[offsets[iHe]];
    }
  }

  // Invoke callbacks
  for (auto& f : edgePermuteCallbackList) {
    f(newIndMap);
  }
}

void SurfaceMesh::compressFaces() {

  // Build the compressing shift
  std::vector<size_t> newIndMap; // maps new ind -> old ind
  std::vector<size_t> oldIndMap; // maps new old -> new ind
  oldIndMap.resize(rawFaces.size());
  for (size_t i = 0; i < rawFaces.size(); i++) {
    if (!rawFaces[i].isDead()) {
      oldIndMap[i] = newIndMap.size();
      newIndMap.push_back(i);
    }
  }

  Face* oldStart = &rawFaces.front();
  std::vector<size_t> offsets(rawHalfedges.size());
  for (size_t iHe = 0; iHe < rawHalfedges.size(); iHe++) {
    if (!rawHalfedges[iHe].isDead()) {
      offsets[iHe] = rawHalfedges[iHe].face - oldStart;
    }
  }

  // Apply the permutation
  rawFaces = applyPermutation(rawFaces, newIndMap);

  Face* newStart = &rawFaces.front();
  for (size_t iHe = 0; iHe < rawHalfedges.size(); iHe++) {
    if (!rawHalfedges[iHe].isDead() && rawHalfedges[iHe].isReal) {
      rawHalfedges[iHe].face = newStart + oldIndMap[offsets[iHe]];
    }
  }

  // Invoke callbacks
  // validateConnectivity();
  for (auto& f : facePermuteCallbackList) {
    f(newIndMap);
  }
}


void SurfaceMesh::compressVertices() {

  // Build the compressing shift
  std::vector<size_t> newIndMap; // maps new ind -> old ind
  std::vector<size_t> oldIndMap; // maps old ind -> new ind
  oldIndMap.resize(rawVertices.size());
  for (size_t i = 0; i < rawVertices.size(); i++) {
    if (!rawVertices[i].isDead()) {
      oldIndMap[i] = newIndMap.size();
      newIndMap.push_back(i);
    }
  }

  Vertex* oldStart = &rawVertices.front();
  std::vector<size_t> offsets(rawHalfedges.size());
  for (size_t iHe = 0; iHe < rawHalfedges.size(); iHe++) {
    if (!rawHalfedges[iHe].isDead()) {
      offsets[iHe] = rawHalfedges[iHe].vertex - oldStart;
    }
  }

  // Apply the permutation
  rawVertices = applyPermutation(rawVertices, newIndMap);

  Vertex* newStart = &rawVertices.front();
  for (size_t iHe = 0; iHe < rawHalfedges.size(); iHe++) {
    if (!rawHalfedges[iHe].isDead()) {
      rawHalfedges[iHe].vertex = newStart + oldIndMap[offsets[iHe]];
    }
  }

  // Invoke callbacks
  for (auto& f : vertexPermuteCallbackList) {
    f(newIndMap);
  }
}

void SurfaceMesh::compress() {

  if (isCompressed()) {
    return;
  }

  compressHalfedges();
  compressEdges();
  compressFaces();
  compressVertices();
  isCompressedFlag = true;
}


*/


} // namespace surface
} // namespace geometrycentral
