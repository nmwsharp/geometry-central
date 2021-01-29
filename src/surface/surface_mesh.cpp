#include "geometrycentral/surface/surface_mesh.h"

#include "geometrycentral/surface/manifold_surface_mesh.h"
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

  // Sanity check to detect unreferenced vertices
#ifndef NGC_SAFETY_CHECKS
  std::vector<char> vertUsed(nVerticesCount, false);
#endif

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

#ifndef NGC_SAFETY_CHECKS
      vertUsed[indTail] = true;
#endif

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

#ifndef NGC_SAFETY_CHECKS
  // Look for any vertices which were unreferenced
  for (size_t iV = 0; iV < nVerticesCount; iV++) {
    GC_SAFETY_ASSERT(vertUsed[iV], "unreferenced vertex " + std::to_string(iV));
  }
#endif

  // === Create edges and hook up twins
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
          heOrientArr[iHe] = true;
          eHalfedgeArr[newEdgeInd] = iHe;
          edgeHistory[eKey] = iHe;
        } else {
          // We're already seen this edge, connect to the previous halfedge incident on the edge
          size_t iPrevHe = edgeHistory[eKey];
          heSiblingArr[iHe] = iPrevHe;
          size_t iE = heEdgeArr[iPrevHe];
          heEdgeArr[iHe] = iE;
          // best we can to is set orientation to match endpoints (need a richer representation to input orientation if
          // endpoints are not unique)
          heOrientArr[iHe] = (heVertexArr[iHe] == heVertexArr[eHalfedgeArr[iE]]);
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
        heSiblingArr[lastHe] = lastHe;
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
    throw std::runtime_error("not implemented");
  }

  initializeHalfedgeNeighbors();

  isCompressedFlag = true;
  // TODO compress here?
}


SurfaceMesh::SurfaceMesh(const std::vector<size_t>& heNextArr_, const std::vector<size_t>& heVertexArr_,
                         const std::vector<size_t>& heFaceArr_, const std::vector<size_t>& vHalfedgeArr_,
                         const std::vector<size_t>& fHalfedgeArr_, const std::vector<size_t>& heSiblingArr_,
                         const std::vector<size_t>& heEdgeArr_, const std::vector<char>& heOrientArr_,
                         const std::vector<size_t>& eHalfedgeArr_, size_t nBoundaryLoopsFillCount_)
    : heNextArr(heNextArr_), heVertexArr(heVertexArr_), heFaceArr(heFaceArr_), vHalfedgeArr(vHalfedgeArr_),
      fHalfedgeArr(fHalfedgeArr_), useImplicitTwinFlag(false), heSiblingArr(heSiblingArr_), heEdgeArr(heEdgeArr_),
      heOrientArr(heOrientArr_), eHalfedgeArr(eHalfedgeArr_) {

  // == Set all counts
  nHalfedgesCount = heNextArr.size();
  nEdgesCount = eHalfedgeArr.size();
  nVerticesCount = vHalfedgeArr.size();
  nFacesCount = fHalfedgeArr.size() - nBoundaryLoopsFillCount_;
  nBoundaryLoopsCount = nBoundaryLoopsFillCount_;
  nVerticesCapacityCount = nVerticesCount;
  nHalfedgesCapacityCount = nHalfedgesCount;
  nEdgesCapacityCount = nEdgesCount;
  nFacesCapacityCount = fHalfedgeArr.size();
  nVerticesFillCount = nVerticesCount;
  nHalfedgesFillCount = nHalfedgesCount;
  nEdgesFillCount = nEdgesCount;
  nFacesFillCount = nFacesCount;
  nBoundaryLoopsFillCount = nBoundaryLoopsFillCount_;

  // Check if its compressed and decrement counts
  isCompressedFlag = true;
  for (size_t iV = 0; iV < nVerticesFillCount; iV++) {
    if (vertexIsDead(iV)) {
      nVerticesCount--;
      isCompressedFlag = false;
    }
  }
  for (size_t iHe = 0; iHe < nHalfedgesFillCount; iHe++) {
    if (halfedgeIsDead(iHe)) {
      nHalfedgesCount--;
      isCompressedFlag = false;
    }
  }
  for (size_t iE = 0; iE < nEdgesFillCount; iE++) {
    if (edgeIsDead(iE)) {
      nEdgesCount--;
      isCompressedFlag = false;
    }
  }
  for (size_t iF = 0; iF < nFacesFillCount; iF++) {
    if (faceIsDead(iF)) {
      nFacesCount--;
      isCompressedFlag = false;
    }
  }
  for (size_t iBl = nFacesFillCount; iBl < nFacesCapacityCount; iBl++) {
    if (faceIsDead(iBl)) {
      nBoundaryLoopsCount--;
      isCompressedFlag = false;
    }
  }

  // Count interior halfedges
  nInteriorHalfedgesCount = 0;
  for (Halfedge he : interiorHalfedges()) {
    nInteriorHalfedgesCount++;
  }

  initializeHalfedgeNeighbors();
}


SurfaceMesh::~SurfaceMesh() {
  for (auto& f : meshDeleteCallbackList) {
    f();
  }
}

void SurfaceMesh::initializeHalfedgeNeighbors() {

  // NOTE it might be nice to maintain an invariant that around manifold vertices, halfedges appear in manifold order
  // (but currently they're in arbitrary order)

  // We will need to loop around vertex below (on a possibly nonmanifold mesh)
  std::vector<size_t> vertexIterationCacheHeIndexIn;
  std::vector<size_t> vertexIterationCacheVertexStartIn;
  generateVertexIterationCache(vertexIterationCacheHeIndexIn, vertexIterationCacheVertexStartIn, true, true);
  std::vector<size_t> vertexIterationCacheHeIndexOut;
  std::vector<size_t> vertexIterationCacheVertexStartOut;
  generateVertexIterationCache(vertexIterationCacheHeIndexOut, vertexIterationCacheVertexStartOut, false, true);

  heVertInNextArr.resize(nHalfedgesCapacityCount);
  heVertInPrevArr.resize(nHalfedgesCapacityCount);
  vHeInStartArr.resize(nVerticesCapacityCount);
  heVertOutNextArr.resize(nHalfedgesCapacityCount);
  heVertOutPrevArr.resize(nHalfedgesCapacityCount);
  vHeOutStartArr.resize(nVerticesCapacityCount);

  for (Vertex v : vertices()) {

    { // Manually traverse the incoming halfedges incident on this vertex
      size_t rangeStart = vertexIterationCacheVertexStartIn[v.getIndex()];
      size_t rangeEnd = vertexIterationCacheVertexStartIn[v.getIndex() + 1];

      vHeInStartArr[v.getIndex()] = vertexIterationCacheHeIndexIn[rangeStart];

      for (size_t traverseInd = rangeStart; traverseInd < rangeEnd; traverseInd++) {
        size_t iHeA = vertexIterationCacheHeIndexIn[traverseInd];
        size_t iHeB =
            vertexIterationCacheHeIndexIn[(((traverseInd - rangeStart) + 1) % (rangeEnd - rangeStart)) + rangeStart];

        heVertInNextArr[iHeA] = iHeB;
        heVertInPrevArr[iHeB] = iHeA;
      }
    }

    { // Manually traverse the outgoing halfedges incident on this vertex
      size_t rangeStart = vertexIterationCacheVertexStartOut[v.getIndex()];
      size_t rangeEnd = vertexIterationCacheVertexStartOut[v.getIndex() + 1];

      vHeOutStartArr[v.getIndex()] = vertexIterationCacheHeIndexOut[rangeStart];

      for (size_t traverseInd = rangeStart; traverseInd < rangeEnd; traverseInd++) {
        size_t iHeA = vertexIterationCacheHeIndexOut[traverseInd];
        size_t iHeB =
            vertexIterationCacheHeIndexOut[(((traverseInd - rangeStart) + 1) % (rangeEnd - rangeStart)) + rangeStart];

        if (heVertexArr[iHeA] != v.getIndex()) throw std::runtime_error("out A problem");
        if (heVertexArr[iHeB] != v.getIndex()) throw std::runtime_error("out B problem");

        heVertOutNextArr[iHeA] = iHeB;
        heVertOutPrevArr[iHeB] = iHeA;
      }
    }
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


bool SurfaceMesh::hasBoundary() {
  for (Edge e : edges()) {
    if (e.isBoundary()) {
      return true;
    }
  }
  return false;
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

bool SurfaceMesh::isEdgeManifold() {
  for (Edge e : edges()) {
    if (!e.isManifold()) return false;
  }
  return true;
}

bool SurfaceMesh::isOriented() {
  for (Edge e : edges()) {
    if (!e.isManifold()) return false;
    if (!e.isOriented()) return false;
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


std::unique_ptr<SurfaceMesh> SurfaceMesh::copy() const { return copyToSurfaceMesh(); }

std::unique_ptr<SurfaceMesh> SurfaceMesh::copyToSurfaceMesh() const {
  SurfaceMesh* newMesh = new SurfaceMesh(false);
  copyInternalFields(*newMesh);
  return std::unique_ptr<SurfaceMesh>(newMesh);
}


std::unique_ptr<ManifoldSurfaceMesh> SurfaceMesh::toManifoldMesh() {
  if (!isManifold()) throw std::runtime_error("must be manifold to create manifold surface mesh");
  if (!isOriented()) throw std::runtime_error("must be oriented to create manifold surface mesh");

  // Construct buffers for connectivity and build a new manifold mesh
  // NOTE: probably could do this much more efficiently by leveraging the internal representation

  std::vector<std::vector<size_t>> polygons = getFaceVertexList();

  HalfedgeData<size_t> iHeInFace(*this, 0);
  FaceData<size_t> faceInd = getFaceIndices();
  for (Face f : faces()) {
    size_t i = 0;
    for (Halfedge he : f.adjacentHalfedges()) {
      iHeInFace[he] = i;
      i++;
    }
  }

  // Build twin array
  std::vector<std::vector<std::tuple<size_t, size_t>>> twins(nFaces());
  for (Face f : faces()) {
    size_t iF = faceInd[f];
    std::vector<std::tuple<size_t, size_t>>& thisTwin = twins[iF];
    thisTwin.resize(polygons[iF].size());

    size_t i = 0;
    for (Halfedge he : f.adjacentHalfedges()) {
      if (he.edge().isBoundary()) {
        thisTwin[i] = std::make_tuple(INVALID_IND, INVALID_IND);
      } else {
        Halfedge heT = he.sibling();
        size_t oF = faceInd[heT.face()];
        size_t heTInd = iHeInFace[heT];
        thisTwin[i] = std::make_tuple(oF, heTInd);
      }
      i++;
    }
  }

  return std::unique_ptr<ManifoldSurfaceMesh>(new ManifoldSurfaceMesh(polygons, twins));
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
  target.heOrientArr = heOrientArr;
  target.eHalfedgeArr = eHalfedgeArr;
  target.heVertInNextArr = heVertInNextArr;
  target.heVertInPrevArr = heVertInPrevArr;
  target.vHeInStartArr = vHeInStartArr;
  target.heVertOutNextArr = heVertOutNextArr;
  target.heVertOutPrevArr = heVertOutPrevArr;
  target.vHeOutStartArr = vHeOutStartArr;

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


Edge SurfaceMesh::connectingEdge(Vertex vA, Vertex vB) {
  for (Edge e : vA.adjacentEdges()) {
    if (e.otherVertex(vA) == vB) {
      return e;
    }
  }
  return Edge();
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


void SurfaceMesh::generateVertexIterationCache(std::vector<size_t>& vertexIterationCacheHeIndex,
                                               std::vector<size_t>& vertexIterationCacheVertexStart, bool incoming,
                                               bool skipDead) {

  // First, count the degree of every vertex
  std::vector<size_t> vDegree(nVerticesFillCount, 0);
  for (size_t iHe = 0; iHe < nHalfedgesFillCount; iHe++) {
    if (skipDead && halfedgeIsDead(iHe)) continue;
    size_t iV = incoming ? heVertex(heNext(iHe)) : heVertex(iHe);
    vDegree[iV]++;
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
    size_t iV = incoming ? heVertex(heNext(iHe)) : heVertex(iHe);
    size_t entryInd = currVertexCacheEntry[iV];
    vertexIterationCacheHeIndex[entryInd] = iHe;
    currVertexCacheEntry[iV]++;
  }
}


// Note that this uses he.tipVertex()  and he.tailVertex(), so heVert and heNext need to be correct when it is called
void SurfaceMesh::removeFromVertexLists(Halfedge he) {
  size_t i = he.getIndex();
  { // incoming array
    size_t iN = heVertInNextArr[i];
    size_t iP = heVertInPrevArr[i];
    heVertInNextArr[iP] = iN;
    heVertInPrevArr[iN] = iP;
    heVertInNextArr[i] = INVALID_IND;
    heVertInPrevArr[i] = INVALID_IND;
    vHeInStartArr[he.tipVertex().getIndex()] = (iP == i) ? INVALID_IND : iP;
  }

  { // outgoing array
    size_t iN = heVertOutNextArr[i];
    size_t iP = heVertOutPrevArr[i];
    heVertOutNextArr[iP] = iN;
    heVertOutPrevArr[iN] = iP;
    heVertOutNextArr[i] = INVALID_IND;
    heVertOutPrevArr[i] = INVALID_IND;
    vHeOutStartArr[he.tailVertex().getIndex()] = (iP == i) ? INVALID_IND : iP;
  }
}
void SurfaceMesh::addToVertexLists(Halfedge he) {
  size_t i = he.getIndex();

  { // incoming array
    // get any vertex in the current list to use as an insertion point
    size_t iV = he.tipVertex().getIndex();

    size_t iN = vHeInStartArr[iV];
    if (iN == INVALID_IND) {
      // this is the only element in the list (rare)
      heVertInPrevArr[i] = i;
      heVertInNextArr[i] = i;
      vHeInStartArr[iV] = i;
    } else {
      size_t iP = heVertInPrevArr[iN];
      heVertInNextArr[iP] = i;
      heVertInPrevArr[i] = iP;
      heVertInNextArr[i] = iN;
      heVertInPrevArr[iN] = i;
    }
  }

  { // outgoing array
    // get any vertex in the current list to use as an insertion point
    size_t iV = he.tailVertex().getIndex();

    size_t iN = vHeOutStartArr[iV];
    if (iN == INVALID_IND) {
      // this is the only element in the list (rare)
      heVertOutPrevArr[i] = i;
      heVertOutNextArr[i] = i;
      vHeOutStartArr[iV] = i;
    } else {
      size_t iP = heVertOutPrevArr[iN];
      heVertOutNextArr[iP] = i;
      heVertOutPrevArr[i] = iP;
      heVertOutNextArr[i] = iN;
      heVertOutPrevArr[iN] = i;
    }
  }
}

void SurfaceMesh::removeFromSiblingList(Halfedge he) {
  Halfedge heNext = he.sibling();
  Halfedge hePrev = heNext;
  while (hePrev.sibling() != he) {
    hePrev = hePrev.sibling();
  }

  heSiblingArr[hePrev.getIndex()] = heNext.getIndex();
}

void SurfaceMesh::invertOrientation(Face f) {
  if (usesImplicitTwin())
    throw std::runtime_error("Cannot invert orientation on oriented surface. Try a general SurfaceMesh.");

  for (Halfedge he : f.adjacentHalfedges()) removeFromVertexLists(he);

  Halfedge firstHe = f.halfedge();
  Halfedge currHe = firstHe;
  Vertex firstVert = currHe.vertex();
  Halfedge prevHe = Halfedge();
  do {
    // gather values
    Halfedge nextHe = currHe.next();
    Vertex nextVert = (nextHe == firstHe) ? firstVert : nextHe.vertex();

    // update
    heVertexArr[currHe.getIndex()] = nextVert.getIndex();
    vHalfedgeArr[nextVert.getIndex()] = currHe.getIndex();
    heOrientArr[currHe.getIndex()] = !heOrientArr[currHe.getIndex()];
    if (prevHe != Halfedge()) {
      heNextArr[currHe.getIndex()] = prevHe.getIndex();
    }

    // continue looping
    prevHe = currHe;
    currHe = nextHe;
  } while (currHe != firstHe);
  heNextArr[firstHe.getIndex()] = prevHe.getIndex();

  for (Halfedge he : f.adjacentHalfedges()) addToVertexLists(he);
  modificationTick++;
}

Face SurfaceMesh::duplicateFace(Face f) {
  if (usesImplicitTwin())
    throw std::runtime_error("Cannot duplicate a face on a manfiold mesh. Try a general SurfaceMesh.");


  Face newFace = getNewFace();
  bool first = true;
  Halfedge prevNewHe, firstNewHe;
  for (Halfedge oldHe : f.adjacentHalfedges()) {
    Halfedge newHe = getNewHalfedge(false);

    if (first) {
      fHalfedgeArr[newFace.getIndex()] = newHe.getIndex();
      firstNewHe = newHe;
      first = false;
    } else {
      heNextArr[prevNewHe.getIndex()] = newHe.getIndex();
    }

    // update the halfedge data
    heVertexArr[newHe.getIndex()] = oldHe.vertex().getIndex();
    heEdgeArr[newHe.getIndex()] = oldHe.edge().getIndex();
    heOrientArr[newHe.getIndex()] = heOrientArr[oldHe.getIndex()];
    heFaceArr[newHe.getIndex()] = newFace.getIndex();

    // insert in to the sibling list
    size_t sibP = oldHe.getIndex();
    size_t sibN = heSiblingArr[sibP];
    heSiblingArr[sibP] = newHe.getIndex();
    heSiblingArr[newHe.getIndex()] = sibN;
    prevNewHe = newHe;
  }
  heNextArr[prevNewHe.getIndex()] = firstNewHe.getIndex();

  for (Halfedge he : newFace.adjacentHalfedges()) {
    addToVertexLists(he);
  }

  modificationTick++;
  return newFace;
}

bool SurfaceMesh::flip(Edge eFlip, bool preventSelfEdges) {
  if (eFlip.isBoundary()) return false;

  // Get halfedges of first face
  Halfedge ha1 = eFlip.halfedge();
  Halfedge ha2 = ha1.next();
  Halfedge ha3 = ha2.next();
  if (ha3.next() != ha1) return false; // not a triangle

  // Get halfedges of second face
  Halfedge hb1 = ha1.sibling();
  Halfedge hb2 = hb1.next();
  Halfedge hb3 = hb2.next();
  if (hb3.next() != hb1) return false; // not a triangle

  if (hb1.sibling() != ha1) return false;     // not manifold
  if (ha2 == hb1 || hb2 == ha1) return false; // incident on degree 1 vertex

  // if the faces have different orientation, temporarily orient and try again
  if (ha1.orientation() == hb1.orientation()) { // same edge orientation means different face orientation
    invertOrientation(ha1.face());
    bool flipResult = flip(eFlip, preventSelfEdges);
    invertOrientation(ha1.face());
    return flipResult;
  }

  // Get vertices and faces
  Vertex va = ha1.vertex();
  Vertex vb = hb1.vertex();
  Vertex vc = ha3.vertex();
  Vertex vd = hb3.vertex();

  if (preventSelfEdges) {
    // If enabled, make sure it is not a duplicate
    for (Vertex v : vc.adjacentVertices()) {
      if(v == vd) return false;
    }
  }
  
  Face fa = ha1.face();
  Face fb = hb1.face();

  if (!usesImplicitTwin()) {
    removeFromVertexLists(ha1);
    removeFromVertexLists(hb1);
  }

  // Update vertex pointers
  if (va.halfedge() == ha1) vHalfedgeArr[va.getIndex()] = hb2.getIndex();
  if (vb.halfedge() == hb1) vHalfedgeArr[vb.getIndex()] = ha2.getIndex();
  // (vc and vd can't be invalidated by the flip)

  // Update edge pointers
  // (e still has the same halfedges)

  // Update face pointers
  fHalfedgeArr[fa.getIndex()] = ha1.getIndex();
  fHalfedgeArr[fb.getIndex()] = hb1.getIndex();

  // Update halfedge pointers
  heNextArr[ha1.getIndex()] = hb3.getIndex();
  heNextArr[hb3.getIndex()] = ha2.getIndex();
  heNextArr[ha2.getIndex()] = ha1.getIndex();
  heNextArr[hb1.getIndex()] = ha3.getIndex();
  heNextArr[ha3.getIndex()] = hb2.getIndex();
  heNextArr[hb2.getIndex()] = hb1.getIndex();

  heVertexArr[ha1.getIndex()] = vc.getIndex();
  heVertexArr[hb1.getIndex()] = vd.getIndex();

  heFaceArr[ha3.getIndex()] = fb.getIndex();
  heFaceArr[hb3.getIndex()] = fa.getIndex();

  if (!usesImplicitTwin()) {
    addToVertexLists(ha1);
    addToVertexLists(hb1);
  }

  modificationTick++;
  return true;
}

Edge SurfaceMesh::separateToNewEdge(Halfedge heA, Halfedge heB) {
  if (usesImplicitTwin())
    throw std::runtime_error(
        "Cannot separate edge from manifold mesh; all are already manifold. Try general SurfaceMesh.");

  if (heA.edge() != heB.edge()) throw std::runtime_error("halfedges must be incident on same edge");
  if (heA == heB) throw std::runtime_error("halfedges must be distinct");

  // If theres <= 2 halfedges incident on the edge, we don't have any work to do.
  Edge e = heA.edge();
  if (e.degree() <= 2) return e;

  Edge newE = getNewEdge();

  // find some other halfedge incident on the old edge, make it e.halfedge()
  for (Halfedge he : e.adjacentHalfedges()) {
    if (he != heA && he != heB) {
      eHalfedgeArr[e.getIndex()] = he.getIndex();
      break;
    }
  }

  removeFromSiblingList(heA);
  removeFromSiblingList(heB);

  eHalfedgeArr[newE.getIndex()] = heA.getIndex();
  heEdgeArr[heA.getIndex()] = newE.getIndex();
  heEdgeArr[heB.getIndex()] = newE.getIndex();
  heSiblingArr[heA.getIndex()] = heB.getIndex();
  heSiblingArr[heB.getIndex()] = heA.getIndex();


  modificationTick++;
  return newE;
}

void SurfaceMesh::separateNonmanifoldEdges() {

  for (Edge e : edges()) {
    while (!e.isManifold()) {
      Halfedge heA = e.halfedge();
      Halfedge heB = heA.sibling();
      separateToNewEdge(heA, heB);
    }
  }

  modificationTick++;
}


VertexData<Vertex> SurfaceMesh::separateNonmanifoldVertices() {

  // Find edge-connected sets of corners
  size_t indMax = nHalfedgesFillCount;
  DisjointSets djSet(indMax);
  for (Edge e : edges()) {
    if (e.isBoundary()) continue;
    if (!e.isManifold()) {
      throw std::runtime_error("mesh must be edge-manifold for separateNonmanifoldVertices()");
    }
    Halfedge heA = e.halfedge();
    Halfedge heB = heA.sibling();

    if (heA.orientation() == heB.orientation()) {
      djSet.merge(heA.corner().getIndex(), heB.corner().getIndex());
      djSet.merge(heA.next().corner().getIndex(), heB.next().corner().getIndex());
    } else {
      djSet.merge(heA.next().corner().getIndex(), heB.corner().getIndex());
      djSet.merge(heA.corner().getIndex(), heB.next().corner().getIndex());
    }
  }

  // Make sure there is a distinct vertex entry for each component
  VertexData<Vertex> parents(*this);
  std::vector<Vertex> vertexEntries(indMax, Vertex());
  VertexData<bool> baseVertexUsed(*this, false);
  for (Corner c : corners()) {
    size_t iComp = djSet.find(c.getIndex());
    Vertex origV = c.vertex();

    // create vertex if needed
    if (vertexEntries[iComp] == Vertex()) {
      if (baseVertexUsed[origV]) {
        Vertex newV = getNewVertex();
        vertexEntries[iComp] = newV;
        parents[newV] = origV;
      } else {
        parents[origV] = origV;
        vertexEntries[iComp] = origV;
        baseVertexUsed[origV] = true;
      }
    }

    // hook up
    Vertex targetV = vertexEntries[iComp];
    Halfedge he = c.halfedge();
    heVertexArr[he.getIndex()] = targetV.getIndex();
    vHalfedgeArr[targetV.getIndex()] = he.getIndex();
  }

  // just rebuild these from scratch, rather than trying to maintain
  initializeHalfedgeNeighbors();

  modificationTick++;
  return parents;
}

void SurfaceMesh::greedilyOrientFaces() {
  // TODO this is only lightly tested. Write some better tests.
  std::vector<Face> toProcess;
  FaceData<double> processed(*this, false);
  for (Face f : faces()) {
    if (processed[f]) continue;

    // start a new search
    toProcess.push_back(f);
    processed[f] = true;

    // cover the whole component
    while (!toProcess.empty()) {
      Face curr = toProcess.back();
      toProcess.pop_back();

      for (Halfedge he : curr.adjacentHalfedges()) {
        if (he.edge().isBoundary() || !he.edge().isManifold()) continue;
        Face fO = he.sibling().face();
        if (processed[fO]) continue;

        if (!he.edge().isOriented()) {
          invertOrientation(fO);
        }

        toProcess.push_back(fO);
        processed[fO] = true;
      }
    }
  }
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

  // Check for overflow / other unreasonable values
  if (nHalfedgesCount > std::numeric_limits<uint64_t>::max() / 2) throw std::logic_error("halfedge count overflow");
  if (nExteriorHalfedges() > std::numeric_limits<uint64_t>::max() / 2)
    throw std::logic_error("exterior halfedge count overflow");
  if (nVerticesCount > std::numeric_limits<uint64_t>::max() / 2) throw std::logic_error("vertex count overflow");
  if (nEdgesCount > std::numeric_limits<uint64_t>::max() / 2) throw std::logic_error("edge count overflow");
  if (nFacesCount > std::numeric_limits<uint64_t>::max() / 2) throw std::logic_error("face count overflow");

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
  // Note: we intentionally mostly avoid using iterators here, because they can be hard to debug when things are broken.
  for (size_t iHe = 0; iHe < nHalfedgesFillCount; iHe++) {
    if (halfedgeIsDead(iHe)) continue;
    validateHalfedge(heTwin(iHe), "he.twin()");
    validateHalfedge(heNextArr[iHe], "he.next()");
    validateVertex(heVertexArr[iHe], "he.vertex()");
    validateEdge(heEdge(iHe), "he.edge()");
    validateFace(heFaceArr[iHe], "he.face()");
    if (!usesImplicitTwin()) {
      validateHalfedge(heVertInNextArr[iHe], "heVertInNextArr");
      validateHalfedge(heVertInPrevArr[iHe], "heVertInPrevArr");
      validateHalfedge(heVertOutNextArr[iHe], "heVertOutNextArr");
      validateHalfedge(heVertOutPrevArr[iHe], "heVertOutPrevArr");
    }
  }
  for (size_t iV = 0; iV < nVerticesFillCount; iV++) {
    if (vertexIsDead(iV)) continue;
    validateHalfedge(vHalfedgeArr[iV], "v.halfedge()");
    if (!usesImplicitTwin()) {
      validateHalfedge(vHeInStartArr[iV], "vHeInStartArr");
      validateHalfedge(vHeOutStartArr[iV], "vHeOutStartArr");
    }
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
    // allowed at boundary of non-implicit
    // if (he == he.sibling()) throw std::logic_error("self-sibling");

    // Check sibling orbit sanity
    Halfedge currHe = he;
    Halfedge firstHe = he;
    size_t count = 0;
    do {
      if (currHe.edge() != he.edge())
        throw std::logic_error("(he sibling) halfedge sibling doesn't have edge == he.edge");
      if (count > nHalfedges()) throw std::logic_error("(he sibling) halfedge sibling doesn't cycle back");
      currHe = currHe.sibling();
      count++;
    } while (currHe != firstHe);
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


    // Check that boundary rules are observed
    if (!he.isInterior()) {
      if (he == he.edge().halfedge()) throw std::logic_error("exterior halfedge is e.halfedge()");
      if (!he.twin().isInterior()) throw std::logic_error("he and he.twin() are both exterior");
    }

    // This can happen in irregular triangulations
    // if (he.vertex == he.next->twin->vertex) throw std::logic_error("halfedge face spur");

    // Check that he.twin().next() locally orbts a vertex
    if (useImplicitTwinFlag) {
      if (he.vertex() != he.twin().next().vertex()) throw std::logic_error("halfedge vertices don't match");

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
  }

  // Check vertex iteration lists
  if (!usesImplicitTwin()) {
    for (Halfedge he : halfedges()) {
      size_t iHe = he.getIndex();
      Vertex thisTail = he.vertex();
      Vertex thisTip = he.next().vertex();
      if (Halfedge(this, heVertOutNextArr[iHe]).vertex() != thisTail)
        throw std::logic_error("heVertOutNextArr is not outgoing from same vert");
      if (Halfedge(this, heVertOutPrevArr[iHe]).vertex() != thisTail)
        throw std::logic_error("heVertOutPrevArr is not outgoing from same vert");
      if (Halfedge(this, heVertInNextArr[iHe]).next().vertex() != thisTip)
        throw std::logic_error("heVertInNextArr is not incoming from same vert");
      if (Halfedge(this, heVertInPrevArr[iHe]).next().vertex() != thisTip)
        throw std::logic_error("heVertInPrevArr is not incoming from same vert");
    }

    for (Vertex v : vertices()) {
      size_t iIn = vHeInStartArr[v.getIndex()];
      if (Halfedge(this, iIn).tipVertex() != v)
        throw std::logic_error("vHeInStartArr[v] is not incoming from this vertex");

      size_t iOut = vHeOutStartArr[v.getIndex()];
      if (Halfedge(this, iOut).tailVertex() != v)
        throw std::logic_error("vHeOutStartArr[v] is not outgoing from this vertex");
    }
  }


  // Check vertex orbit sanity
  for (Vertex v : vertices()) {
    size_t count = 0;
    for (Halfedge currHe : v.outgoingHalfedges()) {
      if (count > nHalfedgesCount) throw std::logic_error("vertex outgoing halfedges has bad cycle");
      if (currHe.vertex() != v) throw std::logic_error("vertex.halfedge doesn't match halfedge.vertex");
      count++;
    }
  }

  // Verify boundary rules are correct (non-implicit doesn't really have any rules)
  if (usesImplicitTwin()) {
    for (Vertex v : vertices()) {

      // Manually check if this is a boundary vertex
      size_t boundaryHeCount = 0;
      for (Halfedge he : v.outgoingHalfedges()) {
        if (!he.isInterior()) boundaryHeCount++;
      }

      if (useImplicitTwinFlag) {
        if (boundaryHeCount > 1) {
          throw std::logic_error("multiple boundaries incident on vertex");
        }
      }
      bool hasBoundaryHe = boundaryHeCount > 0;

      if (hasBoundaryHe) {
        if (useImplicitTwinFlag) {
          if (!v.halfedge().isInterior()) {
            throw std::logic_error("v.halfedge() is exterior");
          }
          if (v.halfedge().twin().isInterior()) {
            throw std::logic_error("v.halfedge() does not border boundary on a boundary vertex");
          }
        }
        if (!v.isBoundary()) {
          throw std::logic_error("computed v.isBoundary is wrong");
        }
      }
    }
  }

  // Check orientation
  if (!usesImplicitTwin()) {
    for (Halfedge he : halfedges()) {
      Halfedge heCan = he.edge().halfedge();

      if (he.orientation() == heCan.orientation()) {
        if (he.tipVertex() != heCan.tipVertex() || he.tailVertex() != heCan.tailVertex())
          throw std::logic_error("orientation is inconsistent with endpoints");
      } else {
        if (he.tipVertex() != heCan.tailVertex() || he.tailVertex() != heCan.tipVertex())
          throw std::logic_error("orientation is inconsistent with endpoints");
      }
    }
  }

  // Make sure that vertex iterators do what you would expect
  if (!usesImplicitTwin()) {

    std::vector<size_t> vertexInCount(nVerticesCapacityCount, 0);
    std::vector<size_t> vertexOutCount(nVerticesCapacityCount, 0);
    for (Halfedge he : halfedges()) {
      vertexOutCount[he.vertex().getIndex()]++;
      vertexInCount[he.next().vertex().getIndex()]++;
    }

    for (Vertex v : vertices()) {
      size_t nOut = 0;
      for (Halfedge he : v.outgoingHalfedges()) {
        nOut++;
      }
      if (nOut != vertexOutCount[v.getIndex()])
        throw std::logic_error(
            "not enough halfedges in vertex outgoing loop, must be disconnected component in heVertOutNextArr/Prev");

      size_t nIn = 0;
      for (Halfedge he : v.incomingHalfedges()) {
        nIn++;
      }
      if (nIn != vertexInCount[v.getIndex()])
        throw std::logic_error(
            "not enough halfedges in vertex incoming loop, must be disconnected component in heVertInNextArr/Prev");
    }
  }

  // Check manifoldness, if the mesh should be manifold
  if (usesImplicitTwin()) {
    // Right now (and for the forseeable future) the usesImplicitTwin is really equivalent to being of class
    // ManifoldSurfaceMesh, so this error prints that, which is easier to understand to the user.
    if (!isManifold()) {
      throw std::logic_error("Mesh with underlying type ManifoldSurfaceMesh is not actually manifold, likely due to a "
                             "nonmanifold vertex with disconnected neighborhoods.");
    }
  }
} // namespace surface


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
    if (!usesImplicitTwin()) {
      vHeInStartArr.resize(newCapacity);
      vHeOutStartArr.resize(newCapacity);
    }

    nVerticesCapacityCount = newCapacity;


    // Invoke relevant callback functions
    for (auto& f : vertexExpandCallbackList) {
      f(newCapacity);
    }
  }

  nVerticesFillCount++;
  nVerticesCount++;

  modificationTick++;
  isCompressedFlag = false;
  return Vertex(this, nVerticesFillCount - 1);
}

Halfedge SurfaceMesh::getNewHalfedge(bool isInterior) {

  if (usesImplicitTwin()) {
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
    if (!usesImplicitTwin()) { // must enter this case, see test above
      heSiblingArr.resize(newHalfedgeCapacity);
      heEdgeArr.resize(newHalfedgeCapacity);
      heOrientArr.resize(newHalfedgeCapacity);
      heVertInNextArr.resize(newHalfedgeCapacity);
      heVertInPrevArr.resize(newHalfedgeCapacity);
      heVertOutNextArr.resize(newHalfedgeCapacity);
      heVertOutPrevArr.resize(newHalfedgeCapacity);
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
  isCompressedFlag = false;
  return Halfedge(this, nHalfedgesFillCount - 1);
}

Edge SurfaceMesh::getNewEdge() {

  if (usesImplicitTwin()) {
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

    if (!usesImplicitTwin()) { // must enter this case, see test above
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
  isCompressedFlag = false;
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
      if (!usesImplicitTwin()) {
        heSiblingArr.resize(newHalfedgeCapacity);
        heEdgeArr.resize(newHalfedgeCapacity);
        heOrientArr.resize(newHalfedgeCapacity);
      }

      nHalfedgesCapacityCount = newHalfedgeCapacity;

      // Invoke relevant callback functions
      for (auto& f : halfedgeExpandCallbackList) {
        f(newHalfedgeCapacity);
      }
    }

    { // expand edges
      nEdgesCapacityCount = newEdgeCapacity;

      if (!usesImplicitTwin()) {
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
  if (!usesImplicitTwin()) {
    heSiblingArr[nHalfedgesFillCount] = nHalfedgesFillCount + 1;
    heSiblingArr[nHalfedgesFillCount + 1] = nHalfedgesFillCount;
    heEdgeArr[nHalfedgesFillCount] = nEdgesFillCount;
    heEdgeArr[nHalfedgesFillCount + 1] = nEdgesFillCount;
    heOrientArr[nHalfedgesFillCount] = true;
    heOrientArr[nHalfedgesFillCount + 1] = false;
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
  isCompressedFlag = false;
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
  isCompressedFlag = false;
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
  isCompressedFlag = false;
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


void SurfaceMesh::deleteEdgeBundle(Edge e) {

  // TODO there are some pitfalls here, because this method needs to test whether the incident halfedges are interior,
  // but that seemingly means adjacent faces must not have been deleted yet. It seems easy to make a mistake when
  // deleting lots of elements. Need to think through this carefully and document the "safe" way to perform deletions.

  std::vector<size_t> heInds;
  for (Halfedge he : e.adjacentHalfedges()) {
    heInds.push_back(he.getIndex());
  }

  // delete all of the incident halfedges
  for (size_t i : heInds) {
    nHalfedgesCount--;
    if (heIsInterior(i)) {
      nInteriorHalfedgesCount--;
    }

    heNextArr[i] = INVALID_IND;
    heVertexArr[i] = INVALID_IND;
    heFaceArr[i] = INVALID_IND;

    if (!usesImplicitTwin()) {
      heSiblingArr[i] = INVALID_IND;
      heEdgeArr[i] = INVALID_IND;
      heVertInNextArr[i] = INVALID_IND;
      heVertInPrevArr[i] = INVALID_IND;
      heVertOutNextArr[i] = INVALID_IND;
      heVertOutPrevArr[i] = INVALID_IND;
    }
  }

  // delete edge stuff
  if (!usesImplicitTwin()) {
    eHalfedgeArr[e.getIndex()] = INVALID_IND;
  }
  nEdgesCount--;

  modificationTick++;
  isCompressedFlag = false;
}

void SurfaceMesh::deleteElement(Halfedge he) {
  GC_SAFETY_ASSERT(!usesImplicitTwin(), "cannot delete a single halfedge with implict twin");

  // delete all of the incident halfedges
  size_t i = he.getIndex();
  heNextArr[i] = INVALID_IND;
  heVertexArr[i] = INVALID_IND;
  heFaceArr[i] = INVALID_IND;
  heSiblingArr[i] = INVALID_IND;
  heEdgeArr[i] = INVALID_IND;
  heOrientArr[i] = false;
  heVertInNextArr[i] = INVALID_IND;
  heVertInPrevArr[i] = INVALID_IND;
  heVertOutNextArr[i] = INVALID_IND;
  heVertOutPrevArr[i] = INVALID_IND;

  nHalfedgesCount--;
  if (heIsInterior(i)) {
    nInteriorHalfedgesCount--;
  }

  modificationTick++;
  isCompressedFlag = false;
}

void SurfaceMesh::deleteElement(Edge e) {
  GC_SAFETY_ASSERT(!usesImplicitTwin(), "cannot delete a single edge with implict twin");

  size_t i = e.getIndex();
  eHalfedgeArr[e.getIndex()] = INVALID_IND;
  nEdgesCount--;

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

void SurfaceMesh::updateValues(std::vector<size_t>& arr, const std::vector<size_t>& oldToNew) {
  for (size_t& x : arr) {
    if (x == INVALID_IND) continue;
    x = oldToNew[x];
  }
}

void SurfaceMesh::compressHalfedges() {

  // Build the compressing shift
  std::vector<size_t> newIndMap;                                   // maps new ind -> old ind
  std::vector<size_t> newIndEdgeMap;                               // maps edge new ind -> old ind
  std::vector<size_t> oldIndMap(nHalfedgesFillCount, INVALID_IND); // maps old ind -> new ind
  for (size_t i = 0; i < nHalfedgesFillCount; i++) {
    if (!halfedgeIsDead(i)) {
      oldIndMap[i] = newIndMap.size();
      newIndMap.push_back(i);

      if (usesImplicitTwin() && i % 2 == 0) {
        size_t iEdge = i / 2;
        newIndEdgeMap.push_back(iEdge);
      }
    }
  }

  // Permute & resize all per-halfedge arrays
  heNextArr = applyPermutation(heNextArr, newIndMap);
  heVertexArr = applyPermutation(heVertexArr, newIndMap);
  heFaceArr = applyPermutation(heFaceArr, newIndMap);
  if (!usesImplicitTwin()) {
    heSiblingArr = applyPermutation(heSiblingArr, newIndMap);
    heEdgeArr = applyPermutation(heEdgeArr, newIndMap);
    heOrientArr = applyPermutation(heOrientArr, newIndMap);
    heVertInNextArr = applyPermutation(heVertInNextArr, newIndMap);
    heVertInPrevArr = applyPermutation(heVertInPrevArr, newIndMap);
    heVertOutNextArr = applyPermutation(heVertOutNextArr, newIndMap);
    heVertOutPrevArr = applyPermutation(heVertOutPrevArr, newIndMap);
  }


  // Update indices in all halfedge-valued arrays
  updateValues(vHalfedgeArr, oldIndMap);
  updateValues(fHalfedgeArr, oldIndMap);
  updateValues(heNextArr, oldIndMap);
  if (!usesImplicitTwin()) {
    updateValues(eHalfedgeArr, oldIndMap);
    updateValues(heSiblingArr, oldIndMap);
    updateValues(heVertInNextArr, oldIndMap);
    updateValues(heVertInPrevArr, oldIndMap);
    updateValues(vHeInStartArr, oldIndMap);
    updateValues(heVertOutNextArr, oldIndMap);
    updateValues(heVertOutPrevArr, oldIndMap);
    updateValues(vHeOutStartArr, oldIndMap);
  }

  // Update counts
  nHalfedgesFillCount = nHalfedgesCount;
  nHalfedgesCapacityCount = nHalfedgesCount;

  // Invoke callbacks
  for (auto& f : halfedgePermuteCallbackList) {
    f(newIndMap);
  }

  // In the implicit-twin case, we also need to update edge data here, because they are always in-sync with halfedges
  if (usesImplicitTwin()) {
    nEdgesFillCount = nEdgesCount;
    nEdgesCapacityCount = nEdgesCount;
    // Invoke callbacks
    for (auto& f : edgePermuteCallbackList) {
      f(newIndEdgeMap);
    }
  }
}


void SurfaceMesh::compressEdges() {

  if (usesImplicitTwin()) {
    // In the implicit-twin case, all updates are handled in the halfedge function (see note there)
    return;
  }

  // Build the compressing shift
  std::vector<size_t> newIndMap;                               // maps new ind -> old ind
  std::vector<size_t> oldIndMap(nEdgesFillCount, INVALID_IND); // maps old ind -> new ind
  for (size_t i = 0; i < nEdgesFillCount; i++) {
    if (!edgeIsDead(i)) {
      oldIndMap[i] = newIndMap.size();
      newIndMap.push_back(i);
    }
  }

  // Permute & resize all per-edge arrays
  eHalfedgeArr = applyPermutation(eHalfedgeArr, newIndMap);

  // Update indices in all edge-valued arrays
  updateValues(heEdgeArr, oldIndMap);

  // Update counts
  nEdgesFillCount = nEdgesCount;
  nEdgesCapacityCount = nEdgesCount;

  // Invoke callbacks
  for (auto& f : edgePermuteCallbackList) {
    f(newIndMap);
  }
}

void SurfaceMesh::compressFaces() {
  // Build the compressing shift
  std::vector<size_t> newIndMap;                                   // maps new ind -> old ind
  std::vector<size_t> newBLIndMap;                                 // maps BL new ind -> old ind
  std::vector<size_t> oldIndMap(nFacesCapacityCount, INVALID_IND); // maps old ind -> new ind
  for (size_t i = 0; i < nFacesCapacityCount; i++) {
    bool isBL = (i >= nFacesCapacityCount - nBoundaryLoopsFillCount);
    if (i < nFacesFillCount || isBL) { // skip gap between faces and BLs
      if (!faceIsDead(i)) {
        oldIndMap[i] = newIndMap.size();
        newIndMap.push_back(i);

        if (isBL) {
          newBLIndMap.push_back(faceIndToBoundaryLoopInd(i));
        }
      }
    }
  }


  // Permute & resize all per-face arrays
  fHalfedgeArr = applyPermutation(fHalfedgeArr, newIndMap);

  // Update indices in all face-valued arrays
  updateValues(heFaceArr, oldIndMap);

  // Update counts
  nFacesFillCount = nFacesCount;
  nFacesCapacityCount = nFacesCount + nBoundaryLoopsCount;
  nBoundaryLoopsFillCount = nBoundaryLoopsCount;

  // Invoke callbacks
  newIndMap.resize(nFacesCount); // truncate all boundary loop entries, so this array now just holds values for faces
  for (auto& f : facePermuteCallbackList) {
    f(newIndMap);
  }
  for (auto& f : boundaryLoopPermuteCallbackList) {
    f(newBLIndMap);
  }
}


void SurfaceMesh::compressVertices() {
  // Build the compressing shift
  std::vector<size_t> newIndMap;                                  // maps new ind -> old ind
  std::vector<size_t> oldIndMap(nVerticesFillCount, INVALID_IND); // maps old ind -> new ind
  for (size_t i = 0; i < nVerticesFillCount; i++) {
    if (!vertexIsDead(i)) {
      oldIndMap[i] = newIndMap.size();
      newIndMap.push_back(i);
    }
  }


  // Permute & resize all per-vertex arrays
  vHalfedgeArr = applyPermutation(vHalfedgeArr, newIndMap);
  if (!usesImplicitTwin()) {
    vHeInStartArr = applyPermutation(vHeInStartArr, newIndMap);
    vHeOutStartArr = applyPermutation(vHeOutStartArr, newIndMap);
  }

  // Update indices in all vertex-valued arrays
  updateValues(heVertexArr, oldIndMap);

  // Update counts
  nVerticesFillCount = nVerticesCount;
  nVerticesCapacityCount = nVerticesCount;

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


} // namespace surface
} // namespace geometrycentral
