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

namespace geometrycentral {
namespace surface {

ManifoldSurfaceMesh::ManifoldSurfaceMesh() : SurfaceMesh(true) {}

ManifoldSurfaceMesh::ManifoldSurfaceMesh(const std::vector<std::vector<size_t>>& polygons) : SurfaceMesh(true) {
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

  // Sanity check to detect unreferenced vertices
#ifndef NGC_SAFETY_CHECKS
  std::vector<char> vertUsed(nVerticesCount, false);
#endif

  // Track halfedges which have already been created
  // TODO replace with compressed list for performance
  std::unordered_map<std::tuple<size_t, size_t>, size_t> createdHalfedges;
  auto createdHeLookup = [&](std::tuple<size_t, size_t> key) -> size_t& {
    if (createdHalfedges.find(key) == createdHalfedges.end()) {
      createdHalfedges[key] = INVALID_IND;
    }
    return createdHalfedges[key];
  };

  // Walk the faces, creating halfedges and hooking up pointers
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
      size_t& halfedgeInd = createdHeLookup(heKey);

      // Some sanity checks
      GC_SAFETY_ASSERT(indTail != indTip,
                       "self-edge in face list " + std::to_string(indTail) + " -- " + std::to_string(indTip));
      GC_SAFETY_ASSERT(halfedgeInd == INVALID_IND,
                       "duplicate edge in list " + std::to_string(indTail) + " -- " + std::to_string(indTip));

      // Find the twin to check if the element is already created
      size_t twinInd = createdHeLookup(heTwinKey);
      if (twinInd == INVALID_IND) {
        // If we haven't seen the twin yet either, create a new edge
        halfedgeInd = getNewEdgeTriple(false).getIndex();

        // Fill arrays with nknown values and placeholders
        heNextArr[halfedgeInd] = INVALID_IND;
        heNextArr[heTwin(halfedgeInd)] = INVALID_IND;
        heVertexArr[halfedgeInd] = indTail;
        heVertexArr[heTwin(halfedgeInd)] = indTip;
        heFaceArr[halfedgeInd] = INVALID_IND;
        heFaceArr[heTwin(halfedgeInd)] = INVALID_IND;
      } else {
        // If the twin has already been created, we have an index for the halfedge
        halfedgeInd = heTwin(twinInd);
      }

      // Hook up a bunch of pointers
      heFaceArr[halfedgeInd] = iFace;
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

  // Ensure that each boundary neighborhood is either a disk or a half-disk. Harder to diagnose if we wait until the
  // boundary walk below.
  {
    std::vector<char> vertexOnBoundary(nVerticesCount, false);
    for (size_t iHe = 0; iHe < nHalfedgesCount; iHe++) {
      if (heNextArr[iHe] == INVALID_IND) {
        size_t v = heVertexArr[iHe];
        GC_SAFETY_ASSERT(!vertexOnBoundary[v],
                         "vertex " + std::to_string(v) + " appears in more than one boundary loop");
        vertexOnBoundary[v] = true;
      }
    }
  }
#endif

  // == Resolve boundary loops
  nInteriorHalfedgesCount = nHalfedgesCount; // will decrement as we find exterior
  for (size_t iHe = 0; iHe < nHalfedgesCount; iHe++) {

    // If the face pointer is invalid, the halfedge must be along an unresolved boundary loop
    if (heFaceArr[iHe] != INVALID_IND) continue;

    // Create the new boundary loop
    size_t boundaryLoopInd = nFacesCount + nBoundaryLoopsCount;
    fHalfedgeArr.push_back(iHe);
    nBoundaryLoopsCount++;

    // = Walk around the loop (CW)
    size_t currHe = iHe;
    size_t prevHe = INVALID_IND;
    size_t loopCount = 0;
    do {

      // The boundary loop is the face for these halfedges
      heFaceArr[currHe] = boundaryLoopInd;

      // currHe.twin() is a boundary interior halfedge, this is a good time to enforce that v.halfedge() is always the
      // boundary interior halfedge for a boundary vertex.
      size_t currHeT = heTwin(currHe);
      vHalfedgeArr[heVertexArr[currHeT]] = currHeT;

      // This isn't an interior halfedge.
      nInteriorHalfedgesCount--;

      // Advance to the next halfedge along the boundary
      prevHe = currHe;
      currHe = heTwin(heNextArr[heTwin(currHe)]);
      size_t loopCountInnter = 0;
      while (heFaceArr[currHe] != INVALID_IND) {
        if (currHe == iHe) break;
        currHe = heTwin(heNextArr[currHe]);
        loopCountInnter++;
        GC_SAFETY_ASSERT(loopCountInnter < nHalfedgesCount, "boundary infinite loop orbit");
      }

      // Set the next pointer around the boundary loop
      heNextArr[currHe] = prevHe;

      // Make sure this loop doesn't infinite-loop. Certainly won't happen for proper input, but might happen for bogus
      // input. I don't _think_ it can happen, but there might be some non-manifold input which manfests failure via an
      // infinte loop here, and such a loop is an inconvenient failure mode.
      loopCount++;
      GC_SAFETY_ASSERT(loopCount < nHalfedgesCount, "boundary infinite loop");
    } while (currHe != iHe);
  }

  // SOMEDAY: could shrink_to_fit() std::vectors here, at the cost of a copy. What's preferable?

  // Set capacities and other properties
  nVerticesCapacityCount = nVerticesCount;
  nHalfedgesCapacityCount = nHalfedgesCount;
  nEdgesCapacityCount = nEdgesCount;
  nFacesCapacityCount = nFacesCount + nBoundaryLoopsCount;
  nVerticesFillCount = nVerticesCount;
  nHalfedgesFillCount = nHalfedgesCount;
  nEdgesFillCount = nEdgesCount;
  nFacesFillCount = nFacesCount;
  nBoundaryLoopsFillCount = nBoundaryLoopsCount;
  isCompressedFlag = true;

#ifndef NGC_SAFETY_CHECKS
  { // Check that the input was manifold in the sense that each vertex has a single connected loop of faces around it.
    std::vector<char> halfedgeSeen(nHalfedgesCount, false);
    for (size_t iV = 0; iV < nVerticesCount; iV++) {

      // For each vertex, orbit around the outgoing halfedges. This _should_ touch every halfedge.
      size_t currHe = vHalfedgeArr[iV];
      size_t firstHe = currHe;
      do {

        GC_SAFETY_ASSERT(!halfedgeSeen[currHe], "somehow encountered outgoing halfedge before orbiting v");
        halfedgeSeen[currHe] = true;

        currHe = heNextArr[heTwinImplicit(currHe)];
      } while (currHe != firstHe);
    }

    // Verify that we actually did touch every halfedge.
    for (size_t iHe = 0; iHe < nHalfedgesCount; iHe++) {
      GC_SAFETY_ASSERT(halfedgeSeen[iHe], "mesh not manifold. Vertex " + std::to_string(heVertexArr[iHe]) +
                                              " has disconnected neighborhoods incident (imagine an hourglass)");
    }
  }
#endif


  // Print some nice statistics
  // printStatistics();
  // std::cout << "Construction took " << pretty_time(FINISH_TIMING(construction)) << std::endl;
}

ManifoldSurfaceMesh::ManifoldSurfaceMesh(const std::vector<std::vector<size_t>>& polygons,
                                         const std::vector<std::vector<std::tuple<size_t, size_t>>>& twins)
    : SurfaceMesh(true) {

  // Assumes that the input index set is dense. This sometimes isn't true of (eg) obj files floating around the
  // internet, so consider removing unused vertices first when reading from foreign sources.

  START_TIMING(construction)

  GC_SAFETY_ASSERT(polygons.size() == twins.size(), "twin list should be same shape as polygon list");

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

  // NOTE IMPORTANT DIFFERENCE: in the first face-only constructor, these keys are (vInd, vInd) pairs, but here they
  // are (fInd, heInFInd) pairs.

  // Track halfedges which have already been created
  // TODO replace with compressed list for performance
  std::unordered_map<std::tuple<size_t, size_t>, size_t> createdHalfedges;
  auto createdHeLookup = [&](std::tuple<size_t, size_t> key) -> size_t& {
    if (createdHalfedges.find(key) == createdHalfedges.end()) {
      createdHalfedges[key] = INVALID_IND;
    }
    return createdHalfedges[key];
  };

  // Walk the faces, creating halfedges and hooking up pointers
  for (size_t iFace = 0; iFace < nFacesCount; iFace++) {
    const std::vector<size_t>& poly = polygons[iFace];
    const std::vector<std::tuple<size_t, size_t>>& polyTwin = twins[iFace];
    GC_SAFETY_ASSERT(poly.size() == polyTwin.size(), "twin list should be same shape as polygon list");

    // Walk around this face
    size_t faceDegree = poly.size();
    size_t prevHeInd = INVALID_IND;
    size_t firstHeInd = INVALID_IND;
    for (size_t iFaceHe = 0; iFaceHe < faceDegree; iFaceHe++) {

      size_t indTail = poly[iFaceHe];
      size_t indTip = poly[(iFaceHe + 1) % faceDegree];

      // Get an index for this halfedge
      std::tuple<size_t, size_t> heKey{iFace, iFaceHe};
      std::tuple<size_t, size_t> heTwinKey{std::get<0>(polyTwin[iFaceHe]), std::get<1>(polyTwin[iFaceHe])};
      size_t& halfedgeInd = createdHeLookup(heKey);

      // Some sanity checks
      GC_SAFETY_ASSERT(indTail != indTip,
                       "self-edge in face list " + std::to_string(indTail) + " -- " + std::to_string(indTip));
      GC_SAFETY_ASSERT(halfedgeInd == INVALID_IND,
                       "duplicate edge in list " + std::to_string(indTail) + " -- " + std::to_string(indTip));

      // Find the twin to check if the element is already created
      size_t twinInd = createdHeLookup(heTwinKey);
      if (twinInd == INVALID_IND) {
        // If we haven't seen the twin yet either, create a new edge
        halfedgeInd = getNewEdgeTriple(false).getIndex();

        // Fill arrays with nknown values and placeholders
        heNextArr[halfedgeInd] = INVALID_IND;
        heNextArr[heTwin(halfedgeInd)] = INVALID_IND;
        heVertexArr[halfedgeInd] = indTail;
        heVertexArr[heTwin(halfedgeInd)] = indTip;
        heFaceArr[halfedgeInd] = INVALID_IND;
        heFaceArr[heTwin(halfedgeInd)] = INVALID_IND;
      } else {
        // If the twin has already been created, we have an index for the halfedge
        halfedgeInd = heTwinImplicit(twinInd);
      }

      // Hook up a bunch of pointers
      heFaceArr[halfedgeInd] = iFace;
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

// Ensure that each boundary neighborhood is either a disk or a half-disk. Harder to diagnose if we wait until the
// boundary walk below.
#ifndef NGC_SAFETY_CHECKS
  {
    std::vector<char> vertexOnBoundary(nVerticesCount, false);
    for (size_t iHe = 0; iHe < nHalfedgesCount; iHe++) {
      if (heNextArr[iHe] == INVALID_IND) {
        size_t v = heVertexArr[iHe];
        GC_SAFETY_ASSERT(!vertexOnBoundary[v],
                         "vertex " + std::to_string(v) + " appears in more than one boundary loop");
        vertexOnBoundary[v] = true;
      }
    }
  }
#endif

  // == Resolve boundary loops
  nInteriorHalfedgesCount = nHalfedgesCount; // will decrement as we find exterior
  for (size_t iHe = 0; iHe < nHalfedgesCount; iHe++) {

    // If the face pointer is invalid, the halfedge must be along an unresolved boundary loop
    if (heFaceArr[iHe] != INVALID_IND) continue;

    // Create the new boundary loop
    size_t boundaryLoopInd = nFacesCount + nBoundaryLoopsCount;
    fHalfedgeArr.push_back(iHe);
    nBoundaryLoopsCount++;

    // = Walk around the loop (CW)
    size_t currHe = iHe;
    size_t prevHe = INVALID_IND;
    size_t loopCount = 0;
    do {

      // The boundary loop is the face for these halfedges
      heFaceArr[currHe] = boundaryLoopInd;

      // currHe.twin() is a boundary interior halfedge, this is a good time to enforce that v.halfedge() is always the
      // boundary interior halfedge for a boundary vertex.
      size_t currHeT = heTwinImplicit(currHe);
      vHalfedgeArr[heVertexArr[currHeT]] = currHeT;

      // This isn't an interior halfedge.
      nInteriorHalfedgesCount--;

      // Advance to the next halfedge along the boundary
      prevHe = currHe;
      currHe = heTwinImplicit(heNextArr[heTwinImplicit(currHe)]);
      size_t loopCountInnter = 0;
      while (heFaceArr[currHe] != INVALID_IND) {
        if (currHe == iHe) break;
        currHe = heTwinImplicit(heNextArr[currHe]);
        loopCountInnter++;
        GC_SAFETY_ASSERT(loopCountInnter < nHalfedgesCount, "boundary infinite loop orbit");
      }

      // Set the next pointer around the boundary loop
      heNextArr[currHe] = prevHe;

      // Make sure this loop doesn't infinite-loop. Certainly won't happen for proper input, but might happen for
      // bogus input. I don't _think_ it can happen, but there might be some non-manifold input which manfests failure
      // via an infinte loop here, and such a loop is an inconvenient failure mode.
      loopCount++;
      GC_SAFETY_ASSERT(loopCount < nHalfedgesCount, "boundary infinite loop");
    } while (currHe != iHe);
  }

  // SOMEDAY: could shrink_to_fit() std::vectors here, at the cost of a copy. What's preferable?

  // Set capacities and other properties
  nVerticesCapacityCount = nVerticesCount;
  nHalfedgesCapacityCount = nHalfedgesCount;
  nFacesCapacityCount = nFacesCount + nBoundaryLoopsCount;
  nVerticesFillCount = nVerticesCount;
  nHalfedgesFillCount = nHalfedgesCount;
  nFacesFillCount = nFacesCount;
  nBoundaryLoopsFillCount = nBoundaryLoopsCount;
  isCompressedFlag = true;

  // Print some nice statistics
  // printStatistics();
  // std::cout << "Construction took " << pretty_time(FINISH_TIMING(construction)) << std::endl;
}

ManifoldSurfaceMesh::ManifoldSurfaceMesh(const std::vector<size_t>& heNextArr_, const std::vector<size_t>& heVertexArr_,
                                         const std::vector<size_t>& heFaceArr_,
                                         const std::vector<size_t>& vHalfedgeArr_,
                                         const std::vector<size_t>& fHalfedgeArr_, size_t nBoundaryLoopsFillCount_)
    : SurfaceMesh(true) {

  heNextArr = heNextArr_;
  heVertexArr = heVertexArr_;
  heFaceArr = heFaceArr_;
  vHalfedgeArr = vHalfedgeArr_;
  fHalfedgeArr = fHalfedgeArr_;

  // == Set all counts
  nHalfedgesCount = heNextArr.size();
  nEdgesCount = nHalfedgesCount / 2;
  nVerticesCount = vHalfedgeArr.size();
  nFacesCount = fHalfedgeArr.size() - nBoundaryLoopsFillCount_;
  nBoundaryLoopsCount = nBoundaryLoopsFillCount_;
  nVerticesCapacityCount = nVerticesCount;
  nHalfedgesCapacityCount = nHalfedgesCount;
  nEdgesCapacityCount = nHalfedgesCapacityCount / 2;
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
}

ManifoldSurfaceMesh::~ManifoldSurfaceMesh() {}


int ManifoldSurfaceMesh::eulerCharacteristic() const {
  // be sure to do intermediate arithmetic with large *signed* integers
  return static_cast<int>(static_cast<long long int>(nVertices()) - static_cast<long long int>(nEdges()) +
                          static_cast<long long int>(nFaces() + nBoundaryLoops()));
}

int ManifoldSurfaceMesh::genus() const {
  int chi = eulerCharacteristic();
  int boundaryLoops = nBoundaryLoops();
  return (2 - boundaryLoops - chi) / 2;
}


// it is an ERROR for for a manifold surface mesh to ever be nonmanifold/unoriented, so don't even bother checking
// if necessary, one can call validateConnectivity() to check manifoldness manually
bool ManifoldSurfaceMesh::isManifold() { return true; }
bool ManifoldSurfaceMesh::isEdgeManifold() { return true; }
bool ManifoldSurfaceMesh::isOriented() { return true; }


// ==========================================================
// ================        Mutation        ==================
// ==========================================================


Halfedge ManifoldSurfaceMesh::insertVertexAlongEdge(Edge e) {

  // == Gather / create elements
  // Faces are identified as 'A', and 'B'
  bool isBoundary = e.isBoundary();

  // Create new elements
  Vertex newV = getNewVertex();
  Halfedge heANew = getNewEdgeTriple(isBoundary);
  Halfedge heBNew = heANew.twin();
  Edge newE = heANew.edge();

  // Gather old elements
  Halfedge heACenter = e.halfedge();
  Halfedge heBCenter = heACenter.twin();
  Halfedge heANext = heACenter.next();
  Halfedge heBNext = heBCenter.next();
  Halfedge heAPrev = heACenter.prevOrbitFace();
  // Halfedge heBPrev = heBCenter.prevOrbitFace();
  Face fA = heACenter.face();
  Face fB = heBCenter.face();
  Vertex oldVBottom = heACenter.vertex();

  // == Hook up all the pointers

  // New vertex
  vHalfedgeArr[newV.getIndex()] = heACenter.getIndex();

  // New halfedge A
  heNextArr[heANew.getIndex()] = heACenter.getIndex();
  heVertexArr[heANew.getIndex()] = oldVBottom.getIndex();
  heFaceArr[heANew.getIndex()] = fA.getIndex();

  // New halfedge B
  heNextArr[heBNew.getIndex()] = heBNext.getIndex();
  heVertexArr[heBNew.getIndex()] = newV.getIndex();
  heFaceArr[heBNew.getIndex()] = fB.getIndex();

  // Fix pointers for old halfedges
  heNextArr[heBCenter.getIndex()] = heBNew.getIndex();
  heNextArr[heAPrev.getIndex()] = heANew.getIndex();
  heVertexArr[heACenter.getIndex()] = newV.getIndex();

  // Only set this if we broke it, to preseve boundary convention
  if (oldVBottom.halfedge() == heACenter) {
    vHalfedgeArr[oldVBottom.getIndex()] = heANew.getIndex();
  }

  modificationTick++;
  return heACenter;
}


Halfedge ManifoldSurfaceMesh::splitEdgeTriangular(Edge e) {

  // Check triangular assumption
  GC_SAFETY_ASSERT(e.halfedge().face().isTriangle(), "splitEdgeTriangular requires triangular faces");
  GC_SAFETY_ASSERT(e.isBoundary() || e.halfedge().twin().face().isTriangle(),
                   "splitEdgeTriangular requires triangular faces");

  // First operation: insert a new vertex along the edge
  Halfedge he = insertVertexAlongEdge(e);

  { // primary face
    Halfedge heOther = he.next().next();
    connectVertices(he, heOther);
  }

  if (he.twin().isInterior()) { // secondary face
    Halfedge heFirst = he.twin().next();
    Halfedge heOther = heFirst.next().next();
    connectVertices(heFirst, heOther);
  }

  modificationTick++;
  return he;
}


Halfedge ManifoldSurfaceMesh::connectVertices(Halfedge heA, Halfedge heB) {

  // Gather a few values
  Halfedge heAPrev = heA.prevOrbitVertex();
  Halfedge heBPrev = heB.prevOrbitVertex();
  Vertex vA = heA.vertex();
  Vertex vB = heB.vertex();
  Face fA = heA.face();

  // Check some sanity
  GC_SAFETY_ASSERT(heA.face() == heB.face(), "connectVertices(): must lie in same face");
  GC_SAFETY_ASSERT(heA != heBPrev && heAPrev != heB, "connectVertices(): must not be adjacent");
  GC_SAFETY_ASSERT(heA != heB, "connectVertices(): cannot connect vertex to itself inside face");


  // Create new elements
  Halfedge heANew = getNewEdgeTriple(false);
  Halfedge heBNew = heANew.twin();
  Edge eNew = heANew.edge();
  Face fB = getNewFace();


  // == Hook up all the pointers

  // Faces
  fHalfedgeArr[fA.getIndex()] = heANew.getIndex();
  fHalfedgeArr[fB.getIndex()] = heBNew.getIndex();

  // Halfedges
  heNextArr[heANew.getIndex()] = heB.getIndex();
  heVertexArr[heANew.getIndex()] = vA.getIndex();
  heFaceArr[heANew.getIndex()] = fA.getIndex();

  heNextArr[heBNew.getIndex()] = heA.getIndex();
  heVertexArr[heBNew.getIndex()] = vB.getIndex();
  heFaceArr[heBNew.getIndex()] = fB.getIndex();

  heNextArr[heAPrev.getIndex()] = heANew.getIndex();
  heNextArr[heBPrev.getIndex()] = heBNew.getIndex();

  // Set all other new .face pointers to fB
  Halfedge currHe = heA;
  while (currHe != heBNew) {
    heFaceArr[currHe.getIndex()] = fB.getIndex();
    currHe = currHe.next();
  }

  modificationTick++;
  return heANew;
}


std::tuple<Halfedge, Halfedge> ManifoldSurfaceMesh::separateEdge(Edge e) {

  // Must not be a boundary edge
  if (e.isBoundary()) {
    throw std::runtime_error("tried to separate boundary edge");
  }

  // Gather values
  Halfedge he = e.halfedge();
  Vertex vA = he.vertex();
  bool vAIsBoundary = vA.isBoundary();
  Vertex vB = he.twin().vertex();
  bool vBIsBoundary = vB.isBoundary();

  // Swap if needed to simplify case 2, so he.vertex() is always on boundary if any vertex is
  bool swapAB = false; // notice: only possibly swap if case 2 below
  if (vBIsBoundary && !vAIsBoundary) {
    swapAB = true;
    he = he.twin();
    std::swap(vA, vB);
    std::swap(vAIsBoundary, vBIsBoundary);
  }

  // Gather some more values
  Halfedge heT = he.twin();
  Halfedge heTNext = heT.next();
  Halfedge heTPrev = heT.prevOrbitFace();
  Face fA = he.face();
  Face fB = heT.face();

  // Gather boundary loops (set to BoundaryLoop() if they don't exist)
  BoundaryLoop boundaryLoopA = BoundaryLoop();
  if (vAIsBoundary) {
    boundaryLoopA = vA.halfedge().twin().face().asBoundaryLoop();
  }
  BoundaryLoop boundaryLoopB = BoundaryLoop();
  if (vBIsBoundary) {
    boundaryLoopB = vB.halfedge().twin().face().asBoundaryLoop();
  }


  // === Case 1: neither vertex is already boundary
  if (!vAIsBoundary && !vBIsBoundary) {

    // = Create a new (two-sided) boundary loop

    // Get new mesh elements
    Halfedge heN1 = getNewEdgeTriple(true);
    Halfedge heN2 = heN1.twin();
    Edge eN = heN1.edge();
    BoundaryLoop blN = getNewBoundaryLoop();

    // Hook up references
    heNextArr[heT.getIndex()] = heN2.getIndex();
    heNextArr[heN2.getIndex()] = heT.getIndex();
    heNextArr[heN1.getIndex()] = heTNext.getIndex();
    heNextArr[heTPrev.getIndex()] = heN1.getIndex();

    heVertexArr[heN1.getIndex()] = vB.getIndex();
    heVertexArr[heN2.getIndex()] = vA.getIndex();

    heFaceArr[heT.getIndex()] = blN.getIndex();
    heFaceArr[heN1.getIndex()] = fB.getIndex();
    heFaceArr[heN2.getIndex()] = blN.getIndex();

    fHalfedgeArr[fB.getIndex()] = heN1.getIndex();
    fHalfedgeArr[blN.getIndex()] = heT.getIndex();

    vHalfedgeArr[vA.getIndex()] = he.getIndex();
    vHalfedgeArr[vB.getIndex()] = heN1.getIndex();

    modificationTick++;
    return std::tuple<Halfedge, Halfedge>{he, heN1};
  }


  // === Case 2: one vertex is already boundary, other is not
  if (vAIsBoundary && !vBIsBoundary) {

    // Gather some more values
    Halfedge heB = vA.halfedge().twin();
    Halfedge heBN = heB.next();
    BoundaryLoop bl = heB.face().asBoundaryLoop();

    // Create a new vertex, join to the existing boundary loop

    Halfedge heN1 = getNewEdgeTriple(true);
    Halfedge heN2 = heN1.twin();
    Edge eN = heN1.edge();
    Vertex vN = getNewVertex();

    // Hook up references
    heNextArr[heT.getIndex()] = heBN.getIndex();
    heNextArr[heN2.getIndex()] = heT.getIndex();
    heNextArr[heN1.getIndex()] = heTNext.getIndex();
    heNextArr[heTPrev.getIndex()] = heN1.getIndex();
    heNextArr[heB.getIndex()] = heN2.getIndex();

    heVertexArr[heN1.getIndex()] = vB.getIndex();
    heVertexArr[heN2.getIndex()] = vA.getIndex();
    Halfedge heCurr = he;
    do { // set new outgoing halfedge from vN
      heVertexArr[heCurr.getIndex()] = vN.getIndex();
      heCurr = heCurr.next().next().twin();
    } while (heCurr != heBN);
    heVertexArr[heCurr.getIndex()] = vN.getIndex();

    heFaceArr[heT.getIndex()] = bl.asFace().getIndex();
    // std::cout << heFaceArr[heT.getIndex()] << std::endl;
    heFaceArr[heN1.getIndex()] = fB.getIndex();
    heFaceArr[heN2.getIndex()] = bl.asFace().getIndex();

    fHalfedgeArr[fB.getIndex()] = heN1.getIndex();

    vHalfedgeArr[vB.getIndex()] = heN1.getIndex();
    vHalfedgeArr[vN.getIndex()] = he.getIndex();

    ensureEdgeHasInteriorHalfedge(he.edge());

    std::tuple<Halfedge, Halfedge> result{he.edge().halfedge(), heN1};
    if (swapAB) {
      std::swap(std::get<0>(result), std::get<1>(result));
    }
    modificationTick++;
    return result;
  }


  // === Case 3: both vertices are distinct boundaries
  // need to merge boundary loops
  // TODO implement
  if (vAIsBoundary && vBIsBoundary && boundaryLoopA != boundaryLoopB) {
    throw std::runtime_error("not implemented: separateEdge() merging distinct boundaries");
    return std::tuple<Halfedge, Halfedge>{Halfedge(), Halfedge()};
  }


  // === Case 4: both vertices are same boundaries
  // need to split off disconnected compoent of surface
  // TODO implement
  if (vAIsBoundary && vBIsBoundary && boundaryLoopA == boundaryLoopB) {
    throw std::runtime_error("not implemented: separateEdge() creating disconnected components");
    return std::tuple<Halfedge, Halfedge>{Halfedge(), Halfedge()};
  }

  throw std::runtime_error("logically unreachable");
  modificationTick++;
  return std::tuple<Halfedge, Halfedge>{Halfedge(), Halfedge()};
}


Halfedge ManifoldSurfaceMesh::switchHalfedgeSides(Edge e) {

  // NOTE: Written to be safe to call even if the invariant that e.halfedge() is interior is violated, so we can use
  // it to impose that invariant.

  // Gather values
  Halfedge he = e.halfedge();
  Halfedge heN = he.next();
  Halfedge heP = he.prevOrbitVertex();

  Halfedge heT = he.twin();
  Halfedge heTN = heT.next();
  Halfedge heTP = heT.prevOrbitVertex();

  Face fA = he.face();  // might be a boundary loop
  Face fB = heT.face(); // might be a boundary loop

  Vertex vA = he.vertex();
  Vertex vB = heT.vertex();

  // Set references
  heNextArr[he.getIndex()] = heTN.getIndex();
  heNextArr[heTP.getIndex()] = he.getIndex();
  heNextArr[heT.getIndex()] = heN.getIndex();
  heNextArr[heP.getIndex()] = heT.getIndex();

  heFaceArr[he.getIndex()] = fB.getIndex();
  heFaceArr[heT.getIndex()] = fA.getIndex();

  heVertexArr[he.getIndex()] = vB.getIndex();
  heVertexArr[heT.getIndex()] = vA.getIndex();

  fHalfedgeArr[fB.getIndex()] = he.getIndex();
  fHalfedgeArr[fA.getIndex()] = heT.getIndex();

  if (fA.isBoundaryLoop() || vB.halfedge() == heT) {
    vHalfedgeArr[vB.getIndex()] = he.getIndex();
  }
  if (fB.isBoundaryLoop() || vA.halfedge() == he) {
    vHalfedgeArr[vA.getIndex()] = heT.getIndex();
  }

  modificationTick++;
  return e.halfedge();
}

bool ManifoldSurfaceMesh::ensureEdgeHasInteriorHalfedge(Edge e) {
  if (!e.halfedge().isInterior()) {
    switchHalfedgeSides(e);
    modificationTick++;
    return true;
  }
  return false;
}

/*

Halfedge ManifoldSurfaceMesh::connectVertices(Face faceIn, Vertex vAIn, Vertex vBIn) {

  // == Find useful halfedges around the face
  Halfedge heANext;
  Halfedge heBNext;
  Halfedge heAPrev;
  Halfedge heBPrev;
  for (Halfedge he : face.adjacentHalfedges()) {
    if (he.vertex() == vA) {
      heANext = he.ptr;
    }
    if (he.vertex() == vB) {
      heBNext = he.ptr;
    }
  }
  for (Halfedge he : face.adjacentHalfedges()) {
    if (he.next().ptr == heANext) {
      heAPrev = he.ptr;
    }
    if (he.next().ptr == heBNext) {
      heBPrev = he.ptr;
    }
  }

  // == Detect bad cases
  if (vA == vB) throw std::logic_error("Tried to connect vertex to self");
  if (heANext == heBPrev || heBNext == heBPrev) throw std::logic_error("Tried to connect adjacent vertices");

  // == Gather other elements
  Face* fA = heBNext->face;
  Vertex* vAp = vA.ptr;
  Vertex* vBp = vB.ptr;

  // == Create new elements
  Halfedge* heANew = getNewHalfedge(true);
  Halfedge* heBNew = getNewHalfedge(true);
  Edge* eNew = getNewEdge();
  Face* fB = getNewFace();


  // == Hook up all the pointers

  // Faces
  fA->halfedge = heBNext;
  fB->halfedge = heANext;

  // Vertices
  vAp->halfedge = heANew;
  vBp->halfedge = heBNew;

  // New edge
  eNew->halfedge = heANew;

  // Halfedges
  heANew->twin = heBNew;
  heANew->next = heBNext;
  heANew->vertex = vAp;
  heANew->edge = eNew;
  heANew->face = fA;

  heBNew->twin = heANew;
  heBNew->next = heANext;
  heBNew->vertex = vBp;
  heBNew->edge = eNew;
  heBNew->face = fB;

  heAPrev->next = heANew;
  heBPrev->next = heBNew;

  // Set all other new .face pointers to fB
  Halfedge* currHe = heANext;
  while (currHe != heBNew) {
    currHe->face = fB;
    currHe = currHe->next;
  }

  return heANew;
}

*/

Vertex ManifoldSurfaceMesh::insertVertex(Face fIn) {

  // Create the new center vertex
  Vertex centerVert = getNewVertex();

  // Count degree to allocate elements
  size_t faceDegree = fIn.degree();

  // == Create new halfedges/edges/faces around the center vertex

  // Create all of the new elements first, then hook them up below
  std::vector<Face> innerFaces;
  std::vector<Halfedge> leadingHalfedges(faceDegree); // the one that points towards the center
  std::vector<Halfedge> trailingHalfedges(faceDegree);
  std::vector<Edge> innerEdges(faceDegree); // aligned with leading he
  for (size_t i = 0; i < faceDegree; i++) {
    // Re-use first face
    if (i == 0) {
      innerFaces.push_back(fIn);
    } else {
      innerFaces.push_back(getNewFace());
    }

    // Get the new edge group
    Halfedge newHe = getNewEdgeTriple(false);

    leadingHalfedges[i] = newHe;
    trailingHalfedges[(i + 1) % faceDegree] = newHe.twin();
    innerEdges[i] = newHe.edge();
  }

  // Form this list before we start, because we're about to start breaking pointers
  std::vector<Halfedge> faceBoundaryHalfedges;
  for (Halfedge he : fIn.adjacentHalfedges()) {
    faceBoundaryHalfedges.push_back(he);
  }

  // Connect up all the pointers
  // Each iteration processes one inner face
  for (size_t i = 0; i < faceDegree; i++) {

    // Gather pointers
    Face f = innerFaces[i];
    Edge e = innerEdges[i];
    Edge prevE = innerEdges[(i + faceDegree - 1) % faceDegree];
    Halfedge leadingHe = leadingHalfedges[i];
    Halfedge trailingHe = trailingHalfedges[i];
    Halfedge boundaryHe = faceBoundaryHalfedges[i];
    Halfedge nextTrailingHe = trailingHalfedges[(i + 1) % faceDegree];
    Halfedge prevLeadingHe = leadingHalfedges[(i + faceDegree - 1) % faceDegree];

    // face
    fHalfedgeArr[f.getIndex()] = boundaryHe.getIndex();

    // leading halfedge
    heNextArr[leadingHe.getIndex()] = trailingHe.getIndex();
    heVertexArr[leadingHe.getIndex()] = boundaryHe.next().vertex().getIndex();
    heFaceArr[leadingHe.getIndex()] = f.getIndex();

    // trailing halfedge
    heNextArr[trailingHe.getIndex()] = boundaryHe.getIndex();
    heVertexArr[trailingHe.getIndex()] = centerVert.getIndex();
    heFaceArr[trailingHe.getIndex()] = f.getIndex();

    // boundary halfedge
    heNextArr[boundaryHe.getIndex()] = leadingHe.getIndex();
    heFaceArr[boundaryHe.getIndex()] = f.getIndex();
  }

  vHalfedgeArr[centerVert.getIndex()] = trailingHalfedges[0].getIndex();

  modificationTick++;
  return centerVert;
}


Vertex ManifoldSurfaceMesh::collapseEdgeTriangular(Edge e) {
  /*  must maintain these
      std::vector<size_t> heNextArr;    // he.next(), forms a circular singly-linked list in each face
      std::vector<size_t> heVertexArr;  // he.vertex()
      std::vector<size_t> heFaceArr;    // he.face()
      std::vector<size_t> vHalfedgeArr; // v.halfedge()
      std::vector<size_t> fHalfedgeArr; // f.halfedge()
  */

  // check triangular
  GC_SAFETY_ASSERT(e.halfedge().face().isTriangle(), "neighborhood must be triangular");
  GC_SAFETY_ASSERT(e.isBoundary() || e.halfedge().twin().face().isTriangle(), "neighborhood must be triangular");

  // assuming on triangle mesh
  if (e.isBoundary()) {
    // Gather some values
    Halfedge heA0 = e.halfedge();
    if (heA0.vertex().degree() == 2) {
      heA0 = heA0.next().next();
      e = heA0.edge();
    }

    // check if edge is part of a 'pinch' triangle
    if (heA0.twin().next().next().next() == heA0.twin()) {
      return Vertex();
    }
    for (Halfedge he1 : heA0.tipVertex().outgoingHalfedges()) {
      for (Halfedge he2 : he1.tipVertex().outgoingHalfedges()) {
        if (!(heA0.next() == he1 && he1.next() == he2) &&
            !(he1.twin().next().twin() == heA0 && he2.twin().next().twin() == he1)) { // not just going around a face
          if (he2.tipVertex() == heA0.vertex()) {
            return Vertex();
          }
        }
      }
    }

    Halfedge heA1 = heA0.next();
    Halfedge heA2 = heA1.next();
    Halfedge heC2 = heA2.twin();
    Halfedge heC0 = heC2.next();
    Halfedge heC1 = heC0.next();

    Face fA = heA0.face();
    Face fC = heC0.face();

    Vertex vA = heA0.vertex();
    Vertex vB = heA1.vertex();
    Vertex vC = heC0.vertex();

    // Special boundary values
    Halfedge heB0 = heA0.twin();
    Halfedge heB0prev;
    for (Halfedge he : vB.incomingHalfedges()) {
      if (!he.isInterior()) {
        heB0prev = he;
        break;
      }
    }
    Halfedge heB0next = heB0.next();
    Face fB = heB0.face(); // boundary loop

    // Reassign connections
    std::vector<Halfedge> toReassignVertex;
    for (Halfedge he : vA.outgoingHalfedges()) {
      toReassignVertex.push_back(he);
    }
    for (Halfedge he : toReassignVertex) {
      heVertexArr[he.getIndex()] = vB.getIndex();
    }
    heNextArr[heC1.getIndex()] = heA1.getIndex();
    heNextArr[heA1.getIndex()] = heC0.getIndex();
    heNextArr[heB0prev.getIndex()] = heB0next.getIndex();
    // do not need to do vB.halfedge
    heFaceArr[heA1.getIndex()] = fC.getIndex();
    if (!vC.isBoundary()) {
      vHalfedgeArr[vC.getIndex()] = heC0.getIndex();
    }
    fHalfedgeArr[fC.getIndex()] = heC0.getIndex();
    fHalfedgeArr[fB.getIndex()] = heB0next.getIndex();

    // Actually delete
    deleteEdgeBundle(e);
    deleteEdgeBundle(heA2.edge());
    deleteElement(vA);
    deleteElement(fA);

    return vB;
  }

  else {
    Halfedge heA0 = e.halfedge();
    if (heA0.vertex().isBoundary() && heA0.twin().vertex().isBoundary()) {
      return Vertex();
    }
    if (heA0.vertex().isBoundary()) {
      heA0 = heA0.twin();
    }

    // check if edge is part of a 'pinch' triangle
    for (Halfedge he1 : heA0.tipVertex().outgoingHalfedges()) {
      for (Halfedge he2 : he1.tipVertex().outgoingHalfedges()) {
        if (!(heA0.next() == he1 && he1.next() == he2) &&
            !(he1.twin().next().twin() == heA0 && he2.twin().next().twin() == he1)) { // not just going around a face
          if (he2.tipVertex() == heA0.vertex()) {
            // std::cerr<<heA0<<" "<<he1<<" "<<he2<<"\n";
            return Vertex();
          }
        }
      }
    }

    if (heA0.vertex().degree() > 3) {
      Halfedge heA1 = heA0.next();
      Halfedge heA2 = heA1.next();
      Halfedge heB0 = heA0.twin();
      Halfedge heB1 = heB0.next();
      Halfedge heB2 = heB1.next();
      Halfedge heC2 = heA2.twin();
      Halfedge heC0 = heC2.next();
      Halfedge heC1 = heC0.next();
      Halfedge heD1 = heB1.twin();
      Halfedge heD2 = heD1.next();
      Halfedge heD0 = heD2.next();

      Face fA = heA0.face();
      Face fB = heB0.face();
      Face fC = heC0.face();
      Face fD = heD0.face();

      Vertex vA = heA0.vertex();
      Vertex vB = heB0.vertex();
      Vertex vC = heC0.vertex();
      Vertex vD = heD1.vertex();

      // Reassign connections
      std::vector<Halfedge> toReassignVertex;
      for (Halfedge he : vA.outgoingHalfedges()) {
        toReassignVertex.push_back(he);
      }
      for (Halfedge he : toReassignVertex) {
        heVertexArr[he.getIndex()] = vB.getIndex();
      }

      heNextArr[heD0.getIndex()] = heB2.getIndex();
      heNextArr[heB2.getIndex()] = heD2.getIndex();
      heNextArr[heC1.getIndex()] = heA1.getIndex();
      heNextArr[heA1.getIndex()] = heC0.getIndex();
      heFaceArr[heB2.getIndex()] = fD.getIndex();
      heFaceArr[heA1.getIndex()] = fC.getIndex();
      fHalfedgeArr[fC.getIndex()] = heC0.getIndex();
      fHalfedgeArr[fD.getIndex()] = heD0.getIndex();
      // boundary vertex already safe with boundary edge
      if (!vB.isBoundary()) {
        vHalfedgeArr[vB.getIndex()] = heA1.getIndex();
      }
      if (!vC.isBoundary()) {
        vHalfedgeArr[vC.getIndex()] = heC0.getIndex();
      }
      if (!vD.isBoundary()) {
        vHalfedgeArr[vD.getIndex()] = heB2.getIndex();
      }

      // Actually delete
      deleteEdgeBundle(e);
      deleteEdgeBundle(heB1.edge());
      deleteEdgeBundle(heA2.edge());
      deleteElement(vA);
      deleteElement(fA);
      deleteElement(fB);

      return vB;
    }

    else if (heA0.vertex().degree() == 3) {
      Halfedge heA1 = heA0.next();
      Halfedge heA2 = heA1.next();
      Halfedge heB0 = heA0.twin();
      Halfedge heB1 = heB0.next();
      Halfedge heB2 = heB1.next();
      Halfedge heC2 = heA2.twin();
      Halfedge heC0 = heC2.next();
      Halfedge heC1 = heC0.next();

      Face fA = heA0.face();
      Face fB = heB0.face();
      Face fC = heC0.face();

      Vertex vA = heA0.vertex();
      Vertex vB = heB0.vertex();
      Vertex vC = heC0.vertex();
      Vertex vD = heC1.vertex();

      heNextArr[heA1.getIndex()] = heC0.getIndex();
      heNextArr[heC0.getIndex()] = heB2.getIndex();
      heNextArr[heB2.getIndex()] = heA1.getIndex();
      heFaceArr[heB2.getIndex()] = fC.getIndex();
      heFaceArr[heA1.getIndex()] = fC.getIndex();
      fHalfedgeArr[fC.getIndex()] = heC0.getIndex();
      // boundary vertex already safe with boundary edge
      if (!vB.isBoundary()) {
        vHalfedgeArr[vB.getIndex()] = heA1.getIndex();
      }
      if (!vC.isBoundary()) {
        vHalfedgeArr[vC.getIndex()] = heC0.getIndex();
      }
      if (!vD.isBoundary()) {
        vHalfedgeArr[vD.getIndex()] = heB2.getIndex();
      }

      // Actually delete
      deleteEdgeBundle(e);
      deleteEdgeBundle(heB1.edge());
      deleteEdgeBundle(heA2.edge());
      deleteElement(vA);
      deleteElement(fA);
      deleteElement(fB);

      return vB;
    }

    else {
      // return heA0.vertex();
      throw std::runtime_error("Collapsing interior vertices with degree < 3 is not supported.");
    }
  }

  modificationTick++;
}



bool ManifoldSurfaceMesh::removeFaceAlongBoundary(Face f) {

  // Find the boundary halfedge
  Halfedge heB;
  int bCount = 0;
  int fCount = 0;
  for (Halfedge he : f.adjacentHalfedges()) {
    if (!he.twin().isInterior()) {
      bCount++;
      heB = he;
    }
    fCount++;
  }
  if (bCount == 0) {
    throw std::runtime_error("called on non-boundary face");
  }
  if (bCount == 1) {
    // Remove a non-ear boundary face with one boundary edge


    // Gather values
    Halfedge heBNext = heB.next();
    Halfedge heBPrev = heB.prevOrbitFace();

    Halfedge heT = heB.twin();
    Halfedge heTNext = heT.next();
    Halfedge heTPrev = heT.prevOrbitVertex();

    Face bLoop = heT.face();


    // Opposite vertex must not be a bounary vertex or this creates a nonmanifold mesh (imagine hourglass)
    if (heBPrev.vertex().isBoundary()) {
      return false;
    }

    // Update refs
    for (Halfedge he : f.adjacentHalfedges()) {
      heFaceArr[he.getIndex()] = bLoop.getIndex();
    }

    // Next refs
    heNextArr[heBPrev.getIndex()] = heTNext.getIndex();
    heNextArr[heTPrev.getIndex()] = heBNext.getIndex();

    // Vertex halfedges
    vHalfedgeArr[heTNext.vertex().getIndex()] = heBPrev.twin().getIndex();
    ensureVertexHasBoundaryHalfedge(heBPrev.vertex());

    fHalfedgeArr[bLoop.getIndex()] = heTNext.getIndex();

    Halfedge currHe = heBNext;
    do {
      Halfedge nextHe = currHe.next();
      ensureEdgeHasInteriorHalfedge(currHe.edge());
      currHe = nextHe;
    } while (currHe != heTNext);

    deleteElement(f);
    deleteEdgeBundle(heB.edge());
    modificationTick++;
    return true;

    /*
    Halfedge* he0 = heBoundary.ptr;
    Halfedge* he0T = he0->twin;
    Halfedge* he1 = he0->next;
    Halfedge* he2 = he1->next;
    Vertex* v0 = he0->vertex;
    Vertex* v1 = he1->vertex;
    Vertex* v2 = he2->vertex;
    Face* fRemove = he0->face;
    Face* bLoop = he0T->face;

    // Vertex halfedges
    v0->halfedge = he2->twin;
    v2->halfedge = he1->twin;

    // Nexts
    he2->next = he0T->next;
    v1->halfedge->twin->next = he1;

    // Faces
    he1->face = bLoop;
    he2->face = bLoop;

    // mark boundary
    v2->isBoundary = true;
    he1->isReal = false;
    he2->isReal = false;

    deleteElement(he0->edge);
    deleteElement(he0);
    deleteElement(he0T);
    deleteElement(fRemove);

    isCanonicalFlag = false;
  modificationTick++;
    return true;
    */

  } else if (bCount == 2) {
    // Remove an "ear" along the boundary

    /*
    // Gather elements
    Halfedge* he0 = f.halfedge().ptr;
    while (!he0->twin->isReal) he0 = he0->next;
    Halfedge* he0T = he0->twin;
    Halfedge* he1 = he0->next;
    Halfedge* he1T = he1->twin;
    Edge* e1 = he1->edge;
    Halfedge* he2 = he1->next;
    Halfedge* he2T = he2->twin;
    Edge* e2 = he2->edge;
    Vertex* v0 = he0->vertex;
    Vertex* v1 = he1->vertex;
    Vertex* v2 = he2->vertex;
    Face* fRemove = he0->face;

    Halfedge* heNextArr = he1T->next;
    Halfedge* hePrev = he0T;
    while (hePrev->isReal) hePrev = hePrev->next->twin;

    // Vertex halfedges
    v0->halfedge = hePrev->twin;
    v1->halfedge = he0T;

    // Nexts
    hePrev->next = heNextArr;

    // Boundary loop
    hePrev->face->halfedge = hePrev;

    // mark boundary
    he0->isReal = false;

    deleteElement(fRemove);
    deleteElement(v2);
    deleteElement(he1);
    deleteElement(he1T);
    deleteElement(e1);
    deleteElement(he2);
    deleteElement(he2T);
    deleteElement(e2);

    isCanonicalFlag = false;
  modificationTick++;
    return true;
    */

    // Not supported yet
    return false;

  } else {
    // Remove entire component

    /*
    Halfedge* he0 = heBoundary.ptr;
    Halfedge* he0T = he0->twin;
    Edge* e0 = he0->edge;
    Halfedge* he1 = he0->next;
    Halfedge* he1T = he1->twin;
    Edge* e1 = he1->edge;
    Halfedge* he2 = he1->next;
    Halfedge* he2T = he2->twin;
    Edge* e2 = he2->edge;
    Vertex* v0 = he0->vertex;
    Vertex* v1 = he1->vertex;
    Vertex* v2 = he2->vertex;
    Face* fFace = he0->face;
    Face* fBound = he0T->face;


    deleteElement(he0);
    deleteElement(he1);
    deleteElement(he2);

    deleteElement(he0T);
    deleteElement(he1T);
    deleteElement(he2T);

    deleteElement(e0);
    deleteElement(e1);
    deleteElement(e2);

    deleteElement(v0);
    deleteElement(v1);
    deleteElement(v2);

    deleteElement(fFace);
    deleteElement(fBound);

    isCanonicalFlag = false;
  modificationTick++;
    return true;
    */

    // The removal/insertion code doesn't support changing boundary structure yet
    return false;
  }
}

Face ManifoldSurfaceMesh::removeVertex(Vertex v) {
  if (v.isBoundary()) {
    throw std::runtime_error("not implemented");
  }

  // Halfedges/edges/faces that will be removed
  // (except first face)
  std::vector<Halfedge> toRemove;
  std::vector<Halfedge> ringHalfedges;
  for (Halfedge he : v.outgoingHalfedges()) {
    toRemove.push_back(he);

    // The one-ring must not contain any other copies of v, or we cannot remove the vertex
    Halfedge oppHe = he.next();
    if (oppHe.vertex() == v || oppHe.twin().vertex() == v) {
      return Face();
    }
    ringHalfedges.push_back(oppHe);
  }

  Face keepFace = toRemove[0].face();

  // Hook up next and face refs for the halfedges along the ring
  size_t N = ringHalfedges.size();
  for (size_t i = 0; i < N; i++) {
    heNextArr[ringHalfedges[(i + 1) % N].getIndex()] = ringHalfedges[i].getIndex(); // since outgoingHalfedges orbits CW
    heFaceArr[ringHalfedges[i].getIndex()] = keepFace.getIndex();

    if (toRemove[i].twin().vertex().halfedge().twin() == toRemove[i]) {
      // only update vHalfedgeArr if needed to avoid disturbing boundary halfedges
      vHalfedgeArr[toRemove[i].twin().vertex().getIndex()] = ringHalfedges[i].getIndex();
    }
  }
  fHalfedgeArr[keepFace.getIndex()] = ringHalfedges[0].getIndex();

  // Actually delete all of the elements
  for (Halfedge he : toRemove) {
    // delete the edge before the face since deleteEdgeBundle() needs to check manifold-ness (see note there)
    Face f = he.face();
    deleteEdgeBundle(he.edge());
    if (f != keepFace) {
      deleteElement(f);
    }
  }
  deleteElement(v);

  modificationTick++;
  return keepFace;
}


void ManifoldSurfaceMesh::ensureVertexHasBoundaryHalfedge(Vertex v) {
  while (true) {
    Halfedge heT = v.halfedge().twin();
    if (!heT.isInterior()) {
      break;
    }
    vHalfedgeArr[v.getIndex()] = heT.next().getIndex();
  }
  modificationTick++;
}

/*

bool ManifoldSurfaceMesh::removeFaceAlongBoundary(Face f) {

  // Find the boundary halfedge
  Halfedge heBoundary;
  int bCount = 0;
  for (Halfedge he : f.adjacentHalfedges()) {
    if (!he.twin().isReal()) {
      bCount++;
      heBoundary = he;
    }
  }
  if (bCount == 0) {
    throw std::runtime_error("called on non-boundary face");
  }
  if (bCount == 1) {
    // Remove a non-ear boundary face with one boundary edge

    Halfedge* he0 = heBoundary.ptr;
    Halfedge* he0T = he0->twin;
    Halfedge* he1 = he0->next;
    Halfedge* he2 = he1->next;
    Vertex* v0 = he0->vertex;
    Vertex* v1 = he1->vertex;
    Vertex* v2 = he2->vertex;
    Face* fRemove = he0->face;
    Face* bLoop = he0T->face;

    // Vertex halfedges
    v0->halfedge = he2->twin;
    v2->halfedge = he1->twin;

    // Nexts
    he2->next = he0T->next;
    v1->halfedge->twin->next = he1;

    // Faces
    he1->face = bLoop;
    he2->face = bLoop;

    // mark boundary
    v2->isBoundary = true;
    he1->isReal = false;
    he2->isReal = false;

    deleteElement(he0->edge);
    deleteElement(he0);
    deleteElement(he0T);
    deleteElement(fRemove);

    return true;

  } else if (bCount == 2) {
    // Remove an "ear" along the boundary

    // Gather elements
    Halfedge* he0 = f.halfedge().ptr;
    while (!he0->twin->isReal) he0 = he0->next;
    Halfedge* he0T = he0->twin;
    Halfedge* he1 = he0->next;
    Halfedge* he1T = he1->twin;
    Edge* e1 = he1->edge;
    Halfedge* he2 = he1->next;
    Halfedge* he2T = he2->twin;
    Edge* e2 = he2->edge;
    Vertex* v0 = he0->vertex;
    Vertex* v1 = he1->vertex;
    Vertex* v2 = he2->vertex;
    Face* fRemove = he0->face;

    Halfedge* heNextArr = he1T->next;
    Halfedge* hePrev = he0T;
    while (hePrev->isReal) hePrev = hePrev->next->twin;

    // Vertex halfedges
    v0->halfedge = hePrev->twin;
    v1->halfedge = he0T;

    // Nexts
    hePrev->next = heNextArr;

    // Boundary loop
    hePrev->face->halfedge = hePrev;

    // mark boundary
    he0->isReal = false;

    deleteElement(fRemove);
    deleteElement(v2);
    deleteElement(he1);
    deleteElement(he1T);
    deleteElement(e1);
    deleteElement(he2);
    deleteElement(he2T);
    deleteElement(e2);

    return true;

  } else {
    // Remove entire component

    // Halfedge* he0 = heBoundary.ptr;
    // Halfedge* he0T = he0->twin;
    // Edge* e0 = he0->edge;
    // Halfedge* he1 = he0->next;
    // Halfedge* he1T = he1->twin;
    // Edge* e1 = he1->edge;
    // Halfedge* he2 = he1->next;
    // Halfedge* he2T = he2->twin;
    // Edge* e2 = he2->edge;
    // Vertex* v0 = he0->vertex;
    // Vertex* v1 = he1->vertex;
    // Vertex* v2 = he2->vertex;
    // Face* fFace = he0->face;
    // Face* fBound = he0T->face;


    // deleteElement(he0);
    // deleteElement(he1);
    // deleteElement(he2);

    // deleteElement(he0T);
    // deleteElement(he1T);
    // deleteElement(he2T);

    // deleteElement(e0);
    // deleteElement(e1);
    // deleteElement(e2);

    // deleteElement(v0);
    // deleteElement(v1);
    // deleteElement(v2);

    // deleteElement(fFace);
    // deleteElement(fBound);

    // return true;

    // The removal/insertion code doesn't support changing boundary structure yet
    return false;
  }
}

*/

std::vector<Face> ManifoldSurfaceMesh::triangulate(Face f) {
  GC_SAFETY_ASSERT(!f.isBoundaryLoop(), "cannot triangulate boundary loop");

  if (f.isTriangle()) {
    return {f};
  }


  std::vector<Halfedge> neighHalfedges;
  for (Halfedge he : f.adjacentHalfedges()) {
    neighHalfedges.emplace_back(he);
  }

  std::vector<Face> allFaces;
  allFaces.emplace_back(f);

  // currently doing a fan triangulation. chould do something better.
  Halfedge connectHe = f.halfedge();
  for (size_t i = 2; i + 1 < neighHalfedges.size(); i++) {
    connectHe = connectVertices(connectHe, neighHalfedges[i]);
    allFaces.emplace_back(connectHe.twin().face());
  }

  modificationTick++;
  return allFaces;
}


// All are automatically true on a manifold mesh
void ManifoldSurfaceMesh::separateNonmanifoldEdges() {}
VertexData<Vertex> ManifoldSurfaceMesh::separateNonmanifoldVertices() {
  VertexData<Vertex> parents(*this);
  for (Vertex v : vertices()) parents[v] = v;
  return parents;
}
void ManifoldSurfaceMesh::greedilyOrientFaces() {}

bool ManifoldSurfaceMesh::hasBoundary() { return nBoundaryLoopsCount > 0; }

std::unique_ptr<ManifoldSurfaceMesh> ManifoldSurfaceMesh::copy() const {
  ManifoldSurfaceMesh* newMesh = new ManifoldSurfaceMesh();
  copyInternalFields(*newMesh);
  return std::unique_ptr<ManifoldSurfaceMesh>(newMesh);
}

std::unique_ptr<SurfaceMesh> ManifoldSurfaceMesh::copyToSurfaceMesh() const {
  std::unique_ptr<ManifoldSurfaceMesh> cMesh = copy();
  return std::unique_ptr<SurfaceMesh>(cMesh.release());
}

} // namespace surface
} // namespace geometrycentral
