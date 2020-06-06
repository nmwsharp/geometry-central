#include "geometrycentral/surface/halfedge_mesh.h"

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

HalfedgeMesh::HalfedgeMesh() {}

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

HalfedgeMesh::HalfedgeMesh(const std::vector<std::vector<size_t>>& polygons, bool requireManifold) {
  if (requireManifold) {
    // constructHalfedgeMeshManifold(polygons);
    constructHalfedgeMeshNonmanifold(polygons, {});
  } else {
    constructHalfedgeMeshNonmanifold(polygons, {});
  }
}

void HalfedgeMesh::constructHalfedgeMeshManifold(const std::vector<std::vector<size_t>>& polygons) {
  useImplicitTwinFlag = true;

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

  // Ensure that each boundary neighborhood is either a disk or a half-disk. Harder to diagnose if we wait until the
  // boundary walk below.
#ifndef NGC_SAFTEY_CHECKS
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

#ifndef NGC_SAFTEY_CHECKS
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

void HalfedgeMesh::constructHalfedgeMeshNonmanifold(const std::vector<std::vector<size_t>>& polygons,
                                                    const std::vector<std::vector<std::tuple<size_t, size_t>>>& twins) {

  useImplicitTwinFlag = false;

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

HalfedgeMesh::HalfedgeMesh(const std::vector<std::vector<size_t>>& polygons,
                           const std::vector<std::vector<std::tuple<size_t, size_t>>>& twins,
                           bool allowVertexNonmanifold) {

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
#ifndef NGC_SAFTEY_CHECKS
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

#ifndef NGC_SAFTEY_CHECKS
  if (!allowVertexNonmanifold) { // Check that the input was manifold in the sense that each vertex has a single
                                 // connected loop of faces around it.
    std::vector<char> halfedgeSeen(nHalfedgesCount, false);
    for (size_t iV = 0; iV < nVerticesCount; iV++) {

      // For each vertex, orbit around the outgoing halfedges. This _should_ touch every halfedge.
      size_t currHe = vHalfedgeArr[iV];
      size_t firstHe = currHe;
      do {

        GC_SAFETY_ASSERT(!halfedgeSeen[currHe], "somehow encountered outgoing halfedge before orbiting v");
        halfedgeSeen[currHe] = true;

        currHe = heNextArr[heTwin(currHe)];
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


HalfedgeMesh::~HalfedgeMesh() {
  for (auto& f : meshDeleteCallbackList) {
    f();
  }
}


// ==========================================================
// ================       Utilities        ==================
// ==========================================================

void HalfedgeMesh::printStatistics() const {
  std::cout << "Halfedge mesh with: " << std::endl;
  std::cout << "    # verts =  " << nVertices() << std::endl;
  std::cout << "    # edges =  " << nEdges() << std::endl;
  std::cout << "    # faces =  " << nFaces() << std::endl;
  std::cout << "    # halfedges =  " << nHalfedges() << "  (" << nInteriorHalfedges() << " interior, "
            << nExteriorHalfedges() << " exterior)" << std::endl;
  std::cout << "      and " << nBoundaryLoops() << " boundary components. " << std::endl;
}

bool HalfedgeMesh::isTriangular() {
  for (Face f : faces()) {
    if (!f.isTriangle()) {
      return false;
    }
  }
  return true;
}

size_t HalfedgeMesh::nConnectedComponents() {
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

size_t HalfedgeMesh::nInteriorVertices() {
  size_t nInteriorVertices = 0;
  for (const Vertex v : vertices()) {
    if (!v.isBoundary()) {
      nInteriorVertices++;
    }
  }
  return nInteriorVertices;
}


int HalfedgeMesh::eulerCharacteristic() const {
  // be sure to do intermediate arithmetic with large *signed* integers
  return static_cast<int>(static_cast<long long int>(nVertices()) - static_cast<long long int>(nEdges()) +
                          static_cast<long long int>(nFaces() + nBoundaryLoops()));
}

int HalfedgeMesh::genus() const {
  int chi = eulerCharacteristic();
  int boundaryLoops = nBoundaryLoops();
  return (2 - boundaryLoops - chi) / 2;
}

VertexData<size_t> HalfedgeMesh::getVertexIndices() {
  VertexData<size_t> indices(*this);
  size_t i = 0;
  for (Vertex v : vertices()) {
    indices[v] = i;
    i++;
  }
  return indices;
}

VertexData<size_t> HalfedgeMesh::getInteriorVertexIndices() {
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

HalfedgeData<size_t> HalfedgeMesh::getHalfedgeIndices() {
  HalfedgeData<size_t> indices(*this);
  size_t i = 0;
  for (Halfedge he : halfedges()) {
    indices[he] = i;
    i++;
  }
  return indices;
}

CornerData<size_t> HalfedgeMesh::getCornerIndices() {
  CornerData<size_t> indices(*this);
  size_t i = 0;
  for (Corner c : corners()) {
    indices[c] = i;
    i++;
  }
  return indices;
}

EdgeData<size_t> HalfedgeMesh::getEdgeIndices() {
  EdgeData<size_t> indices(*this);
  size_t i = 0;
  for (Edge e : edges()) {
    indices[e] = i;
    i++;
  }
  return indices;
}


FaceData<size_t> HalfedgeMesh::getFaceIndices() {
  FaceData<size_t> indices(*this);
  size_t i = 0;
  for (Face f : faces()) {
    indices[f] = i;
    i++;
  }
  return indices;
}

BoundaryLoopData<size_t> HalfedgeMesh::getBoundaryLoopIndices() {
  BoundaryLoopData<size_t> indices(*this);
  size_t i = 0;
  for (BoundaryLoop bl : boundaryLoops()) {
    indices[bl] = i;
    i++;
  }
  return indices;
}


std::unique_ptr<HalfedgeMesh> HalfedgeMesh::copy() const {
  HalfedgeMesh* newMesh = new HalfedgeMesh();


  // == Copy _all_ the fields!

  // Raw data buffers (underlying std::vectors duplicate storage automatically)
  newMesh->heNextArr = heNextArr;
  newMesh->heVertexArr = heVertexArr;
  newMesh->heFaceArr = heFaceArr;
  newMesh->vHalfedgeArr = vHalfedgeArr;
  newMesh->fHalfedgeArr = fHalfedgeArr;
  newMesh->heSiblingArr = heSiblingArr;
  newMesh->heEdgeArr = heEdgeArr;
  newMesh->eHalfedgeArr = eHalfedgeArr;

  // counts and flags
  newMesh->nHalfedgesCount = nHalfedgesCount;
  newMesh->nInteriorHalfedgesCount = nInteriorHalfedgesCount;
  newMesh->nEdgesCount = nEdgesCount;
  newMesh->nVerticesCount = nVerticesCount;
  newMesh->nFacesCount = nFacesCount;
  newMesh->nBoundaryLoopsCount = nBoundaryLoopsCount;
  newMesh->nVerticesCapacityCount = nVerticesCapacityCount;
  newMesh->nHalfedgesCapacityCount = nHalfedgesCapacityCount;
  newMesh->nEdgesCapacityCount = nEdgesCapacityCount;
  newMesh->nFacesCapacityCount = nFacesCapacityCount;
  newMesh->nVerticesFillCount = nVerticesFillCount;
  newMesh->nHalfedgesFillCount = nHalfedgesFillCount;
  newMesh->nEdgesFillCount = nEdgesFillCount;
  newMesh->nFacesFillCount = nFacesFillCount;
  newMesh->nBoundaryLoopsFillCount = nBoundaryLoopsFillCount;

  newMesh->isCompressedFlag = isCompressedFlag;
  newMesh->useImplicitTwinFlag = useImplicitTwinFlag;


  // Note: _don't_ copy callbacks lists! New mesh has new callbacks

  return std::unique_ptr<HalfedgeMesh>(newMesh);
}

std::vector<std::vector<size_t>> HalfedgeMesh::getFaceVertexList() {

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


// ==========================================================
// ================        Mutation        ==================
// ==========================================================


bool HalfedgeMesh::flip(Edge eFlip) {
  if (eFlip.isBoundary()) return false;
  if (!eFlip.isManifold()) return false;

  // Get halfedges of first face
  Halfedge ha1 = eFlip.halfedge();
  Halfedge ha2 = ha1.next();
  Halfedge ha3 = ha2.next();
  if (ha3.next() != ha1) return false; // not a triangle

  // Get halfedges of second face
  Halfedge hb1 = ha1.twin();
  Halfedge hb2 = hb1.next();
  Halfedge hb3 = hb2.next();
  if (hb3.next() != hb1) return false; // not a triangle

  if (ha2 == hb1 || hb2 == ha1) return false; // incident on degree 1 vertex

  // Get vertices and faces
  Vertex va = ha1.vertex();
  Vertex vb = hb1.vertex();
  Vertex vc = ha3.vertex();
  Vertex vd = hb3.vertex();
  Face fa = ha1.face();
  Face fb = hb1.face();

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

  modificationTick++;
  return true;
}


Halfedge HalfedgeMesh::insertVertexAlongEdge(Edge e) {

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


Halfedge HalfedgeMesh::splitEdgeTriangular(Edge e) {

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


Halfedge HalfedgeMesh::connectVertices(Halfedge heA, Halfedge heB) {

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


std::tuple<Halfedge, Halfedge> HalfedgeMesh::separateEdge(Edge e) {

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


Halfedge HalfedgeMesh::switchHalfedgeSides(Edge e) {

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

bool HalfedgeMesh::ensureEdgeHasInteriorHalfedge(Edge e) {
  if (!e.halfedge().isInterior()) {
    switchHalfedgeSides(e);
    modificationTick++;
    return true;
  }
  return false;
}

/*

Halfedge HalfedgeMesh::connectVertices(Face faceIn, Vertex vAIn, Vertex vBIn) {

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

Vertex HalfedgeMesh::insertVertex(Face fIn) {

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


Vertex HalfedgeMesh::collapseEdge(Edge e) {

  // FIXME I think this function is significantly buggy
  throw std::runtime_error("don't trust this function");

  // Is the edge we're collapsing along the boundary
  bool onBoundary = e.isBoundary();

  // Gather some values
  Halfedge heA0 = e.halfedge();
  Halfedge heA1 = heA0.next();
  Halfedge heA2 = heA1.next();
  Face fA = heA0.face();
  Vertex vA = heA0.vertex();
  GC_SAFETY_ASSERT(heA2.next() == heA0, "face must be triangular to collapse")
  Halfedge heA1T = heA1.twin();
  Halfedge heA1TPrev = heA1T.prevOrbitVertex();
  bool vAOnBoundary = vA.isBoundary();

  Halfedge heB0 = heA0.twin();
  Halfedge heB1 = heB0.next();
  Halfedge heB2 = heB1.next();
  Face fB = heB0.face();
  Vertex vB = heB0.vertex();
  GC_SAFETY_ASSERT(heB2.next() == heB0 || onBoundary, "face must be triangular or on boundary to collapse")
  Halfedge heB2T = heB2.twin();
  Halfedge heB2TPrev = heB2T.prevOrbitVertex();
  bool vBOnBoundary = vB.isBoundary();

  // === Check validity

  // Refuse to do a collapse which removes a boundary component
  // (this could be done, but isn't implemented)
  if (onBoundary && heB2.next() == heB0) {
    return Vertex();
  }

  // Refuse to do a collapse which connects separate boundary components (imagine pinching the neck of an hourglass)
  if (!onBoundary && vBOnBoundary && vAOnBoundary) {
    return Vertex();
  }

  // Should be exactly two vertices, the opposite diamond vertices, in the intersection of the 1-rings.
  // Checking this property ensures that triangulations stays a simplicial complex (no new self-edges, etc).
  std::unordered_set<Vertex> vANeighbors;
  for (Vertex vN : Vertex(vA).adjacentVertices()) {
    vANeighbors.insert(vN);
  }
  size_t nShared = 0;
  for (Vertex vN : Vertex(vB).adjacentVertices()) {
    if (vANeighbors.find(vN) != vANeighbors.end()) {
      nShared++;
    }
  }
  if (nShared > 2) {
    return Vertex();
  }


  // TODO degree 2 vertex case?

  // == Fix connections
  //   - the halfedge heA2 will be repurposed as heA1.twin()
  //   - the halfedge heB1 will be repurposed as heB2.twin()

  // Neighbors of vB
  // for(Halfedge he : vB.outgoingHalfedges()) {
  // heVertexArr[he.getIndex()] = vA.getIndex();
  //}

  // == Around face A
  {
    heNextArr[heA1TPrev.getIndex()] = heA2.getIndex();
    heNextArr[heA2.getIndex()] = heNextArr[heA1T.getIndex()];
    heVertexArr[heA2.getIndex()] = heVertexArr[heA1T.getIndex()];
    heFaceArr[heA2.getIndex()] = heFaceArr[heA1T.getIndex()];

    // Vertex connections
    if (heA1T.vertex().halfedge() == heA1T) {
      vHalfedgeArr[heA1T.vertex().getIndex()] = heA2.getIndex();
    }
    if (vA.halfedge() == heA0) {
      vHalfedgeArr[vA.getIndex()] = heA2.twin().getIndex();
    }

    // Face connections
    if (heA1T.face().halfedge() == heA1T) {
      fHalfedgeArr[heA1T.face().getIndex()] = heA2.getIndex();
    }
  }

  // == Around face B
  if (onBoundary) {
    // The case where we're collapsing a boundary halfedge, just need to decrease the degree of the boundary loop

    Halfedge heB0P = heB0.prevOrbitVertex();
    throw std::runtime_error("not quite implemented");

  } else {
    // The normal case where we're not collapsing a boundary halfedge, similar to what we did at face A

    // Handle halfedge connections around heB2
    heNextArr[heB2TPrev.getIndex()] = heB1.getIndex();
    heNextArr[heB1.getIndex()] = heNextArr[heB2T.getIndex()];
    heVertexArr[heB1.getIndex()] = heVertexArr[heB2T.getIndex()];
    heFaceArr[heB1.getIndex()] = heFaceArr[heB2T.getIndex()];

    // Vertex connections
    // Don't need to update this vertex, since heB2T.vertex() is about to be deleted
    // if (heB2T.vertex().halfedge() == heB2T) {
    // vHalfedgeArr[heB2T.vertex().getIndex()] = heB1.getIndex();
    //}

    // Face connections
    if (heB2T.face().halfedge() == heB2T) {
      fHalfedgeArr[heB2T.face().getIndex()] = heB1.getIndex();
    }
  }

  // Make sure we set the "new" vertex to have an acceptable boundary halfedge if we pulled it on to the boundary
  if (vBOnBoundary) {
    ensureVertexHasBoundaryHalfedge(vA);
  }

  // === Delete the actual elements

  deleteEdgeTriple(heA0);
  deleteEdgeTriple(heA1);
  deleteElement(vB);
  deleteElement(fA);
  if (onBoundary) {
    deleteElement(fB);
    deleteEdgeTriple(heB2);
  }


  modificationTick++;
  return vA;
}

bool HalfedgeMesh::removeFaceAlongBoundary(Face f) {

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
    deleteEdgeTriple(heB);
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

Face HalfedgeMesh::removeVertex(Vertex v) {
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
    if (he.face() != keepFace) {
      deleteElement(he.face());
    }
    deleteEdgeTriple(he);
  }
  deleteElement(v);

  modificationTick++;
  return keepFace;
}


void HalfedgeMesh::ensureVertexHasBoundaryHalfedge(Vertex v) {
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

bool HalfedgeMesh::removeFaceAlongBoundary(Face f) {

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

std::vector<Face> HalfedgeMesh::triangulate(Face f) {
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
    allFaces.emplace_back(neighHalfedges[i].face());
  }

  modificationTick++;
  return allFaces;
}

void HalfedgeMesh::populateVertexIterationCache(bool skipDead) {
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

void HalfedgeMesh::validateConnectivity() {

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
}


Vertex HalfedgeMesh::getNewVertex() {

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

Halfedge HalfedgeMesh::getNewHalfedge(bool isInterior) {

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

Edge HalfedgeMesh::getNewEdge() {

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

Halfedge HalfedgeMesh::getNewEdgeTriple(bool onBoundary) {

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


Face HalfedgeMesh::getNewFace() {

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

BoundaryLoop HalfedgeMesh::getNewBoundaryLoop() {

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

void HalfedgeMesh::expandFaceStorage() {
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


void HalfedgeMesh::deleteEdgeTriple(Halfedge he) {
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

void HalfedgeMesh::deleteElement(Vertex v) {
  size_t iV = v.getIndex();

  vHalfedgeArr[iV] = INVALID_IND;
  nVerticesCount--;

  modificationTick++;
  isCompressedFlag = false;
}

void HalfedgeMesh::deleteElement(Face f) {
  size_t iF = f.getIndex();

  fHalfedgeArr[iF] = INVALID_IND;
  nFacesCount--;

  modificationTick++;
  isCompressedFlag = false;
}

void HalfedgeMesh::deleteElement(BoundaryLoop bl) {
  size_t iF = boundaryLoopIndToFaceInd(bl.getIndex());

  fHalfedgeArr[iF] = INVALID_IND;
  nBoundaryLoopsCount--;

  modificationTick++;
  isCompressedFlag = false;
}


/*

void HalfedgeMesh::deleteElement(Halfedge he) {
  he->markDead();
  isCompressedFlag = false;

  if (he.isReal()) {
    nRealHalfedgesCount--;
  } else {
    nImaginaryHalfedgesCount--;
  }
}

void HalfedgeMesh::deleteElement(Edge e) {
  e->markDead();
  isCompressedFlag = false;
  nEdgesCount--;
}

void HalfedgeMesh::deleteElement(Vertex v) {
  v->markDead();
  isCompressedFlag = false;
  nVerticesCount--;
}

void HalfedgeMesh::deleteElement(Face f) {
  f->markDead();
  isCompressedFlag = false;
  nFacesCount--;
}

void HalfedgeMesh::compressHalfedges() {

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

void HalfedgeMesh::compressEdges() {

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

void HalfedgeMesh::compressFaces() {

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


void HalfedgeMesh::compressVertices() {

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

void HalfedgeMesh::compress() {

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
