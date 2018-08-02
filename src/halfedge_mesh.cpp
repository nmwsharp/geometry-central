#include "geometrycentral/halfedge_mesh.h"

#include <algorithm>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

#include "geometrycentral/geometry.h"
#include <geometrycentral/disjoint_sets.h>
#include <geometrycentral/halfedge_mesh_data_transfer.h>
#include <geometrycentral/polygon_soup_mesh.h>
#include <geometrycentral/timing.h>

using std::cout;
using std::endl;

namespace geometrycentral {

// Cache some basic information that may be queried many
// times, but require O(n) computation to determine.
void HalfedgeMesh::cacheInfo() {
  cache_isSimplicial();
  cache_nFacesTriangulation();
  cache_longestBoundaryLoop();
  cache_nInteriorVertices();
  cache_nConnectedComponents();
}

// Returns true if and only if all faces are triangles
bool HalfedgeMesh::isSimplicial() const { return _isSimplicial; }

void HalfedgeMesh::cache_isSimplicial() {
  _isSimplicial = true;
  for (FacePtr f : faces()) {
    if (f.degree() != 3) {
      _isSimplicial = false;
      return;
    }
  }
}

void HalfedgeMesh::cache_nConnectedComponents() {

  VertexData<size_t> vertInd = getVertexIndices();
  DisjointSets dj(nVertices());
  for (EdgePtr e : edges()) {
    dj.merge(vertInd[e.halfedge().vertex()], vertInd[e.halfedge().twin().vertex()]);
  }
  std::unordered_set<size_t> distinctComponents;
  for (size_t i = 0; i < nVertices(); i++) {
    distinctComponents.insert(dj.find(i));
  }
  _nConnectedComponents = distinctComponents.size();
}

void HalfedgeMesh::cache_nInteriorVertices() {
  _nInteriorVertices = 0;
  for (const VertexPtr v : vertices()) {
    if (!v->isBoundary) {
      _nInteriorVertices++;
    }
  }
}

// Counts the total number of faces in a triangulation of this mesh,
// corresponding to the triangulation given by Face::triangulate().
size_t HalfedgeMesh::nFacesTriangulation() const { return _nFacesTriangulation; }

void HalfedgeMesh::cache_nFacesTriangulation() {
  _nFacesTriangulation = 0;
  for (FacePtr f : faces()) {
    _nFacesTriangulation += f.degree() - 2;
  }
}

size_t HalfedgeMesh::longestBoundaryLoop() const { return _longestBoundaryLoop; }

int HalfedgeMesh::eulerCharacteristic() const {
  // be sure to do intermediate arithmetic with large, signed integers
  return static_cast<int>(static_cast<long long int>(nVertices()) - static_cast<long long int>(nEdges()) +
                          static_cast<long long int>(nFaces()));
}

void HalfedgeMesh::cache_longestBoundaryLoop() {
  int max = 0;
  for (size_t i = 0; i < nBoundaryLoops(); i++) {
    // Count halfedges on boundary and check if greater than max
    int n = 0;
    HalfedgePtr he = boundaryLoop(i).halfedge();
    HalfedgePtr h = he;
    do {
      n++;
      h = h.next();
    } while (h != he);

    if (n > max) {
      max = n;
      _longestBoundaryLoop = i;
    }
  }
}

VertexData<size_t> HalfedgeMesh::getVertexIndices() {
  VertexData<size_t> indices(this);
  size_t i = 0;
  for (VertexPtr v : vertices()) {
    indices[v] = i;
    i++;
  }
  return indices;
}

VertexData<size_t> HalfedgeMesh::getInteriorVertexIndices() {
  VertexData<size_t> indices(this);
  size_t i = 0;
  for (VertexPtr v : vertices()) {
    if (v->isBoundary) {
      indices[v] = -7;
    } else {
      indices[v] = i;
      i++;
    }
  }
  return indices;
}

FaceData<size_t> HalfedgeMesh::getFaceIndices() {
  FaceData<size_t> indices(this);
  size_t i = 0;
  for (FacePtr f : faces()) {
    indices[f] = i;
    i++;
  }
  return indices;
}

EdgeData<size_t> HalfedgeMesh::getEdgeIndices() {
  EdgeData<size_t> indices(this);
  size_t i = 0;
  for (EdgePtr e : edges()) {
    indices[e] = i;
    i++;
  }
  return indices;
}

HalfedgeData<size_t> HalfedgeMesh::getHalfedgeIndices() {
  HalfedgeData<size_t> indices(this);
  size_t i = 0;
  for (HalfedgePtr he : allHalfedges()) {
    indices[he] = i;
    i++;
  }
  return indices;
}

CornerData<size_t> HalfedgeMesh::getCornerIndices() {
  CornerData<size_t> indices(this);
  size_t i = 0;
  for (CornerPtr c : corners()) {
    if (c.halfedge().isReal()) {
      indices[c] = i;
      i++;
    }
  }
  return indices;
}

size_t HalfedgeMesh::nInteriorVertices() const { return _nInteriorVertices; }

size_t HalfedgeMesh::nConnectedComponents() const { return _nConnectedComponents; }

HalfedgeMesh::HalfedgeMesh() : _isSimplicial(true), _nFacesTriangulation(0), _longestBoundaryLoop(0) {}

// Helper for below
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

HalfedgeMesh::HalfedgeMesh(const PolygonSoupMesh& input, Geometry<Euclidean>*& geometry) {
  /*   High-level outline of this algorithm:
   *
   *      0. Count how many of each object we will need so we can pre-allocate
   * them
   *
   *      1. Iterate over the faces of the input mesh, creating edge, face, and
   * halfedge objects
   *
   *      2. Walk around boundaries, marking boundary edge and creating
   * imaginary halfedges/faces
   *
   *      3. Copy the vertex positions to the geometry associated with the new
   * mesh
   *
   */

  START_TIMING(construction)

  // === 0. Count needed objects

  // Find which vertices are actually used
  std::vector<bool> usedVerts(input.vertexCoordinates.size());
  std::vector<size_t> usedVertsIndex(input.vertexCoordinates.size());
  std::fill(usedVerts.begin(), usedVerts.end(), false);
  std::fill(usedVertsIndex.begin(), usedVertsIndex.end(), 0);
  size_t nVerts = 0;
  for (auto poly : input.polygons) {
    for (auto i : poly) {
      usedVerts[i] = true;
    }
  }
  for (size_t i = 0; i < usedVertsIndex.size(); i++) {
    if (usedVerts[i]) {
      usedVertsIndex[i] = nVerts++;
    }
  }

  // === 0.5 Efficiently build an adjacent-face lookup table

  size_t INVALID_IND = std::numeric_limits<size_t>::max();

  // Build a sorted list of (directed halfedge) neighbors of each vertex in
  // compressed format.

  // Count neighbors of each vertex
  size_t nDirected = 0;
  size_t maxPolyDegree = 0;
  std::vector<size_t> vertexNeighborsCount(input.vertexCoordinates.size(), 0);
  std::vector<size_t> vertexNeighborsStart(input.vertexCoordinates.size() + 1);
  for (size_t iFace = 0; iFace < input.polygons.size(); iFace++) {
    auto poly = input.polygons[iFace];
    nDirected += poly.size();
    maxPolyDegree = std::max(maxPolyDegree, poly.size());
    for (size_t j : poly) {
      vertexNeighborsCount[j]++;
    }
  }

  // Build a running sum of the number of neighbors to use a compressed list
  vertexNeighborsStart[0] = 0;
  size_t runningSum = 0;
  for (size_t iVert = 0; iVert < input.vertexCoordinates.size(); iVert++) {
    runningSum += vertexNeighborsCount[iVert];
    vertexNeighborsStart[iVert + 1] = runningSum;
  }

  // Populate the compressed list
  // Each vertex's neighbors are stored between vertexNeighborsStart[i] and
  // vertexNeighborsStart[i+1]
  std::vector<size_t> vertexNeighborsInd(vertexNeighborsStart.begin(), vertexNeighborsStart.end() - 1);
  std::vector<size_t> allVertexNeighbors(nDirected);
  for (size_t iFace = 0; iFace < input.polygons.size(); iFace++) {
    auto poly = input.polygons[iFace];
    for (size_t j = 0; j < poly.size(); j++) {
      size_t fromInd = poly[j];
      size_t toInd = poly[(j + 1) % poly.size()];
      allVertexNeighbors[vertexNeighborsInd[fromInd]] = toInd;
      vertexNeighborsInd[fromInd]++;
    }
  }

  // Sort each of the sublists in the compressed list
  for (size_t iVert = 0; iVert < input.vertexCoordinates.size(); iVert++) {
    std::sort(allVertexNeighbors.begin() + vertexNeighborsStart[iVert],
              allVertexNeighbors.begin() + vertexNeighborsStart[iVert + 1]);
  }

  // Count real and imaginary edges and faces and cache adjacent twin indices
  // Note: counting boundary loops is kinda difficult, so we wait to do so until
  // the final step
  size_t nPairedEdges = 0;
  size_t nUnpairedEdges = 0;
  std::vector<size_t> twinInd(allVertexNeighbors.size());
  for (size_t iVert = 0; iVert < input.vertexCoordinates.size(); iVert++) {
    size_t jStart = vertexNeighborsStart[iVert];
    size_t jEnd = vertexNeighborsStart[iVert + 1];
    for (size_t jInd = jStart; jInd < jEnd; jInd++) {
      size_t jVert = allVertexNeighbors[jInd];

      // Search for the j --> i edge
      size_t searchResult =
          halfedgeLookup(allVertexNeighbors, iVert, vertexNeighborsStart[jVert], vertexNeighborsStart[jVert + 1]);
      twinInd[jInd] = searchResult;
      if (searchResult == INVALID_IND) {
        nUnpairedEdges++;
      } else {
        nPairedEdges++;
      }
    }
  }
  nPairedEdges /= 2;

  size_t nTotalEdges = nPairedEdges + nUnpairedEdges;
  nRealHalfedges = 2 * nPairedEdges + nUnpairedEdges;
  size_t nImaginaryHalfedges = nUnpairedEdges;
  size_t nRealFaces = input.polygons.size();

  // Allocate space and construct elements
  rawHalfedges.reserve(nRealHalfedges + nImaginaryHalfedges);
  for (size_t i = 0; i < (nRealHalfedges + nImaginaryHalfedges); i++) getNewHalfedge();
  rawVertices.reserve(nVerts);
  for (size_t i = 0; i < nVerts; i++) getNewVertex();
  rawEdges.reserve(nTotalEdges);
  for (size_t i = 0; i < nTotalEdges; i++) getNewEdge();
  rawFaces.reserve(nRealFaces);
  for (size_t i = 0; i < nRealFaces; i++) getNewFace();

  // === 1. Create faces, edges, and halfedges

  // Keep track of the edges we've already created since we only need one per
  // edge
  std::vector<EdgePtr> sharedEdges(twinInd.size(), EdgePtr());

  // Iterate over faces
  size_t iFace = 0;
  size_t iEdge = 0;
  size_t iHalfedge = 0;
  std::vector<HalfedgePtr> thisFaceHalfedges;
  thisFaceHalfedges.reserve(maxPolyDegree);
  for (auto poly : input.polygons) {
    size_t degree = poly.size();

    // Create a new face object
    FacePtr f{&rawFaces[iFace]};
    f->isReal = true;

    // The halfedges that make up this face
    // std::vector<HalfedgePtr> thisFaceHalfedges(degree);
    thisFaceHalfedges.resize(degree);

    for (size_t iPolyEdge = 0; iPolyEdge < degree; iPolyEdge++) {
      HalfedgePtr he{&rawHalfedges[iHalfedge]};
      iHalfedge++;

      size_t ind1 = poly[iPolyEdge];
      size_t ind2 = poly[(iPolyEdge + 1) % degree];

      // Connect up pointers
      he->vertex = &rawVertices[usedVertsIndex[ind1]];
      he->vertex->halfedge = he.ptr;
      he->face = f.ptr;
      thisFaceHalfedges[iPolyEdge] = he;
      he->twin = nullptr; // ensure this is null so we can detect boundaries below

      // Get a reference to the edge shared by this and its twin, creating the
      // object if needed
      size_t myHeInd =
          halfedgeLookup(allVertexNeighbors, ind2, vertexNeighborsStart[ind1], vertexNeighborsStart[ind1 + 1]);
      size_t twinHeInd = twinInd[myHeInd];
      bool edgeAlreadyCreated = (twinHeInd != INVALID_IND) && (sharedEdges[twinHeInd] != EdgePtr());

      if (edgeAlreadyCreated) {
        EdgePtr sharedEdge = sharedEdges[twinHeInd];
        he->edge = sharedEdge.ptr;
        he->twin = sharedEdge->halfedge;
        he->twin->twin = he.ptr;
      } else {
        EdgePtr sharedEdge{&rawEdges[iEdge]};
        iEdge++;
        sharedEdge->halfedge = he.ptr;
        he->edge = sharedEdge.ptr;
        sharedEdges[myHeInd] = sharedEdge;
      }
    }

    // Do one more lap around the face to set next pointers
    for (size_t iPolyEdge = 0; iPolyEdge < degree; iPolyEdge++) {
      thisFaceHalfedges[iPolyEdge]->next = thisFaceHalfedges[(iPolyEdge + 1) % degree].ptr;
    }

    f->halfedge = thisFaceHalfedges[0].ptr;
    thisFaceHalfedges.clear();
    iFace++;
  }

  // === 2. Walk the boundary to find/create boundary cycles

  // First, do a pre-walk to count the boundary loops we will need and allocate
  // them
  size_t nBoundaryLoops = 0;
  std::set<HalfedgePtr> walkedHalfedges;
  for (size_t iHe = 0; iHe < nHalfedges(); iHe++) {
    if (halfedge(iHe)->twin == nullptr && walkedHalfedges.find(halfedge(iHe)) == walkedHalfedges.end()) {
      nBoundaryLoops++;
      HalfedgePtr currHe = halfedge(iHe);
      walkedHalfedges.insert(currHe);
      size_t walkCount = 0;
      do {
        currHe = currHe->next;
        while (currHe->twin != nullptr) {
          currHe = currHe->twin->next;
          walkCount++;
          if (walkCount > nHalfedges()) {
            throw std::runtime_error(
                "Encountered infinite loop while constructing halfedge mesh. Are you sure the input is manifold?");
          }
        }
        walkedHalfedges.insert(currHe);
      } while (currHe != halfedge(iHe));
    }
  }
  rawBoundaryLoops.resize(nBoundaryLoops);

  // Now do the actual walk in which we construct and connect objects
  size_t iBoundaryLoop = 0;
  size_t iImaginaryHalfedge = 0;
  for (size_t iHe = 0; iHe < nHalfedges(); iHe++) {
    // Note: If distinct holes share a given vertex, this algorithm will see
    // them as a single "figure 8-like"
    // hole and connect them with a single imaginary face.
    // TODO: fix this?

    // If this halfedge doesn't have a twin, it must be on a boundary (or have
    // already been processed while walking a hole)
    if (halfedge(iHe)->twin == nullptr) {
      // Create a boundary loop for this hole
      BoundaryPtr boundaryLoop{&rawBoundaryLoops[iBoundaryLoop]};
      boundaryLoop->isReal = false;

      // Walk around the boundary loop, creating imaginary halfedges
      HalfedgePtr currHe = halfedge(iHe);
      HalfedgePtr prevHe{nullptr};
      bool finished = false;
      while (!finished) {
        // Create a new, imaginary halfedge
        HalfedgePtr newHe{&rawHalfedges[nRealHalfedges + iImaginaryHalfedge]};
        boundaryLoop->halfedge = newHe.ptr;
        iImaginaryHalfedge++;

        // Connect up pointers
        newHe->isReal = false;
        newHe->twin = currHe.ptr;
        currHe->twin = newHe.ptr;
        newHe->face = boundaryLoop.ptr;
        newHe->edge = currHe->edge;
        newHe->vertex = currHe->next->vertex;

        // Some pointers need values only visible from the previous iteration of
        // the loop.
        // The first one we process gets missed, handle it at the end outside
        // the loop.
        if (prevHe == nullptr) {
          newHe->next = nullptr;
        } else {
          newHe->next = prevHe->twin;
        }

        // Set the isBoundary property where appropriate
        currHe->face->isBoundary = true;
        currHe->vertex->isBoundary = true;
        currHe->edge->isBoundary = true;

        // Prepare for the next iteration
        prevHe = currHe;
        currHe = currHe->next;
        while (currHe->twin != nullptr) {
          // When we've finished walking the loop, we'll be able to tell
          // because we spin in circles around the vertex trying to continue
          // the loop. Detect that here and quit.
          if (currHe->twin->next == nullptr) {
            finished = true;
            break;
          }

          currHe = currHe->twin->next;
        }
      }

      // As noted above, the pointers don't get set properly on the first
      // iteration of
      // the loop above because we don't have a reference to prev yet. Fix that
      // here.
      halfedge(iHe)->twin->next = prevHe->twin;

      iBoundaryLoop++;
    }
  }


// When in debug mode, mesh elements know what mesh they are a part of so
// we can do assertions for saftey checks.
#ifndef NDEBUG
  for (HalfedgePtr x : allHalfedges()) {
    x->parentMesh = this;
  }
  for (VertexPtr x : vertices()) {
    x->parentMesh = this;
  }
  for (EdgePtr x : edges()) {
    x->parentMesh = this;
  }
  for (FacePtr x : faces()) {
    x->parentMesh = this;
  }
  for (FacePtr x : boundaryLoops()) {
    x->parentMesh = this;
  }
#endif

  // === 3. Map vertices in the halfedge mesh to the associated vertex
  // coordinates in space
  // Create the vertex objects and build a map to find them
  size_t iVert = 0;
  geometry = new Geometry<Euclidean>(*this);
  for (size_t i = 0; i < input.vertexCoordinates.size(); i++) {
    if (usedVerts[i]) {
      geometry->position(vertex(iVert)) = input.vertexCoordinates[i];
      iVert++;
    }
  }

  // Print some nice statistics
  std::cout << "Constructed halfedge mesh with: " << std::endl;
  std::cout << "    # verts =  " << nVertices() << std::endl;
  std::cout << "    # edges =  " << nEdges() << std::endl;
  std::cout << "    # faces =  " << nFaces() << std::endl;
  std::cout << "    # halfedges =  " << nHalfedges() << std::endl;
  std::cout << "      and " << nBoundaryLoops << " boundary components. " << std::endl;
  std::cout << "Construction took " << pretty_time(FINISH_TIMING(construction)) << std::endl;

  // Compute some basic information about the mesh
  cacheInfo();
}

bool Edge::flip() {
  //                         b b
  //                         * *
  //                        /|\                                                          / \
//                       / | \                                                        /   \
//                      /  |  \                                                      /     \
//                     /   |   \                                                    /       \
//                    /    |    \          ============ FLIP ==========>           /         \
//                   /     |     \                                                /           \
//                  /      |      \                                              /             \
//                 /       |       \                                            /               \
//                /  /     |    |\  \                                          /  /          |\  \
//               /  /      |    | \  \                                        /  /           | \  \
//              /  /       |       \  \                                      /  /               \  \
//             /  /        |        \  \                                    /  /                 \  \
//            /  /         |         \  \                                  /  /        fa         \  \
//           /  /          |          \  \                                /  /                     \  \
//          /  / ha2       |       hb3 \  \                              /  / ha2               hb3 \  \
//         /  /        /|  |  |         \  \                            /  /                         \  \
//        /  /        / |  |  |          \  \                          /  /                           \  \
//       /  /           |  |  |           \  \                        /  /___                          \  \
//      /  /___         |  |  |            \  \                      /                 ha1        \        \
//     /                |  |  |                \                    /     _________________________\        \
//    /                 |  |  |                 \                  /                                         \
// c *      fa      ha1 |  e  | hb1      fb      * d            c *-------------------- e --------------------* d
  //    \                 |  |  |                 /                  \
  //    __________________________         /
  //     \                |  |  |          ___   /                    \     \
  //     hb1                  /
  //      \  \            |  |  |            /  /                      \     \
  //      ___   /
  //       \  \           |  |  |           /  /                        \ /  /
  //        \  \          |  |  | /        /  /                          \  \
  //        /  /
  //         \  \ ha3     |  |  |/    hb2 /  /                            \  \
  //         ha3                 hb2 /  /
  //          \  \           |           /  /                              \  \
  //          /  /
  //           \  \          |          /  /                                \  \
  //           fb          /  /
  //            \  \         |         /  /                                  \
  //            \                   /  /
  //             \  \        |        /  /                                    \
  //             \                 /  /
  //              \  \       |       /  /                                      \
  //              \               /  /
  //               \  \ |    |      /  / \  \             /  /
  //                \  \|    |     /  / \  \ |         /  /
  //                 \       |       / \  \|           /
  //                  \      |      / \             /
  //                   \     |     / \           /
  //                    \    |    / \         /
  //                     \   |   / \       /
  //                      \  |  / \     /
  //                       \ | / \   /
  //                        \|/ \ /
  //                         * *
  //                         a a

  Edge* e = this;

  // Get halfedges of first face
  Halfedge* ha1 = e->halfedge;
  if (!ha1->isReal) return false; // don't flip boundary edges
  Halfedge* ha2 = ha1->next;
  Halfedge* ha3 = ha2->next;
  if (ha3->next != ha1) return false; // not a triangle

  // Get halfedges of second face
  Halfedge* hb1 = ha1->twin;
  if (!hb1->isReal) return false; // don't flip boundary edges
  Halfedge* hb2 = hb1->next;
  Halfedge* hb3 = hb2->next;
  if (hb3->next != hb1) return false; // not a triangle

  // Get vertices and faces
  Vertex* va = ha1->vertex;
  Vertex* vb = hb1->vertex;
  Vertex* vc = ha3->vertex;
  Vertex* vd = hb3->vertex;
  Face* fa = ha1->face;
  Face* fb = hb1->face;

  // Update vertex pointers
  if (va->halfedge == ha1) va->halfedge = hb2;
  if (vb->halfedge == hb1) vb->halfedge = ha2;
  // (vc and vd can't be invalidated by the flip)

  // Update edge pointers
  // (e still has the same halfedges)

  // Update face pointers
  fa->halfedge = ha1;
  fb->halfedge = hb1;

  // Update halfedge pointers
  ha1->next = hb3;
  hb3->next = ha2;
  ha2->next = ha1;
  hb1->next = ha3;
  ha3->next = hb2;
  hb2->next = hb1;
  ha1->vertex = vc;
  hb1->vertex = vd;
  ha3->face = fb;
  hb3->face = fa;

  return true;
}

HalfedgeMesh* HalfedgeMesh::copy() {

  HalfedgeMeshDataTransfer t;
  return copy(t);
}


HalfedgeMesh* HalfedgeMesh::copy(HalfedgeMeshDataTransfer& dataTransfer) {

  HalfedgeMesh* newMesh = new HalfedgeMesh();

  // Copy vectors
  newMesh->rawHalfedges = rawHalfedges;
  newMesh->rawVertices = rawVertices;
  newMesh->rawEdges = rawEdges;
  newMesh->rawFaces = rawFaces;
  newMesh->rawBoundaryLoops = rawBoundaryLoops;
  newMesh->nRealHalfedges = nRealHalfedges;

  // Create the data transfer object
  dataTransfer = HalfedgeMeshDataTransfer(this, newMesh);

  // Build maps
  for (size_t i = 0; i < rawHalfedges.size(); i++) {
    dataTransfer.heMap[halfedge(i)] = &newMesh->rawHalfedges[i];
  }
  for (size_t i = 0; i < rawVertices.size(); i++) {
    dataTransfer.vMap[vertex(i)] = &newMesh->rawVertices[i];
  }
  for (size_t i = 0; i < rawEdges.size(); i++) {
    dataTransfer.eMap[edge(i)] = &newMesh->rawEdges[i];
  }
  for (size_t i = 0; i < rawFaces.size(); i++) {
    dataTransfer.fMap[face(i)] = &newMesh->rawFaces[i];
  }
  for (size_t i = 0; i < rawBoundaryLoops.size(); i++) {
    dataTransfer.fMap[boundaryLoop(i)] = &newMesh->rawBoundaryLoops[i];
  }

  // Shift pointers
  for (size_t i = 0; i < rawHalfedges.size(); i++) {
    newMesh->rawHalfedges[i].next = dataTransfer.heMap[halfedge(i).next()].ptr;
    newMesh->rawHalfedges[i].twin = dataTransfer.heMap[halfedge(i).twin()].ptr;
    newMesh->rawHalfedges[i].vertex = dataTransfer.vMap[halfedge(i).vertex()].ptr;
    newMesh->rawHalfedges[i].edge = dataTransfer.eMap[halfedge(i).edge()].ptr;
    newMesh->rawHalfedges[i].face = dataTransfer.fMap[halfedge(i).face()].ptr;
  }
  for (size_t i = 0; i < rawVertices.size(); i++) {
    newMesh->rawVertices[i].halfedge = dataTransfer.heMap[vertex(i).halfedge()].ptr;
  }
  for (size_t i = 0; i < rawEdges.size(); i++) {
    newMesh->rawEdges[i].halfedge = dataTransfer.heMap[edge(i).halfedge()].ptr;
  }
  for (size_t i = 0; i < rawFaces.size(); i++) {
    newMesh->rawFaces[i].halfedge = dataTransfer.heMap[face(i).halfedge()].ptr;
  }
  for (size_t i = 0; i < rawBoundaryLoops.size(); i++) {
    newMesh->rawBoundaryLoops[i].halfedge = dataTransfer.heMap[boundaryLoop(i).halfedge()].ptr;
  }


#ifndef NDEBUG
  // Set the parent mesh pointers to point to the correct mesh
  for (size_t i = 0; i < newMesh->rawHalfedges.size(); i++) {
    newMesh->rawHalfedges[i].parentMesh = newMesh;
  }
  for (size_t i = 0; i < newMesh->rawVertices.size(); i++) {
    newMesh->rawVertices[i].parentMesh = newMesh;
  }
  for (size_t i = 0; i < newMesh->rawEdges.size(); i++) {
    newMesh->rawEdges[i].parentMesh = newMesh;
  }
  for (size_t i = 0; i < newMesh->rawFaces.size(); i++) {
    newMesh->rawFaces[i].parentMesh = newMesh;
  }
  for (size_t i = 0; i < newMesh->rawBoundaryLoops.size(); i++) {
    newMesh->rawBoundaryLoops[i].parentMesh = newMesh;
  }
#endif

  newMesh->cacheInfo();
  dataTransfer.generateReverseMaps();

  return newMesh;
}

std::vector<std::vector<size_t>> HalfedgeMesh::getPolygonSoupFaces() {

  std::vector<std::vector<size_t>> result;

  VertexData<size_t> vInd = getVertexIndices();
  for (FacePtr f : faces()) {
    std::vector<size_t> faceList;
    for (VertexPtr v : f.adjacentVertices()) {
      faceList.push_back(vInd[v]);
    }
    result.push_back(faceList);
  }

  return result;
}

HalfedgePtr HalfedgeMesh::insertVertexAlongEdge(EdgePtr e) {

  // == Gather / create elements
  // Faces are identified as 'A', and 'B'

  // Create first, because getNew() could invalidate pointers
  Vertex* newV = getNewVertex();
  Edge* newE = getNewEdge();
  Halfedge* heANew = getNewHalfedge();
  Halfedge* heBNew = getNewHalfedge();

  Halfedge* heACenter = e.halfedge().ptr;
  Halfedge* heBCenter = heACenter->twin;
  // Halfedge* heANext = heACenter->next;
  Halfedge* heBNext = heBCenter->next;
  Halfedge* heAPrev = HalfedgePtr{heACenter}.prev().ptr;
  // Halfedge* heBPrev = HalfedgePtr{heBCenter}.prev().ptr;
  Face* fA = heACenter->face;
  Face* fB = heBCenter->face;
  Vertex* oldVBottom = heACenter->vertex;

  // == Hook up all the pointers

  // New vertex
  newV->halfedge = heACenter;

  // New edge
  newE->halfedge = heANew;

  // New halfedge A
  heANew->twin = heBNew;
  heANew->next = heACenter;
  heANew->vertex = oldVBottom;
  heANew->edge = newE;
  heANew->face = fA;

  // New halfedge B
  heBNew->twin = heANew;
  heBNew->next = heBNext;
  heBNew->vertex = newV;
  heBNew->edge = newE;
  heBNew->face = fB;

  // Fix pointers for old halfedges
  heBCenter->next = heBNew;
  heAPrev->next = heANew;
  oldVBottom->halfedge = heBNext;
  heACenter->vertex = newV;


  return HalfedgePtr{heBNew};
}


VertexPtr HalfedgeMesh::splitEdge(EdgePtr e) {

  // Validate that faces are triangular and real
  if (e.isBoundary() || e.halfedge().face().degree() != 3 || e.halfedge().twin().face().degree() != 3 ||
      e.halfedge().face() == e.halfedge().twin().face()) {
    throw std::logic_error("Can only split non-boundary edge which borders two distinct triangular faces");
  }

  // First operation: insert a new vertex along the edge
  VertexPtr newV = insertVertexAlongEdge(e);


  // Second operation: connect both of the new faces
  FacePtr fOppA = newV.halfedge().face();
  VertexPtr vOppA = newV.halfedge().next().next().vertex();
  FacePtr fOppB = newV.halfedge().twin().face();
  VertexPtr vOppB = newV.halfedge().twin().next().next().next().vertex();

  connectVertices(fOppA, vOppA, newV);
  connectVertices(fOppB, vOppB, newV);

  return newV;
}

EdgePtr HalfedgeMesh::connectVertices(VertexPtr vA, VertexPtr vB) {

  // Find the shared face and call the main version
  std::unordered_set<FacePtr> aFaces;
  for (FacePtr f : vA.adjacentFaces()) {
    aFaces.insert(f);
  }
  FacePtr sharedFace = FacePtr();
  for (FacePtr f : vB.adjacentFaces()) {
    if (aFaces.find(f) != aFaces.end()) {
      sharedFace = f;
      break;
    }
  }

  if (sharedFace == FacePtr()) {
    throw std::logic_error("Vertices do not contain shared face");
  }

  return connectVertices(sharedFace, vA, vB);
}

EdgePtr HalfedgeMesh::connectVertices(FacePtr face, VertexPtr vA, VertexPtr vB) {

  // == Create new elements
  Halfedge* heANew = getNewHalfedge();
  Halfedge* heBNew = getNewHalfedge();
  Edge* eNew = getNewEdge();
  Face* fB = getNewFace();

  // == Find useful halfedges around the face
  Halfedge* heANext;
  Halfedge* heBNext;
  Halfedge* heAPrev;
  Halfedge* heBPrev;
  for (HalfedgePtr he : face.adjacentHalfedges()) {
    if (he.vertex() == vA) {
      heANext = he.ptr;
    }
    if (he.vertex() == vB) {
      heBNext = he.ptr;
    }
  }
  for (HalfedgePtr he : face.adjacentHalfedges()) {
    if (he.next().ptr == heANext) {
      heAPrev = he.ptr;
    }
    if (he.next().ptr == heBNext) {
      heBPrev = he.ptr;
    }
  }

  // == Gather other elements
  Face* fA = heBNext->face;
  Vertex* vAp = vA.ptr;
  Vertex* vBp = vB.ptr;

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


  return EdgePtr{eNew};
}


VertexPtr HalfedgeMesh::insertVertex(FacePtr f) {

  // Create the new center vertex
  VertexPtr centerVert{getNewVertex()};
  size_t faceDegree = 0;
  for (HalfedgePtr he : f.adjacentHalfedges()) {
    faceDegree++;
  }

  // == Create new halfedges/edges/faces around the center vertex

  // Create all of the new elements first, then hook them up below, as this can invalidate pointers.
  std::vector<Face*> innerFaces;
  std::vector<Halfedge*> leadingHalfedges; // the one that points towards the center
  std::vector<Halfedge*> trailingHalfedges;
  std::vector<Edge*> innerEdges; // aligned with leading he
  for (size_t i = 0; i < faceDegree; i++) {
    // Re-use first face
    if (i == 0) {
      innerFaces.push_back(f.ptr);
    } else {
      innerFaces.push_back(getNewFace());
    }

    leadingHalfedges.push_back(getNewHalfedge());
    trailingHalfedges.push_back(getNewHalfedge());
    innerEdges.push_back(getNewEdge());
  }

  // Form this list before we start, because we're about to start breaking pointers
  std::vector<Halfedge*> faceBoundaryHalfedges;
  for (HalfedgePtr he : f.adjacentHalfedges()) {
    faceBoundaryHalfedges.push_back(he.ptr);
  }

  // Connect up all the pointers
  // Each iteration processes one inner face
  for (size_t i = 0; i < faceDegree; i++) {

    // Gather pointers
    Face* f = innerFaces[i];
    Edge* e = innerEdges[i];
    Edge* prevE = innerEdges[(i + faceDegree - 1) & faceDegree];
    Halfedge* leadingHe = leadingHalfedges[i];
    Halfedge* trailingHe = trailingHalfedges[i];
    Halfedge* boundaryHe = faceBoundaryHalfedges[i];
    Halfedge* nextTrailingHe = trailingHalfedges[(i + 1) % faceDegree];
    Halfedge* prevLeadingHe = leadingHalfedges[(i + faceDegree - 1) % faceDegree];

    // face
    f->halfedge = boundaryHe;
    f->isReal = true;

    // edge
    e->halfedge = leadingHe;

    // leading halfedge
    leadingHe->twin = nextTrailingHe;
    leadingHe->next = trailingHe;
    leadingHe->vertex = boundaryHe->next->vertex;
    leadingHe->edge = e;
    leadingHe->face = f;

    // trailing halfedge
    trailingHe->twin = prevLeadingHe;
    trailingHe->next = boundaryHe;
    trailingHe->vertex = centerVert.ptr;
    trailingHe->edge = prevE;
    trailingHe->face = f;
  }

  return centerVert;
}

Halfedge* HalfedgeMesh::getNewHalfedge() {

  // The boring case, when no resize is needed
  if (rawHalfedges.size() < rawHalfedges.capacity()) {
    rawHalfedges.emplace_back();
  }
  // The intesting case, where the vector resizes and we need to update pointers.
  else {

    // Create a new halfedge, allowing the list to expand
    Halfedge* oldStart = &rawHalfedges.front();
    rawHalfedges.emplace_back();
    Halfedge* newStart = &rawHalfedges.front();
    std::ptrdiff_t shift = newStart - oldStart;

    // Shift all pointers
    for (Halfedge& he : rawHalfedges) {
      he.twin += shift;
      he.next += shift;
    }
    for (Vertex& v : rawVertices) {
      v.halfedge += shift;
    }
    for (Edge& e : rawEdges) {
      e.halfedge += shift;
    }
    for (Face& f : rawFaces) {
      f.halfedge += shift;
    }
    for (Face& f : rawBoundaryLoops) {
      f.halfedge += shift;
    }
    for (DynamicHalfedgePtr* d : registeredHalfedgePtrs) {
      d->ptr += shift;
    }
  }

  rawHalfedges.back().ID = nextElemID++;
  return &rawHalfedges.back();
}

Vertex* HalfedgeMesh::getNewVertex() {

  // The boring case, when no resize is needed
  if (rawVertices.size() < rawVertices.capacity()) {
    rawVertices.emplace_back();
  }
  // The intesting case, where the vector resizes and we need to update pointers.
  else {

    // Create a new halfedge, allowing the list to expand
    Vertex* oldStart = &rawVertices.front();
    rawVertices.emplace_back();
    Vertex* newStart = &rawVertices.front();
    std::ptrdiff_t shift = newStart - oldStart;

    // Shift all pointers
    for (Halfedge& he : rawHalfedges) {
      he.vertex += shift;
    }
    for (DynamicVertexPtr* d : registeredVertexPtrs) {
      d->ptr += shift;
    }
  }

  rawVertices.back().ID = nextElemID++;
  return &rawVertices.back();
}

Edge* HalfedgeMesh::getNewEdge() {

  // The boring case, when no resize is needed
  if (rawEdges.size() < rawEdges.capacity()) {
    rawEdges.emplace_back();
  }
  // The intesting case, where the vector resizes and we need to update pointers.
  else {

    // Create a new halfedge, allowing the list to expand
    Edge* oldStart = &rawEdges.front();
    rawEdges.emplace_back();
    Edge* newStart = &rawEdges.front();
    std::ptrdiff_t shift = newStart - oldStart;

    // Shift all pointers
    for (Halfedge& he : rawHalfedges) {
      he.edge += shift;
    }
    for (DynamicEdgePtr* d : registeredEdgePtrs) {
      d->ptr += shift;
    }
  }

  rawEdges.back().ID = nextElemID++;
  return &rawEdges.back();
}

Face* HalfedgeMesh::getNewFace() {

  // The boring case, when no resize is needed
  if (rawFaces.size() < rawFaces.capacity()) {
    rawFaces.emplace_back();
  }
  // The intesting case, where the vector resizes and we need to update pointers.
  else {

    // Create a new halfedge, allowing the list to expand
    Face* oldStart = &rawFaces.front();
    rawFaces.emplace_back();
    Face* newStart = &rawFaces.front();
    std::ptrdiff_t shift = newStart - oldStart;

    // Shift all pointers
    for (Halfedge& he : rawHalfedges) {
      he.face += shift;
    }
    for (DynamicFacePtr* d : registeredFacePtrs) {
      d->ptr += shift;
    }
  }

  rawFaces.back().ID = nextElemID++;
  return &rawFaces.back();
}

void HalfedgeMesh::permuteToCanonical() { throw std::logic_error("not implemented"); }


} // namespace geometrycentral
