#include "geometrycentral/surface/common_subdivision.h"

namespace geometrycentral {
namespace surface {

// helpers
namespace {
template <typename E, typename T>
std::vector<T> toStdVector(MeshData<E, T> data) {
  const Eigen::Matrix<T, Eigen::Dynamic, 1>& dataVals = data.raw();
  std::vector<T> vec;
  vec.reserve(dataVals.size());
  for (size_t i = 0; i < (size_t)dataVals.size(); ++i) vec.push_back(dataVals(i));
  return vec;
}

template <typename E, typename T>
MeshData<E, T> fromStdVector(SurfaceMesh& mesh, const std::vector<T>& data) {
  Eigen::Matrix<T, Eigen::Dynamic, 1> dataVals(data.size());
  for (size_t i = 0; i < (size_t)data.size(); ++i) dataVals(i) = data[i];
  return MeshData<E, T>(mesh, dataVals);
}
} // namespace

std::ostream& operator<<(std::ostream& out, const CSIntersectionType& type) {
  switch (type) {
  case CSIntersectionType::VERTEX_VERTEX:
    out << "Vertex-Vertex intersection";
    break;
  case CSIntersectionType::EDGE_TRANSVERSE:
    out << "Edge-Edge intersection (transverse)";
    break;
  case CSIntersectionType::EDGE_PARALLEL:
    out << "Edge-Edge 'intersection' (parallel)";
    break;
  case CSIntersectionType::FACE_VERTEX:
    out << "Face-Vertex intersection ";
    break;
  case CSIntersectionType::EDGE_VERTEX:
    out << "Edge-Vertex intersection ";
    break;
  }
  return out;
}

std::ostream& operator<<(std::ostream& out, const CommonSubdivisionPoint& pt) {
  out << "CommonSubdivisionPoint{ intersectionType: " << pt.intersectionType << "posA: " << pt.posA
      << ", posB: " << pt.posB << ", orientation: " << pt.orientation << "}";
  return out;
}

CommonSubdivision::CommonSubdivision(ManifoldSurfaceMesh& meshA_, ManifoldSurfaceMesh& meshB_)
    : meshA(meshA_), meshB(meshB_) {

  pointsAlongA = EdgeData<std::vector<CommonSubdivisionPoint*>>(meshA);
  pointsAlongB = EdgeData<std::vector<CommonSubdivisionPoint*>>(meshB);
}

// precondition: must lie along e
// must not be a face point
int CommonSubdivision::getOrderAlongEdgeA(const CommonSubdivisionPoint& p, Edge eA) {
  // finds index in pointAlongA[posA.edge]

  const std::vector<CommonSubdivisionPoint*>& path = pointsAlongA[eA];

  for (size_t iP = 0; iP < path.size(); ++iP) {
    if (&p == path[iP]) {
      return iP;
    }
  }

  return -1;
}

int CommonSubdivision::getOrderAlongEdgeB(const CommonSubdivisionPoint& p, Edge eB) {
  const std::vector<CommonSubdivisionPoint*>& path = pointsAlongB[eB];

  for (size_t iP = 0; iP < path.size(); ++iP) {
    if (&p == path[iP]) {
      return iP;
    }
  }

  return -1;
}

int CommonSubdivision::getIndex(const CommonSubdivisionPoint& p) { return getIndex(&p); }

int CommonSubdivision::getIndex(const CommonSubdivisionPoint* p) {
  // finds index in subdivisionPoints
  for (size_t iP = 0; iP < subdivisionPoints.size(); ++iP) {
    if (p == &subdivisionPoints[iP]) {
      return iP;
    }
  }

  return -1;
}

SparseMatrix<double> CommonSubdivision::interpolationMatrixA() {
  // TODO: could do without constructing mesh
  checkMeshConstructed();

  VertexData<size_t> csVertInd = mesh->getVertexIndices();

  SparseMatrix<double> P_A(mesh->nVertices(), meshA.nVertices());
  std::vector<Eigen::Triplet<double>> tripletList;
  VertexData<size_t> AVertInd = meshA.getVertexIndices();

  for (Vertex v : mesh->vertices()) {

    CommonSubdivisionPoint& p = *sourcePoints[v];
    size_t iP = csVertInd[v];

    switch (p.posA.type) {
    case SurfacePointType::Vertex: {
      tripletList.emplace_back(iP, AVertInd[p.posA.vertex], 1.);
      break;
    }
    case SurfacePointType::Edge: {
      Halfedge he = p.posA.edge.halfedge();
      tripletList.emplace_back(iP, AVertInd[he.tailVertex()], 1. - p.posA.tEdge);
      tripletList.emplace_back(iP, AVertInd[he.tipVertex()], p.posA.tEdge);
      break;
    }
    case SurfacePointType::Face: {
      Halfedge he = p.posA.face.halfedge();
      tripletList.emplace_back(iP, AVertInd[he.vertex()], p.posA.faceCoords[0]);
      he = he.next();
      tripletList.emplace_back(iP, AVertInd[he.vertex()], p.posA.faceCoords[1]);
      he = he.next();
      tripletList.emplace_back(iP, AVertInd[he.vertex()], p.posA.faceCoords[2]);
      break;
    }
    }
  }

  P_A.setFromTriplets(tripletList.begin(), tripletList.end());
  return P_A;
}

SparseMatrix<double> CommonSubdivision::interpolationMatrixB() {
  // TODO: could do without constructing mesh
  checkMeshConstructed();

  VertexData<size_t> csVertInd = mesh->getVertexIndices();

  SparseMatrix<double> P_B(mesh->nVertices(), meshB.nVertices());
  std::vector<Eigen::Triplet<double>> tripletList;
  VertexData<size_t> BVertInd = meshB.getVertexIndices();

  for (Vertex v : mesh->vertices()) {

    CommonSubdivisionPoint& p = *sourcePoints[v];
    size_t iP = csVertInd[v];

    switch (p.posB.type) {
    case SurfacePointType::Vertex: {
      tripletList.emplace_back(iP, BVertInd[p.posB.vertex], 1.);
      break;
    }
    case SurfacePointType::Edge: {
      Halfedge he = p.posB.edge.halfedge();
      tripletList.emplace_back(iP, BVertInd[he.tailVertex()], 1. - p.posB.tEdge);
      tripletList.emplace_back(iP, BVertInd[he.tipVertex()], p.posB.tEdge);
      break;
    }
    case SurfacePointType::Face: {
      Halfedge he = p.posB.face.halfedge();
      tripletList.emplace_back(iP, BVertInd[he.vertex()], p.posB.faceCoords[0]);
      he = he.next();
      tripletList.emplace_back(iP, BVertInd[he.vertex()], p.posB.faceCoords[1]);
      he = he.next();
      tripletList.emplace_back(iP, BVertInd[he.vertex()], p.posB.faceCoords[2]);
      break;
    }
    }
  }

  P_B.setFromTriplets(tripletList.begin(), tripletList.end());
  return P_B;
}

EdgeData<double> CommonSubdivision::interpolateEdgeLengthsA(const EdgeData<double>& lengthA) {
  checkMeshConstructed();

  auto fEdgeLengths = [&](Face f) -> Vector3 { // Gather face's edge lengths
    return Vector3{lengthA[f.halfedge().next().edge()], lengthA[f.halfedge().next().next().edge()],
                   lengthA[f.halfedge().edge()]};
  };

  EdgeData<double> lengthCS(*mesh);
  for (Edge e : mesh->edges()) {
    SurfacePoint tail = sourcePoints[e.halfedge().tailVertex()]->posA;
    SurfacePoint tip = sourcePoints[e.halfedge().tipVertex()]->posA;
    Face f = sharedFace(tail, tip);
    GC_SAFETY_ASSERT(f != Face(), "common subdivision edges must be contained in a face");
    tail = tail.inFace(f);
    tip = tip.inFace(f);

    lengthCS[e] = displacementLength(tip.faceCoords - tail.faceCoords, fEdgeLengths(f));
  }

  return lengthCS;
}

EdgeData<double> CommonSubdivision::interpolateEdgeLengthsB(const EdgeData<double>& lengthB) {
  checkMeshConstructed();

  auto fEdgeLengths = [&](Face f) -> Vector3 { // Gather face's edge lengths
    return Vector3{lengthB[f.halfedge().next().edge()], lengthB[f.halfedge().next().next().edge()],
                   lengthB[f.halfedge().edge()]};
  };

  EdgeData<double> lengthCS(*mesh);
  for (Edge e : mesh->edges()) {
    SurfacePoint tail = sourcePoints[e.halfedge().tailVertex()]->posB;
    SurfacePoint tip = sourcePoints[e.halfedge().tipVertex()]->posB;
    Face f = sharedFace(tail, tip);
    GC_SAFETY_ASSERT(f != Face(), "common subdivision edges must be contained in a face");
    tail = tail.inFace(f);
    tip = tip.inFace(f);

    lengthCS[e] = displacementLength(tip.faceCoords - tail.faceCoords, fEdgeLengths(f));
  }

  return lengthCS;
}

SparseMatrix<double> CommonSubdivision::vertexGalerkinMassMatrixFromPositionsA(const VertexData<Vector3>& positionsA) {
  checkMeshConstructed();
  VertexPositionGeometry geo(*mesh, interpolateAcrossA(positionsA));
  geo.requireVertexGalerkinMassMatrix();
  return geo.vertexGalerkinMassMatrix;
}

SparseMatrix<double> CommonSubdivision::vertexGalerkinMassMatrixFromPositionsB(const VertexData<Vector3>& positionsB) {
  checkMeshConstructed();
  VertexPositionGeometry geo(*mesh, interpolateAcrossB(positionsB));
  geo.requireVertexGalerkinMassMatrix();
  return geo.vertexGalerkinMassMatrix;
}

SparseMatrix<double> CommonSubdivision::vertexGalerkinMassMatrixFromLengthsA(const EdgeData<double>& lengthsA) {
  checkMeshConstructed();
  EdgeLengthGeometry geo(*mesh, interpolateEdgeLengthsA(lengthsA));
  geo.requireVertexGalerkinMassMatrix();
  return geo.vertexGalerkinMassMatrix;
}

SparseMatrix<double> CommonSubdivision::vertexGalerkinMassMatrixFromLengthsB(const EdgeData<double>& lengthsB) {
  checkMeshConstructed();
  EdgeLengthGeometry geo(*mesh, interpolateEdgeLengthsB(lengthsB));
  geo.requireVertexGalerkinMassMatrix();
  return geo.vertexGalerkinMassMatrix;
}

size_t CommonSubdivision::nVertices() const {
  size_t nV = meshB.nVertices();
  for (Edge eB : meshB.edges()) nV += intersectionsB(eB);
  return nV;
}

size_t CommonSubdivision::nEdges() const { return std::get<1>(elementCounts()); }

size_t CommonSubdivision::nFaces() const { return std::get<2>(elementCounts()); }

std::tuple<size_t, size_t, size_t> CommonSubdivision::elementCounts() const {
  size_t nV = meshB.nVertices();
  size_t nE = 0;
  size_t nF = 0;

  for (Edge eB : meshB.edges()) {
    size_t n = intersectionsB(eB);
    nV += n;
    nE += n + 1;
  }

  for (Face fB : meshB.faces()) {
    size_t nij = intersectionsB(fB.halfedge().edge());
    size_t njk = intersectionsB(fB.halfedge().next().edge());
    size_t nki = intersectionsB(fB.halfedge().next().next().edge());

    size_t ci = strictCornerCoord(njk, nki, nij);
    size_t cj = strictCornerCoord(nki, nij, njk);
    size_t ck = strictCornerCoord(nij, njk, nki);

    size_t ei = strictDegree(njk, nki, nij);
    size_t ej = strictDegree(nki, nij, njk);
    size_t ek = strictDegree(nij, njk, nki);

    nE += ci + cj + ck + ei + ej + ek;
    nF += ci + cj + ck + ei + ej + ek + 1;
  }


  return std::tuple<size_t, size_t, size_t>{nV, nE, nF};
}

size_t CommonSubdivision::intersectionsA(Edge eA) const {
  if (pointsAlongA[eA].size() == 3 && pointsAlongA[eA][1]->intersectionType == CSIntersectionType::EDGE_PARALLEL) {
    return 0; // shared edge;
  } else {
    return pointsAlongA[eA].size() - 2; // trim off endpoints
  }
}

size_t CommonSubdivision::intersectionsB(Edge eB) const {
  if (pointsAlongB[eB].size() == 3 && pointsAlongB[eB][1]->intersectionType == CSIntersectionType::EDGE_PARALLEL) {
    return 0; // shared edge;
  } else {
    return pointsAlongB[eB].size() - 2; // trim off endpoints
  }
}

std::vector<SurfacePoint> CommonSubdivision::getHalfedgePathAonB(Halfedge heA) {
  std::vector<SurfacePoint> result;

  for (CommonSubdivisionPoint* cs : pointsAlongA[heA.edge()]) {
    result.push_back(cs->posB);
  }


  if (heA != heA.edge().halfedge()) {
    std::reverse(result.begin(), result.end());
  }

  return result;
}

std::vector<SurfacePoint> CommonSubdivision::getHalfedgePathBonA(Halfedge heB) {
  std::vector<SurfacePoint> result;

  for (CommonSubdivisionPoint* cs : pointsAlongB[heB.edge()]) {
    result.push_back(cs->posA);
  }


  if (heB != heB.edge().halfedge()) {
    std::reverse(result.begin(), result.end());
  }

  return result;
}


void CommonSubdivision::constructMesh(bool triangulate, bool skipIfAlreadyConstructed) {

  if (mesh && skipIfAlreadyConstructed) {
    return;
  }

  std::vector<std::vector<size_t>> faces;
  std::vector<CommonSubdivisionPoint*> parents;
  std::vector<Face> sourceFaceA_vec;
  std::vector<Face> sourceFaceB_vec;
  constructMeshData(faces, parents, sourceFaceA_vec, sourceFaceB_vec);

  mesh.reset(new ManifoldSurfaceMesh(faces));
  sourcePoints = fromStdVector<Vertex>(*mesh, parents);
  sourceFaceA = fromStdVector<Face>(*mesh, sourceFaceA_vec);
  sourceFaceB = fromStdVector<Face>(*mesh, sourceFaceB_vec);

  if (triangulate) {
    triangulateMesh();
  }
}

void CommonSubdivision::checkMeshConstructed() const {
  if (!mesh) {
    throw std::runtime_error("common subdivision mesh has not been constructed, call constructMesh()");
  }
}

std::unique_ptr<SimplePolygonMesh> CommonSubdivision::buildSimpleMesh() {
  std::vector<std::vector<size_t>> faces;
  std::vector<CommonSubdivisionPoint*> parents;
  std::vector<Face> sourceFaceA_vec;
  std::vector<Face> sourceFaceB_vec;
  constructMeshData(faces, parents, sourceFaceA_vec, sourceFaceB_vec);

  std::vector<Vector3> dummyPositions(parents.size(), Vector3::zero());
  return std::unique_ptr<SimplePolygonMesh>(new SimplePolygonMesh(faces, dummyPositions));
}


void CommonSubdivision::constructMeshData(std::vector<std::vector<size_t>>& faces_out,
                                          std::vector<CommonSubdivisionPoint*>& parents_out,
                                          std::vector<Face>& sourceFaceA_out, std::vector<Face>& sourceFaceB_out) {

  // Compute element counts to reserve space
  // size_t nV, nE, nF;
  // std::tie(nV, nE, nF) = elementCounts();

  // Id of subdivision points in the parents list
  // TODO: if CommonSubdivisionPoint was hashable, we could just plug it
  // into the hashmap
  std::map<CommonSubdivisionPoint*, size_t> subdivisionPointsId;
  // parents_out.reserve(nE + meshB.nVertices());

  for (CommonSubdivisionPoint& p : subdivisionPoints) {
    if (p.intersectionType != CSIntersectionType::EDGE_PARALLEL) {
      subdivisionPointsId[&p] = parents_out.size();
      parents_out.push_back(&p);
    }
  }

  // Store the vertex id of each crossing along this edge
  // Considers its source vertex as the first crossing and its destination
  // vertex as the final crossing
  EdgeData<std::vector<size_t>> crossingVtxIds(meshB);

  for (Edge eB : meshB.edges()) {
    // crossingVtxIds[eB].reserve(intersectionsB(eB) + 2);

    // Source
    crossingVtxIds[eB].push_back(subdivisionPointsId[pointsAlongB[eB][0]]);

    // Middle points
    for (size_t iC = 1; iC + 1 < pointsAlongB[eB].size(); ++iC) {
      if (pointsAlongB[eB][iC]->intersectionType != CSIntersectionType::EDGE_PARALLEL) {
        crossingVtxIds[eB].push_back(subdivisionPointsId[pointsAlongB[eB][iC]]);
      }
      if (pointsAlongB[eB][iC]->intersectionType == CSIntersectionType::VERTEX_VERTEX) {
        throw std::runtime_error("encountered vertex intersection in the middle of an "
                                 "edge");
      }
    }

    // Dst
    crossingVtxIds[eB].push_back(subdivisionPointsId[pointsAlongB[eB][pointsAlongB[eB].size() - 1]]);
  }


  // faces.reserve(nF);
  // Loop over faces of mesh B and cut along edges of mesh A which cross
  // sourceFaceB_out.reserve(nF);

  bool complained = false;
  for (Face f : meshB.faces()) {
    Halfedge ij = f.halfedge();
    Halfedge jk = ij.next();
    Halfedge ki = jk.next();

    // Get list of crossings along each halfedge
    std::vector<size_t> pij = crossingVtxIds[ij.edge()];
    if (ij != ij.edge().halfedge()) std::reverse(std::begin(pij), std::end(pij));
    std::vector<size_t> pjk = crossingVtxIds[jk.edge()];
    if (jk != jk.edge().halfedge()) std::reverse(std::begin(pjk), std::end(pjk));
    std::vector<size_t> pki = crossingVtxIds[ki.edge()];
    if (ki != ki.edge().halfedge()) std::reverse(std::begin(pki), std::end(pki));

    std::vector<std::vector<size_t>> newFaces = sliceFace(pij, pjk, pki);

    for (const auto& newF : newFaces) {
      // GC_SAFETY_ASSERT(newF.size() > 2,
      //                  "No bigons allowed in common subdivision.");
      if (newF.size() <= 2) {
        if (!complained) {
          complained = true;
          std::cerr << "Error: degree-2 face in common refinement" << std::endl;
        }
        continue;
      }

      faces_out.push_back(newF);
      sourceFaceB_out.push_back(f);

      // Reconstruct the source faces on A. This uses a search around the nearby faces, which is not necessarily an
      // ideal strategy, but should work fine. ONEDAY we might want to to track some extra data about the intersections
      // to immediately recover this relationship.

      // Pick a point, which we will loop over all the neighbors of to identify the shared face.
      // Prefer a face/edge point, because they have 1/2 neighbors to test, vs. a vertex which has many.
      SurfacePoint pSearch;
      std::vector<SurfacePoint> pOthers;
      pOthers.resize(newF.size() - 1);
      for (size_t i = 0; i < newF.size(); i++) {
        SurfacePoint p = parents_out[newF[i]]->posA;
        if (i == 0) {
          pSearch = p;
        } else {
          pOthers[i - 1] = p;
        }

        if (pOthers.back().type == SurfacePointType::Edge && pSearch.type == SurfacePointType::Vertex) {
          std::swap(pOthers.back(), pSearch);
        }
        if (pOthers.back().type == SurfacePointType::Face && pSearch.type != SurfacePointType::Face) {
          std::swap(pOthers.back(), pSearch);
        }
      }


      // Search over the neighboring faces to the point, looking for a shared face
      Face sharedFace;

      // Assemble a list of neighbors to test
      std::vector<Face> testFaces;
      switch (pSearch.type) {
      case SurfacePointType::Vertex:
        for (Face f : pSearch.vertex.adjacentFaces()) {
          testFaces.push_back(f);
        }
        break;
      case SurfacePointType::Edge:
        testFaces.push_back(pSearch.edge.halfedge().face());
        testFaces.push_back(pSearch.edge.halfedge().twin().face());
        break;
      case SurfacePointType::Face:
        testFaces.push_back(pSearch.face);
        break;
      }

      // Test the neibhbors to see if any works
      for (Face f : testFaces) {
        bool isGood = true;
        SurfacePoint testPt(f, Vector3::zero());
        for (SurfacePoint& pO : pOthers) {
          if (!checkAdjacent(testPt, pO)) {
            isGood = false;
            break;
          }
        }
        if (isGood) {
          sharedFace = f;
          break;
        }
      }

      GC_SAFETY_ASSERT(sharedFace != Face(), "could not identify source face on mesh A")
      sourceFaceA_out.push_back(sharedFace);
    }
  }
}


void CommonSubdivision::triangulateMesh() {
  checkMeshConstructed();
  for (Face oldFace : mesh->faces()) {
    std::vector<Face> newFaces = mesh->triangulate(oldFace);

    // fix up data arrays
    for (Face newFace : newFaces) {
      sourceFaceA[newFace] = sourceFaceA[oldFace];
      sourceFaceB[newFace] = sourceFaceB[oldFace];
    }
  }

  mesh->compress();
}

std::vector<std::vector<size_t>> sliceFace(const std::vector<size_t>& pij, const std::vector<size_t>& pjk,
                                           const std::vector<size_t>& pki) {

  if (pij.size() >= pjk.size() && pij.size() >= pki.size()) {
    return sliceNicelyOrderedFace(pij, pjk, pki);
  } else if (pjk.size() >= pki.size() && pjk.size() >= pij.size()) {
    return sliceNicelyOrderedFace(pjk, pki, pij);
  } else {
    return sliceNicelyOrderedFace(pki, pij, pjk);
  }
}


// Precondition: pij.size() >= pjk.size(), pki.size()
std::vector<std::vector<size_t>> sliceNicelyOrderedFace(const std::vector<size_t>& pij, const std::vector<size_t>& pjk,
                                                        const std::vector<size_t>& pki) {
  size_t nij = pij.size() - 2;
  size_t njk = pjk.size() - 2;
  size_t nki = pki.size() - 2;

  // WATCH3(nij, njk, nki);

  auto pji = [&](size_t a) { return pij[pij.size() - 1 - a]; };
  auto pkj = [&](size_t a) { return pjk[pjk.size() - 1 - a]; };
  auto pik = [&](size_t a) { return pki[pki.size() - 1 - a]; };

  auto dropAdjacentDuplicates = [](std::vector<size_t> vec) -> std::vector<size_t> {
    // Weird for loop to erase elements as we go
    // https://stackoverflow.com/a/31329841
    for (size_t i = 0; i < vec.size();) {
      if (vec[i] == vec[(i + 1) % vec.size()]) {
        vec.erase(vec.begin() + i);
      } else {
        i++;
      }
    }
    return vec;
  };

  std::vector<std::vector<size_t>> faces;
  if (nij > njk + nki) {
    // cout << "Fan!" << endl;
    // Fan configuration where vertex k has edges coming out of it

    // Corner coordinates
    // ek must be >= 1, otherwise we're in triforce
    size_t ci = nki;
    size_t cj = njk;
    size_t ek = nij - njk - nki;

    // WATCH3(ci, cj, ek);

    // Make corner faces
    // These go one farther than the corner faces in the triforce case
    // since the go right up to the first edge emanating from k
    for (size_t iC = 0; iC <= ci; ++iC) {
      faces.push_back(dropAdjacentDuplicates({pij[iC], pij[iC + 1], pik(iC + 1), pik(iC)}));
    }
    for (size_t iC = 0; iC <= cj; ++iC) {
      faces.push_back(dropAdjacentDuplicates({pjk[iC], pjk[iC + 1], pji(iC + 1), pji(iC)}));
    }

    // Make central triangles
    for (size_t iE = 0; iE + 1 < ek; iE++) {
      faces.push_back(dropAdjacentDuplicates({pki[0], pij[ci + 1 + iE], pij[ci + 1 + iE + 1]}));
    }

  } else {
    // cout << "Triforce!" << endl;
    // Triforce configuration
    GC_SAFETY_ASSERT((nij + njk + nki) % 2 == 0, "normal coordinates which obey the triangle "
                                                 "inequality must sum to an even number");

    // Corner coordinates
    size_t ci = (nij - njk + nki) / 2;
    size_t cj = (njk - nki + nij) / 2;
    size_t ck = (nki - nij + njk) / 2;

    // cout << "ci: " << ci << "\tcj: " << cj << "\tck: " << ck << endl;

    // Make corner faces
    for (size_t iC = 0; iC + 1 <= ci; ++iC) {
      faces.push_back(dropAdjacentDuplicates({pij[iC], pij[iC + 1], pik(iC + 1), pik(iC)}));
    }
    for (size_t iC = 0; iC + 1 <= cj; ++iC) {
      faces.push_back(dropAdjacentDuplicates({pjk[iC], pjk[iC + 1], pji(iC + 1), pji(iC)}));
      // WATCH(iC);
      // cout << faces[faces.size() - 1] << endl;
    }
    for (size_t iC = 0; iC + 1 <= ck; ++iC) {
      faces.push_back(dropAdjacentDuplicates({pki[iC], pki[iC + 1], pkj(iC + 1), pkj(iC)}));
    }

    // Make central hexagon
    faces.push_back(dropAdjacentDuplicates({pik(ci), pij[ci], pji(cj), pjk[cj], pkj(ck), pki[ck]}));
  }

  return faces;
}

// Write mesh A and common subdivision to obj files
// Vertex positions should be for mesh A
void CommonSubdivision::writeToFile(std::string filename, const VertexData<Vector3>& vertexPositions, int kColors) {
  checkMeshConstructed();

  VertexData<Vector3> subdivisionPositions = interpolateAcrossA(vertexPositions);

  FaceData<double> faceColors = niceColors(meshB, kColors);
  CornerData<Vector2> texCoords(*mesh);
  for (Corner c : mesh->corners()) texCoords[c] = Vector2{faceColors[sourceFaceB[c.face()]], 0.5};

  VertexPositionGeometry inputGeo(meshA, vertexPositions);
  writeSurfaceMesh(meshA, inputGeo, filename + "_input.obj", "obj");

  VertexPositionGeometry subdivisionGeo(*mesh, subdivisionPositions);
  writeSurfaceMesh(*mesh, subdivisionGeo, texCoords, filename + "_intrinsic.obj", "obj");
}

FaceData<double> niceColors(ManifoldSurfaceMesh& mesh, int kColors) {
  // Evenly spaced values on [0,1], avoiding extrema
  std::vector<double> kVals = {1.0 / (2.0 * kColors)};
  for (int i = 1; i < kColors; i++) {
    kVals.push_back(kVals[i - 1] + 1.0 / kColors);
  }
  std::vector<size_t> lastUsedT(kColors, 0);
  size_t timestamp = 0;


  FaceData<double> faceColor(mesh, -1.);
  for (Face f : mesh.faces()) {

    int bestColor = -1;
    size_t oldestColorT = std::numeric_limits<size_t>::max();
    int bestAdjacencyScore = -1; // 0 == used by a neighbor, 1 == used by a neighbor's neighbor

    // Consider all colors
    for (int iCandColor = 0; iCandColor < kColors; iCandColor++) {

      // Check against neighbors
      int adjacencyScore = 2;
      for (Face fN : f.adjacentFaces()) {
        // Single neighbors
        if (faceColor[fN] == kVals[iCandColor]) adjacencyScore = 0;

        // Double neighbors
        for (Face fNN : fN.adjacentFaces()) {
          if (faceColor[fNN] == kVals[iCandColor]) adjacencyScore = std::min(adjacencyScore, 1);
        }
      }

      // Accept if less adjacent then previous best adjacent score
      if (adjacencyScore > bestAdjacencyScore ||
          (adjacencyScore == bestAdjacencyScore && lastUsedT[iCandColor] < oldestColorT)) {
        oldestColorT = lastUsedT[iCandColor];
        bestAdjacencyScore = adjacencyScore;
        bestColor = iCandColor;
      }
    }

    // Should always find something if k > 3
    if (bestColor == -1) throw std::runtime_error("failed to color");

    faceColor[f] = kVals[bestColor];
    lastUsedT[bestColor] = timestamp;

    timestamp++;
  }

  return faceColor;
}


} // namespace surface
} // namespace geometrycentral
