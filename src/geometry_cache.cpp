#include "geometrycentral/geometry_cache.h"

#include "geometrycentral/geometry.h"
#include "geometrycentral/discrete_operators.h"

using std::cout; using std::endl;

namespace geometrycentral {

void DependentQuantity::ensureHaveIfRequired() {
  if(requireCount > 0) {
    ensureHave();
  }
}

void DependentQuantity::ensureHave() {

  // If the quantity is already populated, early out
  if (computed) {
    return;
  }

  // Resolve all of the dependencies
  for (auto x : dependencies) {
    x->ensureHave();
  }

  // Compute this quantity
  evaluateFunc();

  computed = true;
};

void DependentQuantity::require() {
  requireCount++;
  ensureHave();
}

void DependentQuantity::unrequire() {
  requireCount--;

  if (requireCount < 0) {
    throw std::logic_error("Quantity was unrequire()'d more than than it was require()'d");
    requireCount = 0;
  }
}

template <typename G>
GeometryCache<G>::GeometryCache(Geometry<G>* geometry_) : geometry(geometry_) {
  mesh = geometry->getMesh();

  // Helper to add a quantity, binding this instance to its compute function and adding it to the list of all quantities
  auto addQuantity = [&](DependentQuantity& q, std::vector<DependentQuantity*> deps, std::function<void(GeometryCache<G>*)> func) { 
    q = DependentQuantity(deps, std::bind(func, this));
    allQuantities.push_back(&q); 
  };

  // === ALL the quantities
  //          quantity manager            dependencies                                  compute function
  
  // == Basic geometric quantities
  addQuantity(faceAreaNormalsQ,           {},                                           &GeometryCache<G>::computeFaceAreaNormals);
  addQuantity(faceAreasQ,                 {&faceAreaNormalsQ},                          &GeometryCache<G>::computeFaceAreas);
  addQuantity(faceNormalsQ,               {&faceAreaNormalsQ},                          &GeometryCache<G>::computeFaceNormals);
  addQuantity(vertexNormalsQ,             {&faceAreaNormalsQ},                          &GeometryCache<G>::computeVertexNormals);
  addQuantity(vertexDualAreasQ,           {&faceAreasQ},                                &GeometryCache<G>::computeVertexDualAreas);
  addQuantity(halfedgeVectorsQ,           {},                                           &GeometryCache<G>::computeHalfedgeVectors);
  addQuantity(edgeLengthsQ,               {},                                           &GeometryCache<G>::computeEdgeLengths);
  addQuantity(dihedralAnglesQ,            {},                                           &GeometryCache<G>::computeDihedralAngles);
  addQuantity(halfedgeCotanWeightsQ,      {},                                           &GeometryCache<G>::computeHalfedgeCotanWeights);
  addQuantity(edgeCotanWeightsQ,          {},                                           &GeometryCache<G>::computeEdgeCotanWeights);
  addQuantity(vertexAngleDefectsQ,        {&halfedgeOppositeAnglesQ},                   &GeometryCache<G>::computeVertexAngleDefects);
  
  // == Vector fields, angles, and transport
  addQuantity(faceBasesQ,                       {&faceNormalsQ},                                            &GeometryCache<G>::computeFaceBases);
  addQuantity(vertexBasesQ,                     {&vertexNormalsQ},                                          &GeometryCache<G>::computeVertexBases);
  addQuantity(halfedgeFaceCoordsQ,              {&faceBasesQ},                                              &GeometryCache<G>::computeHalfedgeFaceCoords);
  addQuantity(faceTransportCoefsQ,              {&faceBasesQ, &halfedgeFaceCoordsQ},                        &GeometryCache<G>::computeFaceTransportCoefs);
  addQuantity(halfedgeOppositeAnglesQ,          {},                                                         &GeometryCache<G>::computeHalfedgeOppositeAngles);
  addQuantity(halfedgeRescaledOppositeAnglesQ,  {&vertexAngleDefectsQ, &halfedgeOppositeAnglesQ},           &GeometryCache<G>::computeHalfedgeRescaledOppositeAngles);
  addQuantity(halfedgeVertexCoordsQ,            {&halfedgeRescaledOppositeAnglesQ},                         &GeometryCache<G>::computeHalfedgeVertexCoords);
  addQuantity(vertexTransportCoefsQ,            {&vertexBasesQ, &halfedgeVertexCoordsQ},                    &GeometryCache<G>::computeVertexTransportCoefs);
  addQuantity(vertexFaceTransportCoefsQ,        {&halfedgeFaceCoordsQ, &halfedgeVertexCoordsQ},             &GeometryCache<G>::computeVertexFaceTransportCoefs);
  addQuantity(principalDirectionsQ,             {&edgeLengthsQ, &dihedralAnglesQ, &halfedgeVertexCoordsQ},  &GeometryCache<G>::computePrincipalDirections);
  
  // == Indices
  addQuantity(vertexIndicesQ,             {},                                           &GeometryCache<G>::computeVertexIndices);
  addQuantity(interiorVertexIndicesQ,     {},                                           &GeometryCache<G>::computeInteriorVertexIndices);
  addQuantity(faceIndicesQ,               {},                                           &GeometryCache<G>::computeFaceIndices);
  addQuantity(edgeIndicesQ,               {},                                           &GeometryCache<G>::computeEdgeIndices);
  addQuantity(halfedgeIndicesQ,           {},                                           &GeometryCache<G>::computeHalfedgeIndices);

  // == Operators
  addQuantity(basicDECOperatorsQ,         {&vertexDualAreasQ, &edgeCotanWeightsQ},      &GeometryCache<G>::computeBasicDECOperators);
  addQuantity(modifiedDECOperatorsQ,      {&basicDECOperatorsQ},                        &GeometryCache<G>::computeModifiedDECOperators);
  addQuantity(zeroFormWeakLaplacianQ,     {&basicDECOperatorsQ},                        &GeometryCache<G>::computeZeroFormWeakLaplacian);
}

template <typename G>
void GeometryCache<G>::repopulate() {

  for(DependentQuantity* q : allQuantities) {
    q->computed = false;
  }
  for(DependentQuantity* q : allQuantities) {
    q->ensureHaveIfRequired();
  }

}

// Helper
namespace {
void verifyTriangular(HalfedgeMesh* m) {
  if (!m->isSimplicial()) {
    throw std::logic_error("Only implemented for triangular meshes");
  }
}
} // namespace

// === Quantity implementations

// Specialization for euclidean
template <>
void GeometryCache<Euclidean>::computeFaceAreaNormals() {
  faceAreaNormals = FaceData<Vector3>(mesh);
  for (FacePtr f : mesh->faces()) {
    Vector3 AN{0., 0., 0.};
    for (HalfedgePtr h : f.adjacentHalfedges()) {
      Vector3 pi = geometry->position(h.vertex());
      Vector3 pj = geometry->position(h.twin().vertex());
      AN += cross(pi, pj);
    }
    faceAreaNormals[f] = AN * 0.5;
  }

}

template <>
void GeometryCache<Euclidean>::computeFaceAreas() {
  faceAreas = FaceData<double>(mesh);
  for (FacePtr f : mesh->faces()) {
    faceAreas[f] = norm(faceAreaNormals[f]);
  }
}

template <>
void GeometryCache<Euclidean>::computeFaceNormals() {
  faceNormals = FaceData<Vector3>(mesh);
  for (FacePtr f : mesh->faces()) {
    faceNormals[f] = unit(faceAreaNormals[f]);
  }
}

// Area-weighted vertex normals
template <>
void GeometryCache<Euclidean>::computeVertexNormals() {
  vertexNormals = VertexData<Vector3>(mesh);
  for (VertexPtr v : mesh->vertices()) {
    Vector3 n{0., 0., 0.};
    for (FacePtr f : v.adjacentFaces()) {
      n += faceAreaNormals[f];
    }
    vertexNormals[v] = unit(n);
  }
}

template <>
void GeometryCache<Euclidean>::computeVertexDualAreas() {
  vertexDualAreas = VertexData<double>(mesh);
  for (VertexPtr v : mesh->vertices()) {
    double A = 0;
    for (FacePtr f : v.adjacentFaces()) {
      A += faceAreas[f];
    }
    vertexDualAreas[v] = A / 3.0;
  }
}

template <>
void GeometryCache<Euclidean>::computeHalfedgeVectors() {
  halfedgeVectors = HalfedgeData<Vector3>(mesh);
  for(HalfedgePtr he : mesh->halfedges()) {
    halfedgeVectors[he] = geometry->vector(he);
  }
}

template <>
void GeometryCache<Euclidean>::computeEdgeLengths() {
  edgeLengths = EdgeData<double>(mesh);
  for(EdgePtr e : mesh->edges()) {
    edgeLengths[e] = norm(geometry->vector(e.halfedge()));
  }
}

template <>
void GeometryCache<Euclidean>::computeFaceBases() {
  faceBases = FaceData<std::array<Vector3,2>>(mesh);
  for(FacePtr f : mesh->faces()) {
    Vector3 basisX = unit(geometry->vector(f.halfedge()));
    Vector3 basisY = basisX.rotate_around(faceNormals[f], PI/2.0);
    faceBases[f][0] = basisX;
    faceBases[f][1] = basisY;
  }
}

template <>
void GeometryCache<Euclidean>::computeVertexBases() {
  vertexBases = VertexData<std::array<Vector3,2>>(mesh);

  for(VertexPtr v : mesh->vertices()) {

    // Chose any edge (projected in to tangent plane) as one component of basis
    Vector3 basisX = unit(geometry->vector(v.halfedge()));
    basisX = unit(basisX * (1.0 - dot(vertexNormals[v], basisX)));
    Vector3 basisY = basisX.rotate_around(vertexNormals[v], PI/2.0);

    vertexBases[v][0] = basisX;
    vertexBases[v][1] = basisY;
  }
}

template <>
void GeometryCache<Euclidean>::computeHalfedgeFaceCoords() {
  halfedgeFaceCoords = HalfedgeData<Complex>(mesh);

  for (HalfedgePtr he : mesh->halfedges()) {
    FacePtr f = he.face();

    if (f.isReal()) {
      Vector3 heVec = geometry->vector(he);
      halfedgeFaceCoords[he] = Complex(dot(faceBases[f][0], heVec), dot(faceBases[f][1], heVec));
    } else {
      halfedgeFaceCoords[he] =
          std::numeric_limits<double>::quiet_NaN(); // using this basis is never a good idea, so NaN-out
    }
  }
}



template <>
void GeometryCache<Euclidean>::computeFaceTransportCoefs() {

  faceTransportCoefs = HalfedgeData<Complex>(mesh);

  for (HalfedgePtr he : mesh->halfedges()) {
    Complex angleInSource = halfedgeFaceCoords[he];
    Complex desiredAngleInTarget = -halfedgeFaceCoords[he.twin()];
    faceTransportCoefs[he] = desiredAngleInTarget / angleInSource;
  }
}

// Note: This is not the same discretization used in Globally Optimal Direction Fields (etc), due to the difference in 
template <>
void GeometryCache<Euclidean>::computeVertexTransportCoefs() {

  vertexTransportCoefs = HalfedgeData<Complex>(mesh);

  for (HalfedgePtr he : mesh->halfedges()) {
    Complex angleInSource = halfedgeVertexCoords[he];
    Complex desiredAngleInTarget = -halfedgeVertexCoords[he.twin()];
    vertexTransportCoefs[he] = desiredAngleInTarget / angleInSource;
  }
}

template <>
void GeometryCache<Euclidean>::computeVertexFaceTransportCoefs() {

  vertexFaceTransportCoefs = HalfedgeData<Complex>(mesh);

  for (HalfedgePtr he : mesh->halfedges()) {
    Complex angleInSource = halfedgeVertexCoords[he];
    Complex desiredAngleInTarget = halfedgeFaceCoords[he];
    vertexFaceTransportCoefs[he] = desiredAngleInTarget / angleInSource;
  }
}

template <>
void GeometryCache<Euclidean>::computePrincipalDirections() {

  principalDirections = VertexData<Complex>(mesh);

  for (VertexPtr v : mesh->vertices()) {
    Complex principalDir{0.0, 0.0};
    for (HalfedgePtr he : v.outgoingHalfedges()) {
      double len = edgeLengths[he.edge()];
      double alpha = dihedralAngles[he.edge()];
      Complex vec = halfedgeVertexCoords[he];
      principalDir += -vec * vec / len * std::abs(alpha);
    }
    principalDirections[v] = principalDir / 4.0;
  }
}

template <>
void GeometryCache<Euclidean>::computeDihedralAngles() {
  dihedralAngles = EdgeData<double>(mesh);
  for (EdgePtr e : mesh->edges()) {
    dihedralAngles[e] = geometry->dihedralAngle(e);
  }
}

template <>
void GeometryCache<Euclidean>::computeHalfedgeOppositeAngles() {
  verifyTriangular(mesh);

  halfedgeOppositeAngles = HalfedgeData<double>(mesh);
  for (HalfedgePtr he : mesh->halfedges()) {
    if (he.isReal()) {
      halfedgeOppositeAngles[he] = angle(-geometry->vector(he.next()), geometry->vector(he.next().next()));
    } else {
      halfedgeOppositeAngles[he] = std::numeric_limits<double>::quiet_NaN();
    }
  }
}

template <>
void GeometryCache<Euclidean>::computeHalfedgeCotanWeights() {
  halfedgeCotanWeights = HalfedgeData<double>(mesh);
  for (HalfedgePtr he : mesh->halfedges()) {
    halfedgeCotanWeights[he] = geometry->cotan(he);
  }
}

template <>
void GeometryCache<Euclidean>::computeEdgeCotanWeights() {
  edgeCotanWeights = EdgeData<double>(mesh);
  for (EdgePtr e : mesh->edges()) {
    edgeCotanWeights[e] = geometry->cotanWeight(e);
  }
}

template <>
void GeometryCache<Euclidean>::computeVertexAngleDefects() {
  verifyTriangular(mesh);

  vertexAngleDefects = VertexData<double>(mesh);
  for (VertexPtr v : mesh->vertices()) {
    if (v.isBoundary()) {
      // Convention: no curvature at boundary
      vertexAngleDefects[v] = 0;
    } else {
      double sum = 0;
      for (HalfedgePtr he : v.outgoingHalfedges()) {
        sum += halfedgeOppositeAngles[he.next()];
      }
      vertexAngleDefects[v] = 2. * PI - sum;
    }
  }
}

template <>
void GeometryCache<Euclidean>::computeHalfedgeRescaledOppositeAngles() {
  halfedgeRescaledOppositeAngles = HalfedgeData<double>(mesh);
  for (HalfedgePtr he : mesh->halfedges()) {
    double origSum = 2. * PI - vertexAngleDefects[he.vertex()];
    halfedgeRescaledOppositeAngles[he] = halfedgeOppositeAngles[he] * 2. * PI / origSum;
  }
}

template <>
void GeometryCache<Euclidean>::computeHalfedgeVertexCoords() {
  verifyTriangular(mesh);

  halfedgeVertexCoords = HalfedgeData<Complex>(mesh);

  for (VertexPtr v : mesh->vertices()) {

    if (v.isBoundary()) {

      // First, check what angle we associated with the boundary wedge
      // (recall that in a manifold triangle mesh, there can be at most one boundary wedge)
      double angleSum = 0;
      HalfedgePtr afterBoundaryHe;
      for (HalfedgePtr he : v.outgoingHalfedges()) {
        if (he.isReal()) {
          angleSum += halfedgeRescaledOppositeAngles[he.next()];
        }
        if (!he.twin().isReal()) {
          afterBoundaryHe = he;
        }
      }
      double boundaryAngle = 2 * PI - angleSum;

      // Now, loop like in the usual case, but substitute the boundary value when needed
      double coordSum = 0.0;

      // Custom loop to orbit CCW
      HalfedgePtr firstHe = v.halfedge();
      HalfedgePtr currHe = firstHe;
      do {
        halfedgeVertexCoords[currHe] = std::exp(coordSum * IM_I);
        if (currHe.isReal()) {
          coordSum += halfedgeRescaledOppositeAngles[currHe.next()];
          currHe = currHe.next().next().twin();
        } else {
          coordSum += boundaryAngle;
          currHe = afterBoundaryHe;
        }
      } while (currHe != firstHe);


    } else {
      double coordSum = 0.0;

      // Custom loop to orbit CCW
      HalfedgePtr firstHe = v.halfedge();
      HalfedgePtr currHe = firstHe;
      do {
        halfedgeVertexCoords[currHe] = std::exp(coordSum * IM_I);
        coordSum += halfedgeRescaledOppositeAngles[currHe.next()];
        currHe = currHe.next().next().twin();
      } while (currHe != firstHe);
    }
  }
}

template <typename T>
void GeometryCache<T>::computeVertexIndices() {
  vertexIndices = mesh->getVertexIndices();
}

template <typename T>
void GeometryCache<T>::computeInteriorVertexIndices() {
  interiorVertexIndices = mesh->getInteriorVertexIndices();
}

template <typename T>
void GeometryCache<T>::computeFaceIndices() {
  faceIndices = mesh->getFaceIndices();
}

template <typename T>
void GeometryCache<T>::computeEdgeIndices() {
  edgeIndices = mesh->getEdgeIndices();
}


template <typename T>
void GeometryCache<T>::computeHalfedgeIndices() {
  halfedgeIndices = mesh->getHalfedgeIndices();
}

template <>
void GeometryCache<Euclidean>::computeBasicDECOperators() {
  d0 = buildDerivative0(mesh);
  d1 = buildDerivative1(mesh);
  hodge0 = buildHodge0(geometry);
  hodge1 = buildHodge1(geometry);
  hodge2 = buildHodge2(geometry);
}

template <>
void GeometryCache<Euclidean>::computeModifiedDECOperators() {
  hodge0Inv = buildHodge0(geometry).inverse();
  hodge1Inv = buildHodge1(geometry).inverse();
  hodge2Inv = buildHodge2(geometry).inverse();
}


template <>
void GeometryCache<Euclidean>::computeZeroFormWeakLaplacian() {
  zeroFormWeakLaplacian = d0.transpose() * hodge1 * d0;
}


// Explicit template instantions
// Note: inherits problems with Geometry<T>, only works for Euclidean
template class GeometryCache<Euclidean>;
// template class GeometryCache<Planar>;
// template class GeometryCache<Spherical>;


} // namespace geometrycentral
