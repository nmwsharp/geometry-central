#include "geometrycentral/geometry_cache.h"

#include "geometrycentral/geometry.h"

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
  //          quantity manager        dependencies                        compute function
  addQuantity(faceAreaNormalsQ,       {},                                 &GeometryCache<G>::computeFaceAreaNormals);
  addQuantity(faceAreasQ,             {&faceAreaNormalsQ},                &GeometryCache<G>::computeFaceAreas);
  addQuantity(faceNormalsQ,           {&faceAreaNormalsQ},                &GeometryCache<G>::computeFaceNormals);
  addQuantity(vertexNormalsQ,         {&faceAreaNormalsQ},                &GeometryCache<G>::computeVertexNormals);
  addQuantity(vertexDualAreasQ,       {&faceAreasQ},                      &GeometryCache<G>::computeVertexDualAreas);
  addQuantity(edgeLengthsQ,           {},                                 &GeometryCache<G>::computeEdgeLengths);
  addQuantity(faceBasisQ,             {&faceNormalsQ},                    &GeometryCache<G>::computeFaceBasis);
  addQuantity(faceTransportCoefsQ,    {&faceNormalsQ, &faceBasisQ},       &GeometryCache<G>::computeFaceTransportCoefs);
  addQuantity(dihedralAngleQ,         {},                                 &GeometryCache<G>::computeDihedralAngle);
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
void checkTriangular(HalfedgeMesh* m) {
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
void GeometryCache<Euclidean>::computeEdgeLengths() {
  edgeLengths = EdgeData<double>(mesh);
  for(EdgePtr e : mesh->edges()) {
    edgeLengths[e] = norm(geometry->vector(e.halfedge()));
  }
}

template <>
void GeometryCache<Euclidean>::computeFaceBasis() {
  faceBasis = FaceData<std::array<Vector3,2>>(mesh);
  for(FacePtr f : mesh->faces()) {
    Vector3 basisX = unit(geometry->vector(f.halfedge()));
    Vector3 basisY = basisX.rotate_around(faceNormals[f], PI/2.0);
    faceBasis[f][0] = basisX;
    faceBasis[f][1] = basisY;
  }
}

template <>
void GeometryCache<Euclidean>::computeFaceTransportCoefs() {

  faceTransportCoefs = HalfedgeData<Complex>(mesh);

  for(HalfedgePtr he : mesh->halfedges()) {

    // Boundary halfedges don't make sense
    if(!he.twin().isReal() || !he.isReal()) {
      faceTransportCoefs[he] = 0;
      continue;
    }

    FacePtr sourceFace = he.face();
    FacePtr targetFace = he.twin().face();

    // Rotate the neighboring face's basis vector to this one
    Vector3 sourceBasis = faceBasis[sourceFace][0];
    Vector3 targetBasis = faceBasis[targetFace][0];
    Vector3 edgeAxis = unit(geometry->vector(he));
    double edgeAngle = -geometry->dihedralAngle(he.edge());
    Vector3 sourceBasisInFace = sourceBasis.rotate_around(edgeAxis, edgeAngle);

    // Measure the angle between the two vectors now that they are in plane
    // Recall that the transform for data is the opposite of the transform for the axis, hence the negative sign
    double angle = -angleInPlane(sourceBasisInFace, targetBasis, faceNormals[targetFace]);

    // LC connection between the faces
    Complex rBar = std::exp(angle * IM_I);

    faceTransportCoefs[he] = rBar;
  }
}

template <>
void GeometryCache<Euclidean>::computeDihedralAngle() {
  dihedralAngle = EdgeData<double>(mesh);
  for(EdgePtr e : mesh->edges()) {
    dihedralAngle[e] = geometry->dihedralAngle(e);
  }
}



// Explicit template instantions
// Note: inherits problems with Geometry<T>, only works for Euclidean
template class GeometryCache<Euclidean>;
// template class GeometryCache<Planar>;
// template class GeometryCache<Spherical>;


} // namespace geometrycentral
