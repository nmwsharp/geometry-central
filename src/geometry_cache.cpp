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

  // ALL the quantities
  addQuantity(faceAreaNormalsQ,       {},                                 &GeometryCache<G>::computeFaceAreaNormals);
  addQuantity(faceAreasQ,             {&faceAreaNormalsQ},                &GeometryCache<G>::computeFaceAreas);
  addQuantity(faceNormalsQ,           {&faceAreaNormalsQ},                &GeometryCache<G>::computeFaceNormals);
  addQuantity(vertexNormalsQ,         {&faceAreaNormalsQ},                &GeometryCache<G>::computeVertexNormals);
  addQuantity(vertexDualAreasQ,       {&faceAreasQ},                      &GeometryCache<G>::computeVertexDualAreas);
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

// Explicit template instantions
// Note: inherits problems with Geometry<T>, only works for Euclidean
template class GeometryCache<Euclidean>;
// template class GeometryCache<Planar>;
// template class GeometryCache<Spherical>;


} // namespace geometrycentral
