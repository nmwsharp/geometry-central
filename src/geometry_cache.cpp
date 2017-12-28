#include "geometrycentral/geometry_cache.h"

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

GeometryCache::GeometryCache(Geometry<Euclidean>* geometry_) : geometry(geometry_) {
  mesh = geometry->getMesh();

  // Helper to add a quantity, binding this instance to its compute function and adding it to the list of all quantities
  auto addQuantity = [&](DependentQuantity& q, std::vector<DependentQuantity*> deps, std::function<void(GeometryCache*)> func) { 
    q = DependentQuantity(deps, std::bind(func, this));
    allQuantities.push_back(&q); 
  };

  // ALL the quantities
  addQuantity(faceAreaNormalsQ,       {},                                 &GeometryCache::computeFaceAreaNormals);
  addQuantity(faceAreasQ,             {&faceAreaNormalsQ},                &GeometryCache::computeFaceAreas);
  addQuantity(faceNormalsQ,           {&faceAreaNormalsQ},                &GeometryCache::computeFaceNormals);
  addQuantity(vertexNormalsQ,         {&faceAreaNormalsQ},                &GeometryCache::computeVertexNormals);
  addQuantity(vertexDualAreasQ,       {&faceAreasQ},                      &GeometryCache::computeVertexDualAreas);
}

void GeometryCache::repopulate() {

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

void GeometryCache::computeFaceAreaNormals() {
  faceAreaNormalsRaw = FaceData<Vector3>(mesh);
  for (FacePtr f : mesh->faces()) {
    Vector3 AN{0., 0., 0.};
    for (HalfedgePtr h : f.adjacentHalfedges()) {
      Vector3 pi = geometry->position(h.vertex());
      Vector3 pj = geometry->position(h.twin().vertex());
      AN += cross(pi, pj);
    }
    faceAreaNormalsRaw[f] = AN * 0.5;
  }
}

void GeometryCache::computeFaceAreas() {
  faceAreasRaw = FaceData<double>(mesh);
  for (FacePtr f : mesh->faces()) {
    faceAreasRaw[f] = norm(faceAreaNormals[f]);
  }
}

void GeometryCache::computeFaceNormals() {
  faceNormalsRaw = FaceData<Vector3>(mesh);
  for (FacePtr f : mesh->faces()) {
    faceNormalsRaw[f] = unit(faceAreaNormals[f]);
  }
}

// Area-weighted vertex normals
void GeometryCache::computeVertexNormals() {
  vertexNormalsRaw = VertexData<Vector3>(mesh);
  for (VertexPtr v : mesh->vertices()) {
    Vector3 n{0., 0., 0.};
    for (FacePtr f : v.adjacentFaces()) {
      n += faceAreaNormals[f];
    }
    vertexNormalsRaw[v] = unit(n);
  }
}

void GeometryCache::computeVertexDualAreas() {
  vertexDualAreasRaw = VertexData<double>(mesh);
  for (VertexPtr v : mesh->vertices()) {
    double A = 0;
    for (FacePtr f : v.adjacentFaces()) {
      A += faceAreas[f];
    }
    vertexDualAreasRaw[v] = A / 3.0;
  }
}



} // namespace geometrycentral
