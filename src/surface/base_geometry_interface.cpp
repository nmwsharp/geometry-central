#include "geometrycentral/surface/base_geometry_interface.h"

namespace geometrycentral {
namespace surface {


// clang-format off
BaseGeometryInterface::BaseGeometryInterface(SurfaceMesh& mesh_)
    : mesh(mesh_),
      
  // Construct the dependency graph of managed quantities and their callbacks

  vertexIndicesQ           (&vertexIndices,         std::bind(&BaseGeometryInterface::computeVertexIndices, this),          quantities),
  interiorVertexIndicesQ   (&interiorVertexIndices, std::bind(&BaseGeometryInterface::computeInteriorVertexIndices, this),  quantities),
  edgeIndicesQ             (&edgeIndices,           std::bind(&BaseGeometryInterface::computeEdgeIndices, this),            quantities),
  halfedgeIndicesQ         (&halfedgeIndices,       std::bind(&BaseGeometryInterface::computeHalfedgeIndices, this),        quantities),
  cornerIndicesQ           (&cornerIndices,         std::bind(&BaseGeometryInterface::computeCornerIndices, this),          quantities),
  faceIndicesQ             (&faceIndices,           std::bind(&BaseGeometryInterface::computeFaceIndices, this),            quantities),
  boundaryLoopIndicesQ     (&boundaryLoopIndices,   std::bind(&BaseGeometryInterface::computeBoundaryLoopIndices, this),    quantities)

  {
  }
// clang-format on

BaseGeometryInterface::~BaseGeometryInterface() {}

void BaseGeometryInterface::refreshQuantities() {
  for (DependentQuantity* q : quantities) {
    q->computed = false;
  }
  for (DependentQuantity* q : quantities) {
    q->ensureHaveIfRequired();
  }
}

void BaseGeometryInterface::purgeQuantities() {
  for (DependentQuantity* q : quantities) {
    q->clearIfNotRequired();
  }
}

// == Indices

// Vertex indices
void BaseGeometryInterface::computeVertexIndices() { vertexIndices = mesh.getVertexIndices(); }
void BaseGeometryInterface::requireVertexIndices() { vertexIndicesQ.require(); }
void BaseGeometryInterface::unrequireVertexIndices() { vertexIndicesQ.unrequire(); }

// Interior vertex indices
void BaseGeometryInterface::computeInteriorVertexIndices() { interiorVertexIndices = mesh.getInteriorVertexIndices(); }
void BaseGeometryInterface::requireInteriorVertexIndices() { interiorVertexIndicesQ.require(); }
void BaseGeometryInterface::unrequireInteriorVertexIndices() { interiorVertexIndicesQ.unrequire(); }

// Edge indices
void BaseGeometryInterface::computeEdgeIndices() { edgeIndices = mesh.getEdgeIndices(); }
void BaseGeometryInterface::requireEdgeIndices() { edgeIndicesQ.require(); }
void BaseGeometryInterface::unrequireEdgeIndices() { edgeIndicesQ.unrequire(); }

// Halfedge indices
void BaseGeometryInterface::computeHalfedgeIndices() { halfedgeIndices = mesh.getHalfedgeIndices(); }
void BaseGeometryInterface::requireHalfedgeIndices() { halfedgeIndicesQ.require(); }
void BaseGeometryInterface::unrequireHalfedgeIndices() { halfedgeIndicesQ.unrequire(); }

// Corner indices
void BaseGeometryInterface::computeCornerIndices() { cornerIndices = mesh.getCornerIndices(); }
void BaseGeometryInterface::requireCornerIndices() { cornerIndicesQ.require(); }
void BaseGeometryInterface::unrequireCornerIndices() { cornerIndicesQ.unrequire(); }


// Face indices
void BaseGeometryInterface::computeFaceIndices() { faceIndices = mesh.getFaceIndices(); }
void BaseGeometryInterface::requireFaceIndices() { faceIndicesQ.require(); }
void BaseGeometryInterface::unrequireFaceIndices() { faceIndicesQ.unrequire(); }


// Boundary loop indices
void BaseGeometryInterface::computeBoundaryLoopIndices() { boundaryLoopIndices = mesh.getBoundaryLoopIndices(); }
void BaseGeometryInterface::requireBoundaryLoopIndices() { boundaryLoopIndicesQ.require(); }
void BaseGeometryInterface::unrequireBoundaryLoopIndices() { boundaryLoopIndicesQ.unrequire(); }


} // namespace surface
} // namespace geometrycentral
