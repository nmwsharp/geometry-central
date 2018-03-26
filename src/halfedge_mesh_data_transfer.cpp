#include "geometrycentral/halfedge_mesh.h"

namespace geometrycentral {

HalfedgeMeshDataTransfer::HalfedgeMeshDataTransfer(HalfedgeMesh* oldMesh_, HalfedgeMesh* newMesh_)
    : oldMesh(oldMesh_), newMesh(newMesh_) {
  heMap = HalfedgeData<HalfedgePtr>(oldMesh);
  vMap = VertexData<VertexPtr>(oldMesh);
  eMap = EdgeData<EdgePtr>(oldMesh);
  fMap = FaceData<FacePtr>(oldMesh);
}

HalfedgeMeshDataTransfer::HalfedgeMeshDataTransfer() : oldMesh(nullptr), newMesh(nullptr) {}


void HalfedgeMeshDataTransfer::generateReverseMaps() {

  // Halfedges
  heMapBack = HalfedgeData<HalfedgePtr>(newMesh);
  for(HalfedgePtr he : oldMesh->allHalfedges()) {
    heMapBack[heMap[he]] = he;
  }
  
  // Vertices 
  vMapBack = VertexData<VertexPtr>(newMesh);
  for(VertexPtr v : oldMesh->vertices()) {
    vMapBack[vMap[v]] = v;
  }

  // Edges
  eMapBack = EdgeData<EdgePtr>(newMesh);
  for(EdgePtr e : oldMesh->edges()) {
    eMapBack[eMap[e]] = e;
  }
  
  // Faces 
  fMapBack = FaceData<FacePtr>(newMesh);
  for(FacePtr f : oldMesh->faces()) {
    fMapBack[fMap[f]] = f;
  }
  for(FacePtr f : oldMesh->boundaryLoops()) {
    fMapBack[fMap[f]] = f;
  }

}

} // namespace geometrycentral
