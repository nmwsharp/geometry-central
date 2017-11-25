#include "geometrycentral/halfedge_mesh.h"

namespace geometrycentral {

HalfedgeMeshDataTransfer::HalfedgeMeshDataTransfer(HalfedgeMesh* oldMesh_,
                                                   HalfedgeMesh* newMesh_)
    : oldMesh(oldMesh_), newMesh(newMesh_) {
  heMap = HalfedgeData<HalfedgePtr>(oldMesh);
  vMap = VertexData<VertexPtr>(oldMesh);
  eMap = EdgeData<EdgePtr>(oldMesh);
  fMap = FaceData<FacePtr>(oldMesh);
}

HalfedgeMeshDataTransfer::HalfedgeMeshDataTransfer()
    : oldMesh(nullptr), newMesh(nullptr) {}

}  // namespace geometrycentral
