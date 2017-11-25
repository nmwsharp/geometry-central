#include "geometrycentral/halfedge_mesh.h"

namespace geometrycentral {

template<class T> VertexData<T> HalfedgeMeshDataTransfer::transfer(VertexData<T>& inData) {

    VertexData<T> outData(newMesh);
    for(VertexPtr v : oldMesh->vertices()) {
        outData[vMap[v]] = inData[v];
    }

    return outData;
}

template<class T> FaceData<T> HalfedgeMeshDataTransfer::transfer(FaceData<T>& inData) {

    FaceData<T> outData(newMesh);
    for(FacePtr f : oldMesh->faces()) {
        outData[fMap[f]] = inData[f];
    }

    return outData;
}

template<class T> EdgeData<T> HalfedgeMeshDataTransfer::transfer(EdgeData<T>& inData) {

    EdgeData<T> outData(newMesh);
    for(EdgePtr e : oldMesh->edges()) {
        outData[eMap[e]] = inData[e];
    }

    return outData;
}

template<class T> HalfedgeData<T> HalfedgeMeshDataTransfer::transfer(HalfedgeData<T>& inData) {

    HalfedgeData<T> outData(newMesh);
    for(HalfedgePtr he : oldMesh->halfedges()) {
        outData[heMap[he]] = inData[he];
    }

    return outData;
}

template<class T> CornerData<T> HalfedgeMeshDataTransfer::transfer(CornerData<T>& inData) {

    CornerData<T> outData(newMesh);
    for(HalfedgePtr he : oldMesh->halfedges()) {
        outData[heMap[he].corner()] = inData[he.corner()];
    }

    return outData;

}

} // namespace geometry central