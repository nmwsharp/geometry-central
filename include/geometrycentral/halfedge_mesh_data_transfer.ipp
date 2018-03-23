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
    for(FacePtr f : oldMesh->boundaryLoops()) {
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
    for(HalfedgePtr he : oldMesh->allHalfedges()) {
        outData[heMap[he]] = inData[he];
    }

    return outData;
}

template<class T> CornerData<T> HalfedgeMeshDataTransfer::transfer(CornerData<T>& inData) {

    CornerData<T> outData(newMesh);
    for(HalfedgePtr he : oldMesh->allHalfedges()) {
        outData[heMap[he].corner()] = inData[he.corner()];
    }

    return outData;

}


template<class T> VertexData<T> HalfedgeMeshDataTransfer::transferBack(VertexData<T>& inData) {

    VertexData<T> outData(oldMesh);
    for(VertexPtr v : oldMesh->vertices()) {
        outData[v] = inData[vMap[v]];
    }
    return outData;
}


template<class T> FaceData<T> HalfedgeMeshDataTransfer::transferBack(FaceData<T>& inData) {

    FaceData<T> outData(oldMesh);
    for(FacePtr f : oldMesh->faces()) {
        outData[f] = inData[fMap[f]];
    }
    for(FacePtr f : oldMesh->boundaryLoops()) {
        outData[f] = inData[fMap[f]];
    }
    return outData;
}


template<class T> EdgeData<T> HalfedgeMeshDataTransfer::transferBack(EdgeData<T>& inData) {

    EdgeData<T> outData(oldMesh);
    for(EdgePtr e : oldMesh->edges()) {
        outData[e] = inData[eMap[e]];
    }
    return outData;
}

template<class T> HalfedgeData<T> HalfedgeMeshDataTransfer::transferBack(HalfedgeData<T>& inData) {

    HalfedgeData<T> outData(oldMesh);
    for(HalfedgePtr he : oldMesh->allHalfedges()) {
        outData[he] = inData[heMap[he]];
    }
    return outData;
}

template<class T> CornerData<T> HalfedgeMeshDataTransfer::transferBack(CornerData<T>& inData) {

    CornerData<T> outData(oldMesh);
    for(HalfedgePtr he : oldMesh->allHalfedges()) {
        outData[he.corner()] = inData[heMap[he].corner()];
    }
    return outData;
}


} // namespace geometry central
