#pragma once

#include "geometrycentral/surface/common_subdivision.h"

namespace geometrycentral {
namespace surface {

template <typename T>
VertexData<T> CommonSubdivision::interpolateAcrossA(const VertexData<T>& dataA) const {
  checkMeshConstructed();
  VertexData<T> interp(*mesh);
  for (Vertex v : mesh->vertices()) {
    interp[v] = sourcePoints[v]->posA.interpolate(dataA);
  }
  return interp;
}

template <typename T>
VertexData<T> CommonSubdivision::interpolateAcrossB(const VertexData<T>& dataB) const {
  checkMeshConstructed();
  VertexData<T> interp(*mesh);
  for (Vertex v : mesh->vertices()) {
    interp[v] = sourcePoints[v]->posB.interpolate(dataB);
  }
  return interp;
}

template <typename T>
FaceData<T> CommonSubdivision::copyFromA(const FaceData<T>& dataA) const {
  checkMeshConstructed();
  FaceData<T> out(*mesh);
  for (Face f : mesh->faces()) {
    out[f] = dataA[sourceFaceA[f]];
  }
  return out;
}

template <typename T>
FaceData<T> CommonSubdivision::copyFromB(const FaceData<T>& dataB) const {
  checkMeshConstructed();
  FaceData<T> out(*mesh);
  for (Face f : mesh->faces()) {
    out[f] = dataB[sourceFaceB[f]];
  }
  return out;
}


} // namespace surface
} // namespace geometrycentral
