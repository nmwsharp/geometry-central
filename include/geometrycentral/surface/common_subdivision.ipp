#pragma once

#include "geometrycentral/surface/common_subdivision.h"

namespace geometrycentral {
namespace surface {

template <typename T>
VertexData<T> CommonSubdivision::interpolateAcrossA(const VertexData<T>& dataA) const {
  VertexData<T> interp(*mesh);
  for (Vertex v : mesh->vertices()) {
    interp[v] = sourcePoints[v]->posA.interpolate(dataA);
  }
  return interp;
}

template <typename T>
VertexData<T> CommonSubdivision::interpolateAcrossB(const VertexData<T>& dataB) const {
  VertexData<T> interp(*mesh);
  for (Vertex v : mesh->vertices()) {
    interp[v] = sourcePoints[v]->posB.interpolate(dataB);
  }
  return interp;
}


} // namespace surface
} // namespace geometrycentral
