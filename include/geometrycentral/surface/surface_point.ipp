#pragma once

namespace geometrycentral {
namespace surface {

template <typename T>
inline T SurfacePoint::interpolate(const VertexData<T>& data) const {

  switch (type) {
  case SurfacePointType::Vertex: {
    return data[vertex];
    break;
  }
  case SurfacePointType::Edge: {
    T valTail = data[edge.halfedge().vertex()];
    T valTip = data[edge.halfedge().twin().vertex()];
    return (1. - tEdge) * valTail + tEdge * valTip;
    break;
  }
  case SurfacePointType::Face: {
    T valA = data[face.halfedge().vertex()];
    T valB = data[face.halfedge().next().vertex()];
    T valC = data[face.halfedge().next().next().vertex()];

    return (faceCoords.x * valA) + (faceCoords.y * valB) + (faceCoords.z * valC);
    break;
  }
  }

  throw std::logic_error("bad switch");
  return data[vertex];
}


} // namespace surface
} // namespace geometrycentral

