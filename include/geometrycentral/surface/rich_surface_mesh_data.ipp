#pragma once

namespace geometrycentral {
namespace surface {

// The names to use for mesh element types
// clang-format off
template <typename E> std::string plyElementName() { return "X"; }
template<> inline std::string plyElementName<Vertex       >()            { return "vertex";    }
template<> inline std::string plyElementName<Halfedge     >()            { return "halfedge";   }
template<> inline std::string plyElementName<Corner       >()            { return "corner";    }
template<> inline std::string plyElementName<Edge         >()            { return "edge";    }
template<> inline std::string plyElementName<Face         >()            { return "face";    }
template<> inline std::string plyElementName<BoundaryLoop >()            { return "boundaryloop";   }
// clang-format on


// Generic implementations which handle all element types

template <typename E, typename T>
MeshData<E, T> RichSurfaceMeshData::getElementProperty(std::string propertyName) {

  std::string eName = plyElementName<E>();
  std::vector<T> rawData = plyData.getElement(eName).getProperty<T>(propertyName);

  if (rawData.size() != nElements<E>(mesh)) {
    throw std::runtime_error("Property " + propertyName + " does not have size equal to number of " + eName);
  }

  MeshData<E, T> result(*mesh);
  size_t i = 0;
  for (E e : iterateElements<E>(mesh)) {
    result[e] = rawData[i];
    i++;
  }

  return result;
}

template <typename E, typename T>
void RichSurfaceMeshData::addElementProperty(std::string propertyName, const MeshData<E, T>& data) {

  std::string eName = plyElementName<E>();

  // Make sure the element exists
  if (!plyData.hasElement(eName)) {
    plyData.addElement(eName, nElements<E>(mesh));
  }

  std::vector<T> vec;
  vec.reserve(nElements<E>(mesh));
  for (E e : iterateElements<E>(mesh)) {
    vec.push_back(data[e]);
  }

  plyData.getElement(eName).addProperty<T>(propertyName, vec);
}


// == Nicely-named aliases

// = getters

template <class T>
VertexData<T> RichSurfaceMeshData::getVertexProperty(std::string propertyName) {
  return getElementProperty<Vertex, T>(propertyName);
}

template <class T>
HalfedgeData<T> RichSurfaceMeshData::getHalfedgeProperty(std::string propertyName) {
  return getElementProperty<Halfedge, T>(propertyName);
}

template <class T>
CornerData<T> RichSurfaceMeshData::getCornerProperty(std::string propertyName) {
  return getElementProperty<Corner, T>(propertyName);
}

template <class T>
EdgeData<T> RichSurfaceMeshData::getEdgeProperty(std::string propertyName) {
  return getElementProperty<Edge, T>(propertyName);
}

template <class T>
FaceData<T> RichSurfaceMeshData::getFaceProperty(std::string propertyName) {
  return getElementProperty<Face, T>(propertyName);
}

template <class T>
BoundaryLoopData<T> RichSurfaceMeshData::getBoundaryLoopProperty(std::string propertyName) {
  return getElementProperty<BoundaryLoop, T>(propertyName);
}

// = setters

template <class T>
void RichSurfaceMeshData::addVertexProperty(std::string propertyName, const VertexData<T>& data) {
  return addElementProperty<Vertex, T>(propertyName, data);
}

template <class T>
void RichSurfaceMeshData::addHalfedgeProperty(std::string propertyName, const HalfedgeData<T>& data) {
  return addElementProperty<Halfedge, T>(propertyName, data);
}

template <class T>
void RichSurfaceMeshData::addCornerProperty(std::string propertyName, const CornerData<T>& data) {
  return addElementProperty<Corner, T>(propertyName, data);
}

template <class T>
void RichSurfaceMeshData::addEdgeProperty(std::string propertyName, const EdgeData<T>& data) {
  return addElementProperty<Edge, T>(propertyName, data);
}

template <class T>
void RichSurfaceMeshData::addFaceProperty(std::string propertyName, const FaceData<T>& data) {
  return addElementProperty<Face, T>(propertyName, data);
}

template <class T>
void RichSurfaceMeshData::addBoundaryLoopProperty(std::string propertyName, const BoundaryLoopData<T>& data) {
  return addElementProperty<BoundaryLoop, T>(propertyName, data);
}


} // namespace surface
} // namespace geometrycentral
