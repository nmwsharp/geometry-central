#pragma once


#include "geometrycentral/geometry.h"

#include <fstream>
#include <iostream>
#include <string>

#include "happly.h"

// === Reader and writing classes supporting direct interop between the GeometryCentral mesh types and the .ply file
// === format.

namespace geometrycentral {

class PlyHalfedgeMeshData {

public:
  PlyHalfedgeMeshData(std::string filename, bool verbose = false);
  ~PlyHalfedgeMeshData();

  Geometry<Euclidean>* getMesh();

  // === Get properties as geometrycentral types.
  // Will fail with an error if not possible.
  template <class T>
  VertexData<T> getVertexProperty(std::string propertyName);

  template <class T>
  FaceData<T> getFaceProperty(std::string propertyName);

  // === Convenience getters

  // Looks for vertex colors, either as a uchar or a float
  VertexData<Vector3> getVertexColors();

  // === Set properties as geometrycentral types.
  template <class T>
  void addVertexProperty(std::string propertyName, VertexData<T>& vData);

  template <class T>
  void addFaceProperty(std::string propertyName, FaceData<T>& fData);


  // TODO implement creating one of these from a HalfedgeMesh/geometry object

  void write(std::string filename);

private:
  // File data
  std::string filename;
  happly::PLYData* plyData = nullptr;

  bool isBinary;
  float version;


  // Options
  bool verbose;
  std::string vertexName = "vertex";
  std::string faceName = "face";
  std::string edgeName = "edge";
  std::string halfedgeName = "halfedge";
  std::string boundaryLoopName = "boundaryloop";

  // The mesh object corresponding to this file. We keep it around to return VertexData<T> (etc) objects defined on the
  // mesh. These are constructed when getMesh() is called for the first time, and assumed to be valid so long as
  // requests for data are made. User is responsible for memory.
  HalfedgeMesh* mesh = nullptr;
  Geometry<Euclidean>* geometry = nullptr;
};

// === Implementations
template <class T>
VertexData<T> PlyHalfedgeMeshData::getVertexProperty(std::string propertyName) {
  if (mesh == nullptr) {
    getMesh();
  }

  std::vector<T> rawData = plyData->getElement(vertexName).getProperty<T>(propertyName);

  if (rawData.size() != mesh->nVertices()) {
    throw std::runtime_error("Property " + propertyName + " does not have size equal to number of vertices");
  }

  VertexData<T> result(mesh);
  for (size_t i = 0; i < mesh->nVertices(); i++) {
    result[mesh->vertex(i)] = rawData[i];
  }

  return result;
}


template <class T>
FaceData<T> PlyHalfedgeMeshData::getFaceProperty(std::string propertyName) {
  if (mesh == nullptr) {
    getMesh();
  }

  std::vector<T> rawData = plyData->getElement(faceName).getProperty<T>(propertyName);

  if (rawData.size() != mesh->nVertices()) {
    throw std::runtime_error("Property " + propertyName + " does not have size equal to number of vertices");
  }

  VertexData<T> result(mesh);
  for (size_t i = 0; i < mesh->nVertices(); i++) {
    result[mesh->vertex(i)] = rawData[i];
  }

  return result;
}

template <class T>
void PlyHalfedgeMeshData::addVertexProperty(std::string propertyName, VertexData<T>& vData) {

  std::vector<T> vec;
  for (VertexPtr v : mesh->vertices()) {
    vec.push_back(vData[v]);
  }

  plyData->getElement(vertexName).addProperty(new happly::TypedProperty<T>(propertyName, vec));
}


template <class T>
void PlyHalfedgeMeshData::addFaceProperty(std::string propertyName, FaceData<T>& fData) {

  std::vector<T> vec;
  for (FacePtr f : mesh->faces()) {
    vec.push_back(fData[f]);
  }

  plyData->getElement(faceName).addProperty(new happly::TypedProperty<T>(propertyName, vec));
}
}
