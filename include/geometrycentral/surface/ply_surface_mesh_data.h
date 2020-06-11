#pragma once

#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "happly.h"

#include <fstream>
#include <iostream>
#include <string>


namespace geometrycentral {
namespace surface {

// A reader/writer class supporting direct interop between the GeometryCentral mesh types and the .ply file
// format. Allows storing and retriving of VertexData<> etc containers.
// No operations are valid if the mesh is modified after the creation of the reader/writer.
class PlySurfaceMeshData {

  // TODO implement permutation halfedge mesh storage

  // TODO support intrinsic geometry

  // TODO give helpers to store Vector2 and Vector3 types easily, rather than having to store each component manually.

public:
  // Construct by reading from file, mapping the elements on to an existing mesh.
  // To simultaneously read the mesh encoded by the file, see the static method loadMeshAndData below.
  PlySurfaceMeshData(SurfaceMesh& mesh_, std::string filename);

  // Construct a data object. Connectivity will be added automatically, geometry and other data fields can be added.
  PlySurfaceMeshData(SurfaceMesh& mesh_);

  // Convenience factory method to simultaneously read the mesh from a file and
  static std::tuple<std::unique_ptr<SurfaceMesh>, std::unique_ptr<PlySurfaceMeshData>>
  loadMeshAndData(std::string filename);

  // The mesh on which the properties in this file are presumed to exist
  SurfaceMesh& mesh;

  // The underlying reader/writer object
  happly::PLYData plyData;

  // Write this object out to file
  void write(std::string filename);

  // Options
  happly::DataFormat outputFormat = happly::DataFormat::Binary;

  // === Get properties as geometrycentral types.
  // Will fail with an error if not possible.

  // The generic getter
  // (all the nicely-named versions below are just wrappers around this)
  template <typename E, typename T>
  MeshData<E, T> getElementProperty(std::string propertyName);

  template <class T>
  VertexData<T> getVertexProperty(std::string propertyName);

  template <class T>
  HalfedgeData<T> getHalfedgeProperty(std::string propertyName);

  template <class T>
  CornerData<T> getCornerProperty(std::string propertyName);

  template <class T>
  EdgeData<T> getEdgeProperty(std::string propertyName);

  template <class T>
  FaceData<T> getFaceProperty(std::string propertyName);

  template <class T>
  BoundaryLoopData<T> getBoundaryLoopProperty(std::string propertyName);


  // = Convenience getters

  // Creates vertex posititions from a geometry object with the file
  std::unique_ptr<VertexPositionGeometry> getGeometry();

  // Looks for vertex colors, either as a uchar or a float
  VertexData<Vector3> getVertexColors();


  // === Set properties as geometrycentral types.

  // The generic setter
  // (all the nicely-named versions below are just wrappers around this)
  template <typename E, typename T>
  void addElementProperty(std::string propertyName, const MeshData<E, T>& data);

  template <class T>
  void addVertexProperty(std::string propertyName, const VertexData<T>& data);

  template <class T>
  void addHalfedgeProperty(std::string propertyName, const HalfedgeData<T>& data);

  template <class T>
  void addCornerProperty(std::string propertyName, const CornerData<T>& data);

  template <class T>
  void addEdgeProperty(std::string propertyName, const EdgeData<T>& data);

  template <class T>
  void addFaceProperty(std::string propertyName, const FaceData<T>& data);

  template <class T>
  void addBoundaryLoopProperty(std::string propertyName, const BoundaryLoopData<T>& fData);

  // = Convenience setters

  // Registers vertex posititions from a geometry object with the file
  void addGeometry(EmbeddedGeometryInterface& geometry);


private:
  // Set the names to use for various element types
  std::string vertexName = "vertex";
  std::string faceName = "face";
  std::string edgeName = "edge";
  std::string halfedgeName = "halfedge";
  std::string boundaryLoopName = "boundaryloop";
};


} // namespace surface
} // namespace geometrycentral


#include "geometrycentral/surface/ply_surface_mesh_data.ipp"
