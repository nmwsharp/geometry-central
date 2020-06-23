#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/edge_length_geometry.h"

#include "happly.h"

#include <fstream>
#include <iostream>
#include <string>


namespace geometrycentral {
namespace surface {

// A reader/writer class supporting direct interop between the GeometryCentral mesh types and the .ply file
// format. Allows storing and retriving of VertexData<> etc containers.
// No operations are valid if the mesh is modified after the creation of the reader/writer.
class RichSurfaceMeshData {

  // TODO give helpers to store Vector2 and Vector3 types easily, rather than having to store each component manually.

public:
  // Construct by reading from file, mapping the elements on to an existing mesh.
  // To simultaneously read the mesh encoded by the file, see the static method readMeshAndData below.
  RichSurfaceMeshData(SurfaceMesh& mesh_, std::string filename);
  RichSurfaceMeshData(SurfaceMesh& mesh_, std::istream& in);

  // Construct a data object. Connectivity will be added automatically, geometry and other data fields can be added.
  RichSurfaceMeshData(SurfaceMesh& mesh_);

  // Convenience factory method to simultaneously read the mesh from a file and load the ply object
  static std::tuple<std::unique_ptr<SurfaceMesh>, std::unique_ptr<RichSurfaceMeshData>>
  readMeshAndData(std::string filename);
  static std::tuple<std::unique_ptr<SurfaceMesh>, std::unique_ptr<RichSurfaceMeshData>>
  readMeshAndData(std::istream& in);
  static std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<RichSurfaceMeshData>>
  readManifoldMeshAndData(std::string filename);
  static std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<RichSurfaceMeshData>>
  readManifoldMeshAndData(std::istream& in);

  // The underlying reader/writer object
  happly::PLYData plyData;

  // Write this object out to file
  void write(std::string filename);
  void write(std::ostream& out);

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

  // Creates a geometry object from vertex positions in the file
  std::unique_ptr<VertexPositionGeometry> getGeometry();
  
  // Creates a geometry object from edge lengths in the file
  std::unique_ptr<EdgeLengthGeometry> getIntrinsicGeometry();

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

  // Add the connectivity of the wrapped mesh to the file
  void addMeshConnectivity();

  // Registers vertex posititions from a geometry object with the file
  void addGeometry(EmbeddedGeometryInterface& geometry);
  
  // Registers edge lengths from a geometry object with the file
  void addIntrinsicGeometry(IntrinsicGeometryInterface& geometry);


private:
  // Create the mesh while constructing
  RichSurfaceMeshData(std::string filename);
  RichSurfaceMeshData(std::istream& in);
  void loadMeshFromFile();

  // The mesh on which the properties in this file are presumed to exist
  SurfaceMesh* mesh = nullptr;

  // Set the names to use for various element types
  std::string vertexName = "vertex";
  std::string faceName = "face";
  std::string edgeName = "edge";
  std::string halfedgeName = "halfedge";
  std::string boundaryLoopName = "boundaryloop";
};


} // namespace surface
} // namespace geometrycentral


#include "geometrycentral/surface/rich_surface_mesh_data.ipp"
