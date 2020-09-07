#pragma once

#include "geometrycentral/utilities/utilities.h"
#include "geometrycentral/utilities/vector2.h"
#include "geometrycentral/utilities/vector3.h"

#include <cstdlib>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace geometrycentral {
namespace surface {

class SimplePolygonMesh {
public:
  SimplePolygonMesh();
  SimplePolygonMesh(std::string meshFilename, std::string type = "");
  SimplePolygonMesh(std::istream& in, std::string type);
  SimplePolygonMesh(const std::vector<std::vector<size_t>>& polygons_, const std::vector<Vector3>& vertexCoordinates_);
  SimplePolygonMesh(const std::vector<std::vector<size_t>>& polygons_, const std::vector<Vector3>& vertexCoordinates_,
                    const std::vector<std::vector<Vector2>>& paramCoordinates_);

  // == Mesh data
  std::vector<std::vector<size_t>> polygons;
  std::vector<Vector3> vertexCoordinates;
  std::vector<std::vector<Vector2>> paramCoordinates; // optional UV coords, in correspondence with polygons array

  // == Accessors
  inline size_t nFaces() const { return polygons.size(); }
  inline size_t nVertices() const { return vertexCoordinates.size(); }
  inline bool hasParameterization() const { return !paramCoordinates.empty(); }

  // == Mutators

  // Mutate this mesh by merging vertices with identical floating point positions.
  // Useful for loading .stl files, which don't contain information about which
  // triangle corners meet at vertices.
  void mergeIdenticalVertices();

  // Mutate this mesh by removing any entries in vertexCoordinates which appear in any polygon. Update polygon indexing
  // accordingly.
  // Returns a index translation vector mapping old indices to new, such that vec[ind_old] == ind_new. Holds INVALID_IND
  // for removed vertices.
  std::vector<size_t> stripUnusedVertices();

  // Mutate this mesh by removing any faces with repeated vertices.
  void stripFacesWithDuplicateVertices();

  // Mutate this mesh and by naively triangulating polygons
  void triangulate();

  // Empty all data arrays
  void clear();


  // === Input & ouput

  void readMeshFromFile(std::istream& in, std::string type);
  void readMeshFromFile(std::string filename, std::string type = "");
  void readMeshFromFile(std::string filename, std::string type,
                        std::string& detectedType); // also returns type string for filetype that was used
  void writeMesh(std::ostream& out, std::string type);
  void writeMesh(std::string filename, std::string type = "");


private:
  std::string detectFileType(std::string filename);

  // Read helpers
  void readMeshFromObjFile(std::istream& in);
  void readMeshFromPlyFile(std::istream& in);
  void readMeshFromStlFile(std::istream& in);
  void readMeshFromOffFile(std::istream& in);
  void readMeshFromAsciiStlFile(std::istream& in);
  void readMeshFromBinaryStlFile(std::istream& in);

  // Write helpers
  void writeMeshObj(std::ostream& out);
};

std::unique_ptr<SimplePolygonMesh> unionMeshes(const std::vector<SimplePolygonMesh>& meshes);

} // namespace surface
} // namespace geometrycentral
