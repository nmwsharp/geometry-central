#pragma once

#include <cstdlib>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include <geometrycentral/utilities/utilities.h>
#include <geometrycentral/utilities/vector2.h>
#include <geometrycentral/utilities/vector3.h>

namespace geometrycentral {

class PolygonSoupMesh {
public:
  PolygonSoupMesh();
  PolygonSoupMesh(std::string meshFilename, std::string type = "");
  PolygonSoupMesh(const std::vector<std::vector<size_t>>& polygons_, const std::vector<Vector3>& vertexCoordinates_);

  // Mutate this mesh and by naively triangulating polygons
  void triangulate();

  // Mesh data
  std::vector<std::vector<size_t>> polygons;
  std::vector<Vector3> vertexCoordinates;

  // Mutate this mesh by merging vertices with identical floating point positions.
  // Useful for loading .stl files, which don't contain information about which
  // triangle corners meet at vertices.
  void mergeIdenticalVertices();

private:
  void readMeshFromObjFile(std::string filename);
  void readMeshFromStlFile(std::string filename);
  void readMeshFromAsciiStlFile(std::ifstream& in);
  void readMeshFromBinaryStlFile(std::ifstream in);
};

} // namespace geometrycentral
