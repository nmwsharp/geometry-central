#pragma once

#include <cstdlib>
#include <vector>
#include <memory>
#include <fstream>
#include <stdexcept>

#include <vector3.h>
#include <utilities.h>
#include <geometry.h>

namespace geometrycentral {

class PolygonSoupMesh {
 public:
  PolygonSoupMesh();
  PolygonSoupMesh(std::string meshFilename);
  PolygonSoupMesh(const std::vector<std::vector<size_t>>& polygons_,
                  const std::vector<Vector3>& vertexCoordinates_);

  // Mutate this mesh and by naively triangulating polygons
  void triangulate();

  // Mesh data
  std::vector<std::vector<size_t>> polygons;
  std::vector<Vector3> vertexCoordinates;

 private:
  void readMeshFromFile(std::string filename);
};

} // namespace geometrycentral