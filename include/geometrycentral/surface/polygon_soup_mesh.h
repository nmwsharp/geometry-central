#pragma once

#include <cstdlib>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <vector>

#include "geometrycentral/utilities/vector3.h"
#include <geometrycentral/utilities/utilities.h>

namespace geometrycentral {
namespace surface {

class PolygonSoupMesh {
public:
  PolygonSoupMesh();
  PolygonSoupMesh(std::string meshFilename);
  PolygonSoupMesh(const std::vector<std::vector<size_t>>& polygons_, const std::vector<Vector3>& vertexCoordinates_);

  // Mutate this mesh and by naively triangulating polygons
  void triangulate();

  // Mesh data
  std::vector<std::vector<size_t>> polygons;
  std::vector<Vector3> vertexCoordinates;

private:
  void readMeshFromFile(std::string filename);
};

std::unique_ptr<PolygonSoupMesh> unionMeshes(const std::vector<PolygonSoupMesh>& soups);

} // namespace surface
} // namespace geometrycentral
