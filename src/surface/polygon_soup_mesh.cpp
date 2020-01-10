#include <map>
#include <set>
#include <sstream>
#include <string>

#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/polygon_soup_mesh.h"

namespace geometrycentral {
namespace surface {

PolygonSoupMesh::PolygonSoupMesh() {}

PolygonSoupMesh::PolygonSoupMesh(std::string meshFilename) { readMeshFromFile(meshFilename); }

PolygonSoupMesh::PolygonSoupMesh(const std::vector<std::vector<size_t>>& polygons_,
                                 const std::vector<Vector3>& vertexCoordinates_)
    : polygons(polygons_), vertexCoordinates(vertexCoordinates_) {}

// String manipulation helpers to parse .obj files
// See http://stackoverflow.com/a/236803
std::vector<std::string>& split(const std::string& s, char delim, std::vector<std::string>& elems) {
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
  return elems;
}
std::vector<std::string> split(const std::string& s, char delim) {
  std::vector<std::string> elems;
  split(s, delim, elems);
  return elems;
}

class Index {
public:
  Index() {}

  Index(long long int v, long long int vt, long long int vn) : position(v), uv(vt), normal(vn) {}

  bool operator<(const Index& i) const {
    if (position < i.position) return true;
    if (position > i.position) return false;
    if (uv < i.uv) return true;
    if (uv > i.uv) return false;
    if (normal < i.normal) return true;
    if (normal > i.normal) return false;

    return false;
  }

  long long int position = -1;
  long long int uv = -1;
  long long int normal = -1;
};

Index parseFaceIndex(const std::string& token) {
  std::stringstream in(token);
  std::string indexString;
  int indices[3] = {1, 1, 1};

  int i = 0;
  while (std::getline(in, indexString, '/')) {
    if (indexString != "\\") {
      std::stringstream ss(indexString);
      ss >> indices[i++];
    }
  }

  // decrement since indices in OBJ files are 1-based
  return Index(indices[0] - 1, indices[1] - 1, indices[2] - 1);
}

// Read a .obj file containing a polygon mesh
void PolygonSoupMesh::readMeshFromFile(std::string filename) {
  // std::cout << "Reading mesh from file: " << filename << std::endl;

  polygons.clear();
  vertexCoordinates.clear();
  cornerCoords.clear();

  // corner UV coords, unpacked below
  std::vector<Vector2> coords;
  std::vector<std::vector<size_t>> polygonCoordInds;

  // Open the file
  std::ifstream in(filename);
  if (!in) throw std::invalid_argument("Could not open mesh file " + filename);

  // parse obj format
  std::string line;
  while (getline(in, line)) {
    std::stringstream ss(line);
    std::string token;

    ss >> token;

    if (token == "v") {
      double x, y, z;
      ss >> x >> y >> z;

      vertexCoordinates.push_back(Vector3{x, y, z});

    } else if (token == "vt") {
      double u, v;
      ss >> u >> v;

      coords.push_back(Vector2{u, v});


    } else if (token == "vn") {
      // Do nothing

    } else if (token == "f") {
      std::vector<size_t> face;
      std::vector<size_t> faceCoordInds;
      while (ss >> token) {
        Index index = parseFaceIndex(token);
        if (index.position < 0) {
          getline(in, line);
          size_t i = line.find_first_not_of("\t\n\v\f\r ");
          index = parseFaceIndex(line.substr(i));
        }

        face.push_back(index.position);

        if (index.uv != -1) {
          faceCoordInds.push_back(index.uv);
        }
      }

      polygons.push_back(face);
      if (!faceCoordInds.empty()) {
        polygonCoordInds.push_back(faceCoordInds);
      }
    }
  }

  // If we got uv coords, unpack them in to per-corner values
  if (!polygonCoordInds.empty()) {
    for (std::vector<size_t>& faceCoordInd : polygonCoordInds) {
      cornerCoords.emplace_back();
      std::vector<Vector2>& faceCoord = cornerCoords.back();
      for (size_t i : faceCoordInd) {
        if (i < coords.size()) faceCoord.push_back(coords[i]);
      }
    }
  }
}

void PolygonSoupMesh::triangulate() {
  std::vector<std::vector<size_t>> newPolygons;

  for (auto poly : polygons) {
    if (poly.size() <= 2) {
      throw std::runtime_error("ERROR: PolygonSoupMesh has degree < 3 polygon");
    }

    for (size_t i = 2; i < poly.size(); i++) {
      std::vector<size_t> tri = {poly[0], poly[i - 1], poly[i]};
      newPolygons.push_back(tri);
    }
  }

  polygons = newPolygons;
}

std::unique_ptr<PolygonSoupMesh> unionMeshes(const std::vector<PolygonSoupMesh>& soups) {

  std::vector<std::vector<size_t>> unionFaces;
  std::vector<Vector3> unionVerts;

  for (const PolygonSoupMesh& soup : soups) {

    size_t offset = unionVerts.size();
    for (Vector3 v : soup.vertexCoordinates) {
      unionVerts.push_back(v);
    }

    for (std::vector<size_t> f : soup.polygons) {
      for (size_t& i : f) i += offset;
      unionFaces.push_back(f);
    }
  }

  return std::unique_ptr<PolygonSoupMesh>(new PolygonSoupMesh(unionFaces, unionVerts));
}

} // namespace surface
} // namespace geometrycentral
