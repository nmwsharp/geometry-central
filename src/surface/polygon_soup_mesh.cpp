#include <map>
#include <set>
#include <sstream>
#include <string>

#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/polygon_soup_mesh.h"

namespace geometrycentral {

PolygonSoupMesh::PolygonSoupMesh() {}

PolygonSoupMesh::PolygonSoupMesh(std::string meshFilename, std::string type) {

  // Attempt to detect filename
  bool typeGiven = type != "";
  std::string::size_type sepInd = meshFilename.rfind('.');
  if (!typeGiven) {
    if (sepInd != std::string::npos) {
      std::string extension;
      extension = meshFilename.substr(sepInd + 1);

      // Convert to all lowercase
      std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
      type = extension;
    }
  }

  if (type == "obj") {
    readMeshFromObjFile(meshFilename);
  } else if (type == "stl") {
    readMeshFromStlFile(meshFilename);
  } else {
    if (typeGiven) {
      throw std::runtime_error("Did not recognize mesh file type " + type);
    } else {
      throw std::runtime_error("Could not detect file type to load mesh from " + meshFilename + ". (Found type " +
                               type + ", but cannot load this)");
    }
  }
}

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

  Index(int v, int vt, int vn) : position(v), uv(vt), normal(vn) {}

  bool operator<(const Index& i) const {
    if (position < i.position) return true;
    if (position > i.position) return false;
    if (uv < i.uv) return true;
    if (uv > i.uv) return false;
    if (normal < i.normal) return true;
    if (normal > i.normal) return false;

    return false;
  }

  int position;
  int uv;
  int normal;
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
void PolygonSoupMesh::readMeshFromObjFile(std::string filename) {
  // std::cout << "Reading mesh from file: " << filename << std::endl;

  polygons.clear();
  vertexCoordinates.clear();

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
      Vector3 position;
      ss >> position;

      vertexCoordinates.push_back(position);

    } else if (token == "vt") {
      // Do nothing
    } else if (token == "vn") {
      // Do nothing

    } else if (token == "f") {
      std::vector<size_t> face;
      while (ss >> token) {
        Index index = parseFaceIndex(token);
        if (index.position < 0) {
          getline(in, line);
          size_t i = line.find_first_not_of("\t\n\v\f\r ");
          index = parseFaceIndex(line.substr(i));
        }

        face.push_back(index.position);
      }

      polygons.push_back(face);
    }
  }
}

// Assumes that first line has already been consumed
void PolygonSoupMesh::readMeshFromAsciiStlFile(std::ifstream& in) {
  std::string line;
  std::stringstream ss;
  size_t lineNum = 1;

  auto assertToken = [&](const std::string& expected) {
    std::string token;
    ss >> token;
    if (token != expected) {
      std::ostringstream errorMessage;
      errorMessage << "Failed to parse ASCII stl file." << std::endl
                   << "Error on line " << lineNum << ". Expected \"" << expected << "\" but token \"" << token << "\""
                   << std::endl
                   << "Full line: \"" << line << "\"" << std::endl;
      throw std::runtime_error(errorMessage.str());
    }
  };

  auto nextLine = [&]() {
    if (!getline(in, line)) {
      return false;
    }

    ss = std::stringstream(line);
    lineNum++;
    return true;
  };

  auto startsWithToken = [](const std::string& str, const std::string& prefix) {
    std::stringstream ss(str);
    std::string token;
    ss >> token;
    return token == prefix;
  };

  // Parse STL file
  while (nextLine() && !startsWithToken(line, "endsolid")) {
    assertToken("facet");
    assertToken("normal");

    // TODO: store this normal?
    Vector3 normal;
    ss >> normal;

    nextLine();

    assertToken("outer");
    assertToken("loop");

    std::vector<size_t> face;
    while (nextLine() && !startsWithToken(line, "endloop")) {
      assertToken("vertex");

      Vector3 position;
      ss >> position;
      vertexCoordinates.push_back(position);

      face.push_back(vertexCoordinates.size() - 1);
    }

    nextLine();
    assertToken("endfacet");

    // Orient face using normal
    Vector3 faceNormal = cross(vertexCoordinates[face[1]] - vertexCoordinates[face[0]],
                               vertexCoordinates[face[2]] - vertexCoordinates[face[0]]);
    if (dot(faceNormal, normal) < 0) {
      std::reverse(std::begin(face), std::end(face));
    }

    polygons.push_back(face);
  }
}

void PolygonSoupMesh::readMeshFromBinaryStlFile(std::ifstream in) {
  auto parseVector3 = [&](std::ifstream& in) {
    char buffer[3 * sizeof(float)];
    in.read(buffer, 3 * sizeof(float));
    float* fVec = (float*)buffer;
    return Vector3{fVec[0], fVec[1], fVec[2]};
  };

  char header[80];
  char nTriangleChars[4];
  in.read(header, 80);
  in.read(nTriangleChars, 4);
  unsigned int* intPtr = (unsigned int*)nTriangleChars;
  size_t nTriangles = *intPtr;

  for (size_t iT = 0; iT < nTriangles; ++iT) {
    // TODO: store this normal?
    Vector3 normal = parseVector3(in);
    std::vector<size_t> face;
    for (size_t iV = 0; iV < 3; ++iV) {
      vertexCoordinates.push_back(parseVector3(in));
      face.push_back(vertexCoordinates.size() - 1);
    }

    // Orient face using normal
    Vector3 faceNormal = cross(vertexCoordinates[face[1]] - vertexCoordinates[face[0]],
                               vertexCoordinates[face[2]] - vertexCoordinates[face[0]]);
    if (dot(faceNormal, normal) < 0) {
      std::reverse(std::begin(face), std::end(face));
    }

    polygons.push_back(face);
    char dummy[2];
    in.read(dummy, 2);
  }
}

// Read a .stl file containing a polygon mesh
void PolygonSoupMesh::readMeshFromStlFile(std::string filename) {
  // TODO: read stl file name
  polygons.clear();
  vertexCoordinates.clear();

  // Open the file
  std::ifstream in(filename);
  if (!in) throw std::invalid_argument("Could not open mesh file " + filename);

  // parse stl format
  std::string line;
  getline(in, line);
  std::stringstream ss(line);
  std::string token;
  ss >> token;
  if (token == "solid") {
    readMeshFromAsciiStlFile(in);
  } else {
    readMeshFromBinaryStlFile(std::ifstream(filename, std::ios::in | std::ios::binary));
  }
}

// Mutate this mesh by merging vertices with identical floating point positions.
// Useful for loading .stl files, which don't contain information about which
// triangle corners meet at vertices.
void PolygonSoupMesh::mergeIdenticalVertices() {
  std::vector<Vector3> compressedPositions;
  // Store mapping from original vertex index to merged vertex index
  std::vector<size_t> compressVertex;
  compressVertex.reserve(vertexCoordinates.size());

  std::unordered_map<Vector3, size_t> canonicalIndex;

  for (size_t iV = 0; iV < vertexCoordinates.size(); ++iV) {
    Vector3 v = vertexCoordinates[iV];
    auto it = canonicalIndex.find(v);

    // Check if vertex exists in map or not
    if (it == canonicalIndex.end()) {
      compressedPositions.push_back(v);
      size_t vecIndex = compressedPositions.size() - 1;
      canonicalIndex[v] = vecIndex;
      compressVertex.push_back(vecIndex);
    } else {
      size_t vecIndex = it->second;
      compressVertex.push_back(vecIndex);
    }
  }

  vertexCoordinates = std::move(compressedPositions);

  // Update face indices
  for (std::vector<size_t>& face : polygons) {
    for (size_t& iV : face) {
      iV = compressVertex[iV];
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

} // namespace geometrycentral
