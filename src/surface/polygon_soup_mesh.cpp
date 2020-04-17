#include <map>
#include <set>
#include <sstream>
#include <string>

#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/polygon_soup_mesh.h"

#include "happly.h"

namespace geometrycentral {
namespace surface {

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
  } else if (type == "ply") {
    readMeshFromPlyFile(meshFilename);
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
// TODO namespace / move to utility
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
void PolygonSoupMesh::readMeshFromObjFile(std::string filename) {
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
      Vector3 position;
      ss >> position;

      vertexCoordinates.push_back(position);

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

// Read a .ply file containing a polygon mesh
void PolygonSoupMesh::readMeshFromPlyFile(std::string filename) {
  polygons.clear();
  vertexCoordinates.clear();

  happly::PLYData plyIn(filename);


  std::vector<std::array<double, 3>> vPos = plyIn.getVertexPositions();
  vertexCoordinates.resize(vPos.size());
  for (size_t iV = 0; iV < vPos.size(); iV++) {
    for (int j = 0; j < 3; j++) {
      vertexCoordinates[iV][j] = vPos[iV][j];
    }
  }

  polygons = plyIn.getFaceIndices<size_t>();
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


void PolygonSoupMesh::stripUnusedVertices() {

  // Check which indices are used
  size_t nV = vertexCoordinates.size();
  std::vector<char> vertexUsed(nV, false);
  for (auto poly : polygons) {
    for (auto i : poly) {
      GC_SAFETY_ASSERT(i < nV,
                       "polygon list has index " + std::to_string(i) + " >= num vertices " + std::to_string(nV));
      vertexUsed[i] = true;
    }
  }


  // Re-index
  std::vector<size_t> newInd(nV, INVALID_IND);
  std::vector<Vector3> newVertexCoordinates(nV);
  size_t nNewV = 0;
  for (size_t iOldV = 0; iOldV < nV; iOldV++) {
    if (!vertexUsed[iOldV]) continue;
    size_t iNewV = nNewV++;
    newInd[iOldV] = iNewV;
    newVertexCoordinates[iNewV] = vertexCoordinates[iOldV];
  }
  vertexCoordinates = newVertexCoordinates;

  // Translate the polygon listing
  for (auto& poly : polygons) {
    for (auto& i : poly) {
      i = newInd[i];
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

    // TODO support corner coords
  }

  polygons = newPolygons;
}

void PolygonSoupMesh::writeMesh(std::string filename, std::string type) {

  // Attempt to detect filename
  bool typeGiven = type != "";
  std::string::size_type sepInd = filename.rfind('.');
  if (!typeGiven) {
    if (sepInd != std::string::npos) {
      std::string extension;
      extension = filename.substr(sepInd + 1);

      // Convert to all lowercase
      std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
      type = extension;
    }
  }

  if (type == "obj") {
    return writeMeshObj(filename);
  } else {
    if (typeGiven) {
      throw std::runtime_error("Can not write mesh file type " + type);
    } else {
      throw std::runtime_error("Could not detect file type to write mesh to " + filename);
    }
  }
}

void PolygonSoupMesh::writeMeshObj(std::string filename) {

  std::ofstream out(filename);
  if (!out) {
    throw std::runtime_error("failed to create output stream " + filename);
  }

  // Write header
  out << "# Mesh exported from GeometryCentral" << std::endl;
  out << "#  vertices: " << vertexCoordinates.size() << std::endl;
  out << "#     faces: " << polygons.size() << std::endl;
  out << "#     texture coordinates: NO" << std::endl;
  out << std::endl;

  // Write vertices
  for (Vector3 p : vertexCoordinates) {
    out << "v " << p.x << " " << p.y << " " << p.z << std::endl;
  }

  // Write texture coords (if present)
  for (std::vector<Vector2>& coords : cornerCoords) {
    for (Vector2 c : coords) {
      out << "vt " << c.x << " " << c.y << std::endl;
    }
  }

  // Write faces
  size_t iC = 0;
  for (size_t iF = 0; iF < polygons.size(); iF++) {
    std::vector<size_t>& face = polygons[iF];

    out << "f";
    for (size_t i = 0; i < face.size(); i++) {
      size_t ind = face[i];
      out << " " << (ind + 1);

      if (!cornerCoords.empty()) {
        out << "/" << (iC + 1);
        iC++;
      }
    }
    out << std::endl;
  }
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
