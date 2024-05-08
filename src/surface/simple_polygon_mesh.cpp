#include "geometrycentral/surface/simple_polygon_mesh.h"

#include "happly.h"
#include "nanort/nanort.h"

#include <algorithm>
#include <deque>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <unordered_set>
// For strncmp
#include <string.h>


namespace geometrycentral {
namespace surface {

SimplePolygonMesh::SimplePolygonMesh() {}

SimplePolygonMesh::SimplePolygonMesh(std::string meshFilename, std::string type) {
  readMeshFromFile(meshFilename, type);
}

SimplePolygonMesh::SimplePolygonMesh(std::istream& in, std::string type) { readMeshFromFile(in, type); }

SimplePolygonMesh::SimplePolygonMesh(const std::vector<std::vector<size_t>>& polygons_,
                                     const std::vector<Vector3>& vertexCoordinates_)
    : polygons(polygons_), vertexCoordinates(vertexCoordinates_) {}

SimplePolygonMesh::SimplePolygonMesh(const std::vector<std::vector<size_t>>& polygons_,
                                     const std::vector<Vector3>& vertexCoordinates_,
                                     const std::vector<std::vector<Vector2>>& paramCoordinates_)
    : polygons(polygons_), vertexCoordinates(vertexCoordinates_), paramCoordinates(paramCoordinates_) {}


namespace { // helpers for parsing

// String manipulation helpers to parse .obj files
// See http://stackoverflow.com/a/236803
// ONEDAY: move to utility?
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

std::vector<std::string> supportedMeshTypes = {"obj", "ply", "stl", "off"};

} // namespace


std::string SimplePolygonMesh::detectFileType(std::string filename) {
  std::string::size_type sepInd = filename.rfind('.');
  std::string type;

  if (sepInd != std::string::npos) {
    std::string extension;
    extension = filename.substr(sepInd + 1);

    // Convert to all lowercase
    std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
    type = extension;
  } else {
    throw std::runtime_error("Could not auto-detect file type to load mesh from " + filename);
  }

  // Check if this is one of the filetypes we're aware of
  if (std::find(std::begin(supportedMeshTypes), std::end(supportedMeshTypes), type) == std::end(supportedMeshTypes)) {
    throw std::runtime_error("Detected file type " + type + " to load mesh from " + filename +
                             ". This is not a supported file type.");
  }

  return type;
}

void SimplePolygonMesh::readMeshFromFile(std::string filename, std::string type) {
  std::string unused;
  readMeshFromFile(filename, type, unused);
}

void SimplePolygonMesh::readMeshFromFile(std::string filename, std::string type, std::string& detectedType) {

  // Attempt to detect filename
  bool typeGiven = type != "";
  if (!typeGiven) {
    type = detectFileType(filename);
  }

  // == Open the file and load it
  // NOTE: Intentionally always open the stream as binary, even though some of the subsequent formats are plaintext and
  // others are binary.  The only real difference is that non-binary mode performs automatic translation of line ending
  // characters (e.g. \r\n --> \n from DOS). However, this behavior is platform-dependent and having platform-dependent
  // behavior seems more confusing then just handling the newlines properly in the parsers.
  std::ifstream inStream(filename, std::ios::binary);
  if (!inStream) throw std::runtime_error("couldn't open file " + filename);
  readMeshFromFile(inStream, type);

  detectedType = type;
}

void SimplePolygonMesh::readMeshFromFile(std::istream& in, std::string type) {

  if (type == "obj") {
    readMeshFromObjFile(in);
  } else if (type == "stl") {
    readMeshFromStlFile(in);
  } else if (type == "ply") {
    readMeshFromPlyFile(in);
  } else if (type == "off") {
    readMeshFromOffFile(in);
  } else {
    throw std::runtime_error("Did not recognize mesh file type " + type);
  }
}

// Read a .obj file containing a polygon mesh
void SimplePolygonMesh::readMeshFromObjFile(std::istream& in) {
  clear();

  // corner UV coords, unpacked below
  std::vector<Vector2> coords;
  std::vector<std::vector<size_t>> polygonCoordInds;

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
      paramCoordinates.emplace_back();
      std::vector<Vector2>& faceCoord = paramCoordinates.back();
      for (size_t i : faceCoordInd) {
        if (i < coords.size()) faceCoord.push_back(coords[i]);
      }
    }
  }
}

// Assumes that first line has already been consumed
void SimplePolygonMesh::readMeshFromAsciiStlFile(std::istream& in) {
  clear();

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
  // Eat the header line
  nextLine();
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

void SimplePolygonMesh::readMeshFromBinaryStlFile(std::istream& in) {
  clear();

  auto parseVector3 = [&](std::istream& in) {
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
void SimplePolygonMesh::readMeshFromStlFile(std::istream& in) {
  clear();

  // Parse the STL format by looking for the keyword "solid"
  // as the first 5 bytes of the file.  If found, this is
  // our indication we are reading an ASCII format STL file.
  std::array<char, 16> buffer;
  std::fill(begin(buffer), end(buffer), 0);
  in.read(buffer.data(), 5);
  // Conver to lower case for string comparison
  std::transform(begin(buffer), end(buffer), begin(buffer), [](char c) -> char { return std::tolower(c); });
  // In both cases, go ahead and rewind the file
  // to the beginning.  We will handle the STL
  // header in each of the specialized read functions.
  in.seekg(-5, std::ios::cur);
  if (strncmp("solid", buffer.data(), 5) == 0) {
    readMeshFromAsciiStlFile(in);
  } else {
    readMeshFromBinaryStlFile(in);
  }
}

void SimplePolygonMesh::readMeshFromOffFile(std::istream& in) {
  clear();

  // == Parse
  auto getNextLine = [&]() {
    std::string line;
    do {
      if (!std::getline(in, line)) {
        throw std::runtime_error("ran out of lines while parsing off file");
      }
    } while (line.size() == 0 || line[0] == '#');
    return line;
  };

  // header
  std::string headerLine = getNextLine();
  if (headerLine.rfind("OFF", 0) != 0) throw std::runtime_error("does not seem to be valid OFF file");

  // counts
  size_t nVert, nFace;
  std::string countLine = getNextLine();
  std::stringstream countStream(countLine);
  countStream >> nVert >> nFace; // ignore nEdges, if present

  // parse vertices
  vertexCoordinates.resize(nVert);
  for (size_t iV = 0; iV < nVert; iV++) {
    std::string vertLine = getNextLine();
    std::stringstream vertStream(vertLine);
    Vector3 p;
    vertStream >> p.x >> p.y >> p.z; // ignore color etc, if present
    vertexCoordinates[iV] = p;
  }

  // = Get face indices
  polygons.resize(nFace);
  for (size_t iF = 0; iF < nFace; iF++) {
    std::string faceLine = getNextLine();
    std::stringstream faceStream(faceLine);

    size_t degree;
    faceStream >> degree;
    std::vector<size_t>& face = polygons[iF];
    for (size_t i = 0; i < degree; i++) {
      size_t ind;
      faceStream >> ind;
      face.push_back(ind);
    }
  }
}


// Read a .ply file containing a polygon mesh
void SimplePolygonMesh::readMeshFromPlyFile(std::istream& in) {
  clear();

  happly::PLYData plyIn(in);

  std::vector<std::array<double, 3>> vPos = plyIn.getVertexPositions();
  vertexCoordinates.resize(vPos.size());
  for (size_t iV = 0; iV < vPos.size(); iV++) {
    for (int j = 0; j < 3; j++) {
      vertexCoordinates[iV][j] = vPos[iV][j];
    }
  }

  polygons = plyIn.getFaceIndices();
}

// Mutate this mesh by merging vertices with identical floating point positions.
// Useful for loading .stl files, which don't contain information about which
// triangle corners meet at vertices.
void SimplePolygonMesh::mergeIdenticalVertices() {
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


std::vector<size_t> SimplePolygonMesh::stripUnusedVertices() {

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
  std::vector<Vector3> newVertexCoordinates;
  size_t nNewV = 0;
  for (size_t iOldV = 0; iOldV < nV; iOldV++) {
    if (!vertexUsed[iOldV]) continue;
    size_t iNewV = nNewV++;
    newInd[iOldV] = iNewV;
    newVertexCoordinates.push_back(vertexCoordinates[iOldV]);
  }
  vertexCoordinates = newVertexCoordinates;

  // Translate the polygon listing
  for (auto& poly : polygons) {
    for (auto& i : poly) {
      i = newInd[i];
    }
  }

  return newInd;
}

void SimplePolygonMesh::clear() {
  polygons.clear();
  vertexCoordinates.clear();
  paramCoordinates.clear();
}


void SimplePolygonMesh::stripFacesWithDuplicateVertices() { stripInvalidFaceWorker(false, false, true); }

void SimplePolygonMesh::stripInvalidFaces() { stripInvalidFaceWorker(true, true, true); }

void SimplePolygonMesh::stripInvalidFaceWorker(bool removeWithBadInds, bool removeLowDegree,
                                               bool removeWithRepeatedInds) {

  std::vector<std::vector<size_t>> newFaces;
  std::vector<std::vector<Vector2>> newParamCoordinates;

  for (size_t iF = 0; iF < nFaces(); iF++) {

    const std::vector<size_t>& face = polygons[iF];

    bool cullFace = false;

    // Test out-of-bounds vertex indices
    if (removeWithBadInds) {
      for (const size_t ind : face) {
        if (ind >= nVertices()) {
          cullFace = true;
        }
      }
    }

    // Test faces with degree < 3
    if (removeWithBadInds) {
      if (face.size() < 3) {
        cullFace = true;
      }
    }

    // Test out-of-bounds vertex indices
    if (removeWithRepeatedInds) {

      // Generally use a simple search
      size_t D = face.size();
      if (D < 8) {
        for (size_t i = 0; i < D; i++) {
          for (size_t j = i + 1; j < D; j++) {
            if (face[i] == face[j]) cullFace = true;
          }
        }
      }
      // Use a hashset to avoid n^2 for big faces
      else {
        std::unordered_set<size_t> inds;
        for (size_t ind : face) {
          if (inds.find(ind) != inds.end()) cullFace = true;
          inds.insert(ind);
        }
      }
    }


    if (!cullFace) {
      newFaces.push_back(face);
      if (hasParameterization()) {
        newParamCoordinates.push_back(paramCoordinates[iF]);
      }
    }
  }

  // Update the face set for this mesh to be the newly created face set
  polygons = newFaces;
  if (hasParameterization()) {
    paramCoordinates = newParamCoordinates;
  }
}

void SimplePolygonMesh::stripDuplicateFaces() {

  // Build a list of each polygon, with its indices sorted and an index
  // to the original triangle
  struct IndPoly {
    std::vector<size_t> inds;
    size_t orig_iF;
  };
  std::vector<IndPoly> sortpolys(nFaces());
  for (size_t iF = 0; iF < nFaces(); iF++) {
    sortpolys[iF] = IndPoly{polygons[iF], iF};
    std::sort(sortpolys[iF].inds.begin(), sortpolys[iF].inds.end());
  }

  // Sort the overall list lexicographically
  std::sort(sortpolys.begin(), sortpolys.end(), [](const IndPoly& a, const IndPoly& b) { return a.inds < b.inds; });


  // Keep only the first entry from each run of sorted ones
  std::vector<std::vector<size_t>> newPolygons;
  std::vector<std::vector<Vector2>> newParamCoordinates;
  for (size_t iF = 0; iF < nFaces(); /* pass */) {

    // for the first instance of this vertex set, keep it
    size_t keep_iF = sortpolys[iF].orig_iF;
    newPolygons.push_back(polygons[keep_iF]);
    if (hasParameterization()) {
      newParamCoordinates.push_back(paramCoordinates[keep_iF]);
    }

    // for any subsequent entries, mark them for removal
    size_t start_iF = iF;
    iF++;
    while (iF < nFaces() && sortpolys[start_iF].inds == sortpolys[iF].inds) {
      iF++;
    }
  }

  polygons = newPolygons;
  paramCoordinates = newParamCoordinates;
}

void SimplePolygonMesh::triangulate() {
  std::vector<std::vector<size_t>> newPolygons;

  for (auto poly : polygons) {
    if (poly.size() <= 2) {
      throw std::runtime_error("ERROR: SimplePolygonMesh has degree < 3 polygon");
    }

    for (size_t i = 2; i < poly.size(); i++) {
      std::vector<size_t> tri = {poly[0], poly[i - 1], poly[i]};
      newPolygons.push_back(tri);
    }
  }

  polygons = newPolygons;
}

void SimplePolygonMesh::consistentlyOrientFaces(bool outwardOrient) {


  // If we are going to want it, build a raytracing acceleration structure to do the outward orientation step.
  int32_t N_RAYS_PER_FACE = 4;
  std::vector<double> rawPositions;
  std::vector<unsigned int> rawFaces;
  nanort::BVHAccel<double, nanort::TriangleMesh<double>, nanort::TriangleSAHPred<double>,
                   nanort::TriangleIntersector<double>>
      accel;
  double lengthScale = -1;
  auto traceDist = [&](Vector3 root, Vector3 dir) {
    nanort::Ray<double> ray;
    ray.min_t = 1e-6 * lengthScale;
    ray.max_t = 1e1 * lengthScale;
    for (int i = 0; i < 3; i++) ray.org[i] = root[i];
    for (int i = 0; i < 3; i++) ray.dir[i] = dir[i];
    nanort::BVHTraceOptions trace_options;
    nanort::TriangleIntersector<double> triangle_intersector(rawPositions.data(), rawFaces.data(), sizeof(double) * 3);
    bool hit = accel.Traverse(ray, trace_options, triangle_intersector);
    return triangle_intersector.intersection.t;
  };
  std::random_device rd;
  std::mt19937 gen(rd()); // we'll use this to generate rays on triangles
  std::uniform_real_distribution<double> realDist(0., 1.);
  if (outwardOrient) {

    double INF = std::numeric_limits<double>::infinity();
    Vector3 bboxMin{INF, INF, INF};
    Vector3 bboxMax{-INF, -INF, -INF};
    for (Vector3 v : vertexCoordinates) {
      rawPositions.push_back(v.x);
      rawPositions.push_back(v.y);
      rawPositions.push_back(v.z);
      bboxMin = componentwiseMin(bboxMin, v);
      bboxMax = componentwiseMax(bboxMax, v);
    }
    lengthScale = norm(bboxMax - bboxMin);
    for (const std::vector<size_t>& poly : polygons) {
      if (poly.size() != 3) {
        throw std::runtime_error("consistentlyOrientFaces() with outwardOrient=true only supports triangular meshes");
      }
      rawFaces.push_back(poly[0]);
      rawFaces.push_back(poly[1]);
      rawFaces.push_back(poly[2]);
    }


    nanort::TriangleMesh<double> nanortMesh(rawPositions.data(), rawFaces.data(), sizeof(double) * 3);
    nanort::TriangleSAHPred<double> nanortMeshPred(rawPositions.data(), rawFaces.data(), sizeof(double) * 3);
    nanort::BVHBuildOptions<double> options; // Use default options
    bool ret = accel.Build(nFaces(), options, nanortMesh, nanortMeshPred);
    if (!ret) {
      throw std::runtime_error("BVH construction failed");
    }
  }

  // pre-build a connectivity lookup map
  struct FaceSideEntry {
    size_t indA;
    size_t indB;
    size_t iF;
    size_t localInd;
    bool isCanonical;
  };
  struct FaceSideNeighbor {
    bool hasNeighbor;
    size_t neighF;
    bool neighSameOrientation;
  };

  // build out entries for each faceside
  std::vector<FaceSideEntry> entries;
  for (size_t iF = 0; iF < polygons.size(); iF++) {
    const std::vector<size_t>& face = polygons[iF];
    size_t D = face.size();
    for (size_t j = 0; j < D; j++) {
      size_t iV = face[j];
      size_t iVn = face[(j + 1) % D];
      FaceSideEntry entry{std::min(iV, iVn), std::max(iV, iVn), iF, j, iV < iVn};
      entries.push_back(entry);
    }
  }

  // sort so that adjacnet facesides appear next to each other in order
  std::sort(entries.begin(), entries.end(), [](const FaceSideEntry& a, const FaceSideEntry& b) {
    if (a.indA != b.indA) {
      return a.indA < b.indA;
    }
    return a.indB < b.indB;
  });


  // pre-allocate a list of neighbor entries
  std::vector<std::vector<FaceSideNeighbor>> neighbors(nFaces());
  for (size_t iF = 0; iF < polygons.size(); iF++) {
    const std::vector<size_t>& face = polygons[iF];
    size_t D = face.size();
    neighbors[iF] = std::vector<FaceSideNeighbor>(D, FaceSideNeighbor{false, INVALID_IND, true});
  }

  // walk the list matching up pairs of neighboring facesides
  for (size_t i = 0; i < entries.size(); /* pass */) {
    size_t iStart = i;

    const FaceSideEntry& entryStart = entries[i];
    size_t iEnd = i;
    iEnd++;
    // advance the range to cover all matching entries
    while (iEnd < entries.size() && entries[iEnd].indA == entryStart.indA && entries[iEnd].indB == entryStart.indB) {
      iEnd++;
    }

    if (iEnd - iStart == 2) { // exactly two entries --> manifold
      const FaceSideEntry& entryOther = entries[iStart + 1];

      FaceSideNeighbor& nStart = neighbors[entryStart.iF][entryStart.localInd];
      FaceSideNeighbor& nEnd = neighbors[entryOther.iF][entryOther.localInd];

      bool sameOrient = entryStart.isCanonical != entryOther.isCanonical;

      nStart.hasNeighbor = true;
      nStart.neighF = entryOther.iF;
      nStart.neighSameOrientation = sameOrient;

      nEnd.hasNeighbor = true;
      nEnd.neighF = entryStart.iF;
      nEnd.neighSameOrientation = sameOrient;

    } else {
      // nonmanifold --> no neighbor, which is the default for the preallocated entries
    }

    i = iEnd;
  }

  // greedly process manifold patches of faces, flipping faces to be consistently oriented
  std::vector<char> processed(nFaces(), false);
  for (size_t startFace = 0; startFace < polygons.size(); startFace++) {
    if (processed[startFace]) continue;

    // Start processing a new patch
    std::vector<size_t> allFacesInThisPatch;
    std::deque<std::tuple<size_t, bool>> processQueue;
    processQueue.push_back({startFace, false});

    while (!processQueue.empty()) {
      std::tuple<size_t, bool> thisTup = processQueue.front();
      processQueue.pop_front();
      size_t iF = std::get<0>(thisTup);
      bool needsFlip = std::get<1>(thisTup);

      if (processed[iF]) continue;
      processed[iF] = true;

      if (outwardOrient) {
        // we only need this if doing the outward orientation step
        allFacesInThisPatch.push_back(iF);
      }

      if (needsFlip) {
        // actually flip the face
        // (the entries in the neighbor list point at this face are wrong now, but that won't matter)
        std::reverse(polygons[iF].begin(), polygons[iF].end());
        if (hasParameterization()) {
          std::reverse(paramCoordinates[iF].begin(), paramCoordinates[iF].end());
        }
      }

      // add the neighbors for prcessing
      for (const FaceSideNeighbor& n : neighbors[iF]) {
        if (!n.hasNeighbor) {
          continue;
        }
        size_t jF = n.neighF;
        if (processed[jF]) {
          continue;
        }

        bool neighNeedsFlip = (n.neighSameOrientation == needsFlip);
        processQueue.push_back({jF, neighNeedsFlip});
      }
    }

    // we finished processing a maximal patch, now flip it for outwardness
    if (outwardOrient) {

      int32_t outCount = 0;
      for (size_t iF : allFacesInThisPatch) {
        for (int32_t iRay = 0; iRay < N_RAYS_PER_FACE; iRay++) {

          // Generate a point in the face
          double r1 = realDist(gen);
          double r2 = realDist(gen);
          Vector3 bcoord{1. - std::sqrt(r1), std::sqrt(r1) * (1. - r2), std::sqrt(r1) * r2};
          Vector3 p0 = vertexCoordinates[polygons[iF][0]];
          Vector3 p1 = vertexCoordinates[polygons[iF][1]];
          Vector3 p2 = vertexCoordinates[polygons[iF][2]];
          Vector3 root = bcoord[0] * p0 + bcoord[1] * p1 + bcoord[2] * p2;
          Vector3 outwardNormal = unit(cross(p1 - p0, p2 - p0));

          // Trace in both directions
          double outwardHitDist = traceDist(root, outwardNormal);
          double inwardHitDist = traceDist(root, -outwardNormal);

          // Whichever ray goes farther gives a vote for that outward orientation
          if (outwardHitDist > inwardHitDist) {
            outCount++;
          } else {
            outCount--;
          }
        }
      }

      // if there were more votes for inward than outward, flip the whole patch
      if (outCount < 0) {
        for (size_t iF : allFacesInThisPatch) {
          std::reverse(polygons[iF].begin(), polygons[iF].end());
          if (hasParameterization()) {
            std::reverse(paramCoordinates[iF].begin(), paramCoordinates[iF].end());
          }
        }
      }
    }
  }
}

void SimplePolygonMesh::writeMesh(std::string filename, std::string type) {

  // Auto-detect type if needed
  bool typeGiven = type != "";
  if (!typeGiven) {
    type = detectFileType(filename);
  }

  // NOTE if/when we ever start writing binary file formats, will need to open in binary mode for those formats
  std::ofstream outStream(filename);
  if (!outStream) throw std::runtime_error("couldn't open output file " + filename);
  writeMesh(outStream, type);
}

void SimplePolygonMesh::writeMesh(std::ostream& out, std::string type) {
  if (type == "obj") {
    return writeMeshObj(out);
  } else {
    throw std::runtime_error("Write mesh file type " + type + " not supported");
  }
}

void SimplePolygonMesh::writeMeshObj(std::ostream& out) {

  // Make sure we write out at full precision
  out << std::setprecision(std::numeric_limits<double>::max_digits10);

  // Write header
  out << "# Mesh exported from geometry-central" << std::endl;
  out << "#  vertices: " << vertexCoordinates.size() << std::endl;
  out << "#     faces: " << polygons.size() << std::endl;
  out << std::endl;

  // Write vertices
  for (Vector3 p : vertexCoordinates) {
    out << "v " << p.x << " " << p.y << " " << p.z << std::endl;
  }


  // Write texture coords (if present)
  for (std::vector<Vector2>& coords : paramCoordinates) {
    for (Vector2 c : coords) {
      out << "vt " << c.x << " " << c.y << std::endl;
    }
  }

  // Write faces
  size_t iC = 0;
  for (std::vector<size_t>& face : polygons) {
    out << "f";
    for (size_t ind : face) {
      out << " " << (ind + 1);

      if (!paramCoordinates.empty()) {
        out << "/" << (iC + 1);
        iC++;
      }
    }
    out << std::endl;
  }
}

std::unique_ptr<SimplePolygonMesh> unionMeshes(const std::vector<SimplePolygonMesh>& meshes) {

  std::vector<std::vector<size_t>> unionFaces;
  std::vector<std::vector<Vector2>> unionCoords;
  std::vector<Vector3> unionVerts;

  bool keepCoords = true;
  for (const SimplePolygonMesh& mesh : meshes) {
    if (!mesh.hasParameterization()) {
      keepCoords = false;
    }
  }

  for (const SimplePolygonMesh& mesh : meshes) {

    size_t offset = unionVerts.size();
    for (Vector3 v : mesh.vertexCoordinates) {
      unionVerts.push_back(v);
    }

    for (std::vector<size_t> f : mesh.polygons) {
      for (size_t& i : f) i += offset;
      unionFaces.push_back(f);
    }

    if (keepCoords) {
      for (std::vector<Vector2> c : mesh.paramCoordinates) {
        unionCoords.push_back(c);
      }
    }
  }

  return std::unique_ptr<SimplePolygonMesh>(new SimplePolygonMesh(unionFaces, unionVerts, unionCoords));
}

} // namespace surface
} // namespace geometrycentral
