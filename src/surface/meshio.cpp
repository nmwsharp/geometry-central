#include "geometrycentral/surface/meshio.h"

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/simple_polygon_mesh.h"
#include "geometrycentral/surface/surface_mesh_factories.h"

#include "happly.h"

#include <iostream>
#include <limits>

using std::cout;
using std::endl;

namespace geometrycentral {
namespace surface {

// ======= Input =======

// Helpers for mesh loading
namespace {
void processLoadedMesh(SimplePolygonMesh& mesh, std::string type) {
  mesh.stripUnusedVertices();

  // Apply any special processing for particular filetypes
  if (type == "stl") {
    mesh.mergeIdenticalVertices();
  }
}


std::vector<Vector3> geometryToStdVector(SurfaceMesh& mesh, EmbeddedGeometryInterface& geometry) {
  geometry.requireVertexPositions();
  std::vector<Vector3> verts(mesh.nVertices());
  size_t i = 0;
  for (Vertex v : mesh.vertices()) {
    verts[i] = geometry.vertexPositions[v];
    i++;
  }
  geometry.unrequireVertexPositions();
  return verts;
}


std::vector<std::vector<Vector2>> paramToStdVector(SurfaceMesh& mesh, CornerData<Vector2>& param) {
  std::vector<std::vector<Vector2>> uv(mesh.nFaces());
  size_t i = 0;
  for (Face f : mesh.faces()) {
    for (Corner c : f.adjacentCorners()) {
      uv[i].push_back(param[c]);
    }
    i++;
  }
  return uv;
}

} // namespace

std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>> loadMesh(std::string filename,
                                                                                                   std::string type) {
  return readManifoldSurfaceMesh(filename, type);
}


// Load a general surface mesh, which might or might not be manifold
std::tuple<std::unique_ptr<SurfaceMesh>, std::unique_ptr<VertexPositionGeometry>> readSurfaceMesh(std::string filename,
                                                                                                  std::string type) {
  std::string loadType;
  SimplePolygonMesh simpleMesh;
  simpleMesh.readMeshFromFile(filename, type, loadType);
  processLoadedMesh(simpleMesh, loadType);
  return makeSurfaceMeshAndGeometry(simpleMesh.polygons, simpleMesh.vertexCoordinates);
}
std::tuple<std::unique_ptr<SurfaceMesh>, std::unique_ptr<VertexPositionGeometry>> readSurfaceMesh(std::istream& in,
                                                                                                  std::string type) {
  std::string loadType = type;
  SimplePolygonMesh simpleMesh;
  simpleMesh.readMeshFromFile(in, type);
  processLoadedMesh(simpleMesh, loadType);
  return makeSurfaceMeshAndGeometry(simpleMesh.polygons, simpleMesh.vertexCoordinates);
}

// Load a manifold surface mesh; an exception will by thrown if the mesh is not manifold.
std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
readManifoldSurfaceMesh(std::string filename, std::string type) {
  std::string loadType;
  SimplePolygonMesh simpleMesh;
  simpleMesh.readMeshFromFile(filename, type, loadType);
  processLoadedMesh(simpleMesh, loadType);
  auto lvals = makeManifoldSurfaceMeshAndGeometry(simpleMesh.polygons, simpleMesh.vertexCoordinates);
  return std::tuple<std::unique_ptr<ManifoldSurfaceMesh>,
                    std::unique_ptr<VertexPositionGeometry>>(std::move(std::get<0>(lvals)),  // mesh
                                                             std::move(std::get<1>(lvals))); // geometry
}
std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
readManifoldSurfaceMesh(std::istream& in, std::string type) {
  std::string loadType = type;
  SimplePolygonMesh simpleMesh;
  simpleMesh.readMeshFromFile(in, type);
  processLoadedMesh(simpleMesh, loadType);

  auto lvals = makeManifoldSurfaceMeshAndGeometry(simpleMesh.polygons, simpleMesh.vertexCoordinates);
  return std::tuple<std::unique_ptr<ManifoldSurfaceMesh>,
                    std::unique_ptr<VertexPositionGeometry>>(std::move(std::get<0>(lvals)),  // mesh
                                                             std::move(std::get<1>(lvals))); // geometry
}

// Load a mesh with UV coordinates, which will be stored as data at triangle
// corners (to allow for UVs that are discontinuous across edges, e.g., at cuts)
std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>,
           std::unique_ptr<CornerData<Vector2>>>
readParameterizedManifoldSurfaceMesh(std::string filename, std::string type) {
  std::string loadType;
  SimplePolygonMesh simpleMesh;
  simpleMesh.readMeshFromFile(filename, type, loadType);
  processLoadedMesh(simpleMesh, loadType);

  return makeManifoldSurfaceMeshAndGeometry(simpleMesh.polygons, {}, simpleMesh.vertexCoordinates,
                                            simpleMesh.paramCoordinates);
}

// Load a mesh with UV coordinates, which will be stored as data at triangle
// corners (to allow for UVs that are discontinuous across edges, e.g., at cuts)
std::tuple<std::unique_ptr<SurfaceMesh>, std::unique_ptr<VertexPositionGeometry>, std::unique_ptr<CornerData<Vector2>>>
readParameterizedSurfaceMesh(std::string filename, std::string type) {
  std::string loadType;
  SimplePolygonMesh simpleMesh;
  simpleMesh.readMeshFromFile(filename, type, loadType);
  processLoadedMesh(simpleMesh, loadType);

  return makeSurfaceMeshAndGeometry(simpleMesh.polygons, {}, simpleMesh.vertexCoordinates, simpleMesh.paramCoordinates);
}

// ======= Output =======


void writeSurfaceMesh(SurfaceMesh& mesh, EmbeddedGeometryInterface& geometry, std::string filename, std::string type) {
  SimplePolygonMesh simpleMesh(mesh.getFaceVertexList(), geometryToStdVector(mesh, geometry));
  simpleMesh.writeMesh(filename, type);
}

void writeSurfaceMesh(SurfaceMesh& mesh, EmbeddedGeometryInterface& geometry, CornerData<Vector2>& texCoords,
                      std::string filename, std::string type) {
  SimplePolygonMesh simpleMesh(mesh.getFaceVertexList(), geometryToStdVector(mesh, geometry),
                               paramToStdVector(mesh, texCoords));
  simpleMesh.writeMesh(filename, type);
}

void writeSurfaceMesh(SurfaceMesh& mesh, EmbeddedGeometryInterface& geometry, std::ostream& out, std::string type) {
  SimplePolygonMesh simpleMesh(mesh.getFaceVertexList(), geometryToStdVector(mesh, geometry));
  simpleMesh.writeMesh(out, type);
}
void writeSurfaceMesh(SurfaceMesh& mesh, EmbeddedGeometryInterface& geometry, CornerData<Vector2>& texCoords,
                      std::ostream& out, std::string type) {
  SimplePolygonMesh simpleMesh(mesh.getFaceVertexList(), geometryToStdVector(mesh, geometry),
                               paramToStdVector(mesh, texCoords));
  simpleMesh.writeMesh(out, type);
}


CornerData<Vector2> packToParam(SurfaceMesh& mesh, VertexData<double>& vals) {
  CornerData<Vector2> out(mesh);
  for (Corner c : mesh.corners()) {
    Vertex v = c.vertex();
    out[c] = Vector2{vals[v], 0.};
  }
  return out;
}
CornerData<Vector2> packToParam(SurfaceMesh& mesh, VertexData<double>& valsX, VertexData<double>& valsY) {
  CornerData<Vector2> out(mesh);
  for (Corner c : mesh.corners()) {
    Vertex v = c.vertex();
    out[c] = Vector2{valsX[v], valsY[v]};
  }
  return out;
}


bool WavefrontOBJ::write(std::string filename, EmbeddedGeometryInterface& geometry) {
  std::ofstream out;
  if (!openStream(out, filename)) return false;

  writeHeader(out, geometry);
  out << "# texture coordinates: NO" << endl;
  cout << endl;

  writeVertices(out, geometry);

  bool useTexCoords = false;
  bool useNormals = false;
  writeFaces(out, geometry, useTexCoords, useNormals);

  return true;
}

bool WavefrontOBJ::write(std::string filename, EmbeddedGeometryInterface& geometry, CornerData<Vector2>& texcoords) {
  std::ofstream out;
  if (!openStream(out, filename)) return false;

  writeHeader(out, geometry);
  out << "# texture coordinates: YES" << endl;
  cout << endl;

  writeVertices(out, geometry);
  writeTexCoords(out, geometry, texcoords);

  bool useTexCoords = true;
  writeFaces(out, geometry, useTexCoords);

  return true;
}

bool WavefrontOBJ::write(std::string filename, EmbeddedGeometryInterface& geometry, CornerData<Vector3>& normals) {
  std::ofstream out;
  if (!openStream(out, filename)) return false;

  writeHeader(out, geometry);
  out << "# texture coordinates: YES" << endl;
  cout << endl;

  writeVertices(out, geometry);
  writeNormals(out, geometry, normals);

  bool useTexCoords = false;
  bool useNormals = true;
  writeFaces(out, geometry, useTexCoords, useNormals);

  return true;
}

bool WavefrontOBJ::write(std::string filename, EmbeddedGeometryInterface& geometry, CornerData<Vector2>& texCoords,
                         CornerData<Vector3>& normals) {
  std::ofstream out;
  if (!openStream(out, filename)) return false;

  writeHeader(out, geometry);
  out << "# texture coordinates: YES" << endl;
  cout << endl;

  writeVertices(out, geometry);
  writeNormals(out, geometry, normals);

  bool useTexCoords = true;
  bool useNormals = true;
  writeFaces(out, geometry, useTexCoords, useNormals);

  return true;
}

bool WavefrontOBJ::openStream(std::ofstream& out, std::string filename) {
  out.open(filename);

  if (!out.is_open()) {
    return false;
  }

  // Use full precision---yes, this makes files bigger, but it also
  // means that saving and then re-loading files doesn't result in
  // unexpected behavior.
  out.precision(std::numeric_limits<double>::max_digits10);

  return true;
}

void WavefrontOBJ::writeHeader(std::ofstream& out, EmbeddedGeometryInterface& geometry) {
  out << "# Mesh exported from GeometryCentral" << endl;
  out << "#  vertices: " << geometry.mesh.nVertices() << endl;
  out << "#     edges: " << geometry.mesh.nEdges() << endl;
  out << "#     faces: " << geometry.mesh.nFaces() << endl;
}

void WavefrontOBJ::writeVertices(std::ofstream& out, EmbeddedGeometryInterface& geometry) {
  SurfaceMesh& mesh(geometry.mesh);
  geometry.requireVertexPositions();

  for (Vertex v : mesh.vertices()) {
    Vector3 p = geometry.vertexPositions[v];
    out << "v " << p.x << " " << p.y << " " << p.z << endl;
  }
}

void WavefrontOBJ::writeTexCoords(std::ofstream& out, EmbeddedGeometryInterface& geometry,
                                  CornerData<Vector2>& texcoords) {
  SurfaceMesh& mesh(geometry.mesh);

  for (Corner c : mesh.corners()) {
    Vector2 z = texcoords[c];
    out << "vt " << z.x << " " << z.y << endl;
  }
}

void WavefrontOBJ::writeNormals(std::ofstream& out, EmbeddedGeometryInterface& geometry, CornerData<Vector3>& normals) {
  SurfaceMesh& mesh(geometry.mesh);

  for (Corner c : mesh.corners()) {
    Vector3 n = normals[c];
    out << "vn " << n.x << " " << n.y << " " << n.z << endl;
  }
}

void WavefrontOBJ::writeFaces(std::ofstream& out, EmbeddedGeometryInterface& geometry, bool useTexCoords,
                              bool useNormals) {
  SurfaceMesh& mesh(geometry.mesh);

  // Get vertex indices
  VertexData<size_t> indices = mesh.getVertexIndices();
  CornerData<size_t> cIndices = mesh.getCornerIndices();

  auto indexFn = [&](Corner c) {
    std::string texCoordString = (useTexCoords) ? std::to_string(cIndices[c] + 1) : "";
    std::string normalString = (useNormals) ? std::to_string(cIndices[c] + 1) : "";
    return " " + std::to_string(indices[c.vertex()] + 1) + "/" + texCoordString + "/" + normalString;
  };

  for (Face f : mesh.faces()) {
    out << "f";
    for (Corner c : f.adjacentCorners()) {
      out << indexFn(c);
    }
    out << endl;
  }
}

std::array<std::pair<std::vector<size_t>, size_t>, 5> polyscopePermutations(SurfaceMesh& mesh) {
  std::array<std::pair<std::vector<size_t>, size_t>, 5> result;

  // This works because of the iteration order that we we know these iterators obey. If iteration orders ever change,
  // this will be broken.

  // Note: vertices & faces are generally in the polyscope default ordering, but this is still useful in hte case of a
  // non-compressed mesh.

  { // Vertices
    std::vector<size_t>& vertPerm = result[0].first;
    vertPerm.resize(mesh.nVertices());
    result[0].second = mesh.nVerticesCapacity();

    size_t i = 0;
    for (Vertex v : mesh.vertices()) {
      vertPerm[i++] = v.getIndex();
    }
  }

  { // Faces
    std::vector<size_t>& facePerm = result[1].first;
    facePerm.resize(mesh.nFaces());
    result[1].second = mesh.nFacesCapacity();

    size_t i = 0;
    for (Face f : mesh.faces()) {
      facePerm[i++] = f.getIndex();
    }
  }

  { // Edges
    std::vector<size_t>& edgePerm = result[2].first;
    edgePerm.resize(mesh.nEdges());
    result[2].second = mesh.nEdgesCapacity();

    EdgeData<char> edgeSeen(mesh, false);
    size_t i = 0;
    for (Face f : mesh.faces()) {
      for (Edge e : f.adjacentEdges()) {
        if (!edgeSeen[e]) {
          edgePerm[i++] = e.getIndex();
          edgeSeen[e] = true;
        }
      }
    }
  }

  { // Halfedges
    std::vector<size_t>& halfedgePerm = result[3].first;
    halfedgePerm.resize(mesh.nInteriorHalfedges());
    result[3].second = mesh.nHalfedgesCapacity();

    size_t i = 0;
    for (Face f : mesh.faces()) {
      for (Halfedge he : f.adjacentHalfedges()) {
        halfedgePerm[i++] = he.getIndex();
      }
    }
  }

  { // Corners
    std::vector<size_t>& cornerPerm = result[4].first;
    cornerPerm.resize(mesh.nCorners());
    result[4].second = mesh.nHalfedgesCapacity();

    size_t i = 0;
    for (Face f : mesh.faces()) {
      for (Corner c : f.adjacentCorners()) {
        cornerPerm[i++] = c.getIndex();
      }
    }
  }

  return result;
}

EdgeData<char> polyscopeEdgeOrientations(SurfaceMesh& mesh) {

  EdgeData<char> result(mesh);
  VertexData<size_t> vertexIndices = mesh.getVertexIndices();

  for (Edge e : mesh.edges()) {
    result[e] = vertexIndices[e.halfedge().vertex()] < vertexIndices[e.halfedge().next().vertex()];
  }

  return result;
}

} // namespace surface
} // namespace geometrycentral
