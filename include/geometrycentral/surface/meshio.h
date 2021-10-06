#pragma once

// The MeshIO class provides a variety of methods for mesh input/output.

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include <fstream>
#include <string>

namespace geometrycentral {
namespace surface {

// Loads a halfedge mesh and its geometry from file.
// Specify a type like "ply" or "obj", if no type is specified, attempts to infer from extension.

// Load a general surface mesh, which might or might not be manifold
std::tuple<std::unique_ptr<SurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
readSurfaceMesh(std::string filename, std::string type = "");
std::tuple<std::unique_ptr<SurfaceMesh>, std::unique_ptr<VertexPositionGeometry>> readSurfaceMesh(std::istream& in,
                                                                                                  std::string type);

// Load a manifold surface mesh; an exception will by thrown if the mesh is not manifold.
std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
readManifoldSurfaceMesh(std::string filename, std::string type = "");
std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
readManifoldSurfaceMesh(std::istream& in, std::string type);

// Load a mesh with UV coordinates, which will be stored as data at triangle
// corners (to allow for UVs that are discontinuous across edges, e.g., at cuts)
std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>,
           std::unique_ptr<CornerData<Vector2>>>
readParameterizedManifoldSurfaceMesh(std::string filename, std::string type = "");
std::tuple<std::unique_ptr<SurfaceMesh>, std::unique_ptr<VertexPositionGeometry>, std::unique_ptr<CornerData<Vector2>>>
readParameterizedSurfaceMesh(std::string filename, std::string type = "");

// Legacy method, prefer one of the variants above
std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
loadMesh(std::string filename, std::string type = "");


// Write a surface mesh
void writeSurfaceMesh(SurfaceMesh& mesh, EmbeddedGeometryInterface& geometry, std::string filename,
                      std::string type = "");
void writeSurfaceMesh(SurfaceMesh& mesh, EmbeddedGeometryInterface& geometry, CornerData<Vector2>& texCoords,
                      std::string filename, std::string type = "");
void writeSurfaceMesh(SurfaceMesh& mesh, EmbeddedGeometryInterface& geometry, std::ostream& out, std::string type);
void writeSurfaceMesh(SurfaceMesh& mesh, EmbeddedGeometryInterface& geometry, CornerData<Vector2>& texCoords,
                      std::ostream& out, std::string type);

// Helpers to map vertex data
CornerData<Vector2> packToParam(SurfaceMesh& mesh, VertexData<double>& vals); // 0 in Y coord
CornerData<Vector2> packToParam(SurfaceMesh& mesh, VertexData<double>& valsX, VertexData<double>& valsY);

// Legacy writer, prefer on of the variants above
class WavefrontOBJ {
public:
  static bool write(std::string filename, EmbeddedGeometryInterface& geometry);
  static bool write(std::string filename, EmbeddedGeometryInterface& geometry, CornerData<Vector2>& texcoords);
  static bool write(std::string filename, EmbeddedGeometryInterface& geometry, CornerData<Vector3>& normals);
  static bool write(std::string filename, EmbeddedGeometryInterface& geometry, CornerData<Vector2>& texcoords,
                    CornerData<Vector3>& normals);

protected:
  static bool openStream(std::ofstream& out, std::string filename);
  static void writeHeader(std::ofstream& out, EmbeddedGeometryInterface& geometry);
  static void writeVertices(std::ofstream& out, EmbeddedGeometryInterface& geometry);
  static void writeTexCoords(std::ofstream& out, EmbeddedGeometryInterface& geometry, CornerData<Vector2>& texcoords);
  static void writeNormals(std::ofstream& out, EmbeddedGeometryInterface& geometry, CornerData<Vector3>& normals);
  static void writeFaces(std::ofstream& out, EmbeddedGeometryInterface& geometry, bool useTexCoords = false,
                         bool useNormals = false);
};


// TODO write halfedge mesh as a permutation, in binary format (for quicker loading/smaller files)


// === Integrations with other libraries and formats

// Generate the permutations which map geometry-central's indexing order to the Polyscope indexing order
std::array<std::pair<std::vector<size_t>, size_t>, 5> polyscopePermutations(SurfaceMesh& mesh);

// Generate booleans to communicate the canonical orientations of edges for 1form-valued data
EdgeData<char> polyscopeEdgeOrientations(SurfaceMesh& mesh);


} // namespace surface
} // namespace geometrycentral
