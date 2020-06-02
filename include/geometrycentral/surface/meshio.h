#pragma once

// The MeshIO class provides a variety of methods for mesh input/output.

#include "geometrycentral/surface/vertex_position_geometry.h"

#include <fstream>
#include <string>

namespace geometrycentral {
namespace surface {

// Loads a halfedge mesh and its geometry from file.
// Specify a type like "ply" or "obj", if no type is specified, attempts to infer from extension.
std::tuple<std::unique_ptr<HalfedgeMesh>, std::unique_ptr<VertexPositionGeometry>>
loadMesh(std::string filename, std::string type = "");

class WavefrontOBJ {
public:
  static bool write(std::string filename, EmbeddedGeometryInterface& geometry);
  static bool write(std::string filename, EmbeddedGeometryInterface& geometry, CornerData<Vector2>& texcoords);

protected:
  static bool openStream(std::ofstream& out, std::string filename);
  static void writeHeader(std::ofstream& out, EmbeddedGeometryInterface& geometry);
  static void writeVertices(std::ofstream& out, EmbeddedGeometryInterface& geometry);
  static void writeTexCoords(std::ofstream& out, EmbeddedGeometryInterface& geometry, CornerData<Vector2>& texcoords);
  static void writeFaces(std::ofstream& out, EmbeddedGeometryInterface& geometry, bool useTexCoords = false);
};


// TODO write halfedge mesh as a permutation, in binary format (for quicker loading/smaller files)


// === Integrations with other libraries and formats

// Generate the permutations which map geometry-central's indexing order to the Polyscope indexing order
std::array<std::pair<std::vector<size_t>, size_t>, 5> polyscopePermutations(HalfedgeMesh& mesh);

// Generate an booleans to communicate the canonical orientations of edges for 1form-valued data
EdgeData<char> polyscopeEdgeOrientations(HalfedgeMesh& mesh);


} // namespace surface
} // namespace geometrycentral
