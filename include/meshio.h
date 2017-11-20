#pragma once

// The MeshIO class provides a variety of methods for mesh input/output.

#include <string>
#include <fstream>
#include <geometry.h>

class WavefrontOBJ {
 public:
  static bool write(std::string filename, Geometry<Euclidean>& geometry);
  static bool write(std::string filename, Geometry<Euclidean>& geometry,
                    CornerData<Vector2>& texcoords);

 protected:
  static bool openStream(std::ofstream& out, std::string filename);
  static void writeHeader(std::ofstream& out, Geometry<Euclidean>& geometry);
  static void writeVertices(std::ofstream& out, Geometry<Euclidean>& geometry);
  static void writeTexCoords(std::ofstream& out, Geometry<Euclidean>& geometry,
                             CornerData<Vector2>& texcoords);
  static void writeFaces(std::ofstream& out, Geometry<Euclidean>& geometry,
                         bool useTexCoords = false);
};

class PLY {
 public:
  static bool write(std::string filename, Geometry<Euclidean>& geometry,
                    VertexData<Vector3> colors);
};

// TODO write halfedge mesh as a permutation, in binary format (for quicker
// loading/smaller files)
