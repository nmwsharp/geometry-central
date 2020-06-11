#include "geometrycentral/surface/ply_surface_mesh_data.h"

#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/polygon_soup_mesh.h"

#include <cstring>
#include <fstream>
#include <iostream>


using std::cerr;
using std::cout;
using std::endl;

namespace geometrycentral {
namespace surface {


PlySurfaceMeshData::PlySurfaceMeshData(SurfaceMesh& mesh_, std::string filename_)
    : mesh(mesh_), plyData(filename_) {}

PlySurfaceMeshData::PlySurfaceMeshData(SurfaceMesh& mesh_)
    : mesh(mesh_), plyData() {
  // Write connectiviy as indices
  std::vector<std::vector<size_t>> faceIndices = mesh.getFaceVertexList();
  plyData.addFaceIndices(faceIndices);
}

std::tuple<std::unique_ptr<SurfaceMesh>, std::unique_ptr<PlySurfaceMeshData>>
PlySurfaceMeshData::loadMeshAndData(std::string filename) {

  // First, just load the connectivity
  SimplePolygonMesh simpleMesh(filename, "ply");
  std::unique_ptr<SurfaceMesh> mesh(new SurfaceMesh(simpleMesh.polygons));

  // Now, open file on that mesh
  std::unique_ptr<PlySurfaceMeshData> data(new PlySurfaceMeshData(*mesh, filename));

  return std::make_tuple(std::move(mesh), std::move(data));
}

void PlySurfaceMeshData::addGeometry(EmbeddedGeometryInterface& geometry) {

  geometry.requireVertexPositions();

  // separate x/y/z coordinates
  VertexData<double> x(mesh);
  VertexData<double> y(mesh);
  VertexData<double> z(mesh);

  for (Vertex v : mesh.vertices()) {
    Vector3 p = geometry.vertexPositions[v];
    x[v] = p.x;
    y[v] = p.y;
    z[v] = p.z;
  }

  addVertexProperty("x", x);
  addVertexProperty("y", y);
  addVertexProperty("z", z);
}


std::unique_ptr<VertexPositionGeometry> PlySurfaceMeshData::getGeometry() {

  // Get x/y/z coordinates
  VertexData<double> x = getVertexProperty<double>("x");
  VertexData<double> y = getVertexProperty<double>("y");
  VertexData<double> z = getVertexProperty<double>("z");

  VertexData<Vector3> positions(mesh);
  for (Vertex v : mesh.vertices()) {
    positions[v] = Vector3{x[v], y[v], z[v]};
  }

  // Return a new geometry object
  std::unique_ptr<VertexPositionGeometry> geom(new VertexPositionGeometry(mesh, positions));
  return geom;
}


VertexData<Vector3> PlySurfaceMeshData::getVertexColors() {

  VertexData<Vector3> color(mesh);

  try {
    // Try uchar first
    VertexData<unsigned char> r = getVertexProperty<unsigned char>("red");
    VertexData<unsigned char> g = getVertexProperty<unsigned char>("green");
    VertexData<unsigned char> b = getVertexProperty<unsigned char>("blue");
    for (Vertex v : mesh.vertices()) {
      color[v][0] = r[v] / 255.0;
      color[v][1] = g[v] / 255.0;
      color[v][2] = b[v] / 255.0;
    }
    return color;

  } catch (std::runtime_error orig_e) {

    // If that doesn't work, try float
    try {
      VertexData<double> r = getVertexProperty<double>("red");
      VertexData<double> g = getVertexProperty<double>("green");
      VertexData<double> b = getVertexProperty<double>("blue");
      for (Vertex v : mesh.vertices()) {
        color[v][0] = r[v];
        color[v][1] = g[v];
        color[v][2] = b[v];
      }
      return color;
    } catch (std::runtime_error second_e) {
      throw std::runtime_error("Could not find vertex colors in PLY file, as uchar or float");
    }
  }
}

void PlySurfaceMeshData::write(std::string filename) { plyData.write(filename, outputFormat); }

} // namespace surface
} // namespace geometrycentral
