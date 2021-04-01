#include "geometrycentral/pointcloud/point_cloud_io.h"

#include "geometrycentral/surface/simple_polygon_mesh.h"


#include "happly.h"

namespace geometrycentral {
namespace pointcloud {

// Anonymous helpers
namespace {

std::vector<std::string> supportedPointCloudTypes = {"obj", "ply"};

std::string typeFromFilename(std::string filename) {

  std::string::size_type sepInd = filename.rfind('.');
  std::string type;

  if (sepInd != std::string::npos) {
    std::string extension;
    extension = filename.substr(sepInd + 1);

    // Convert to all lowercase
    std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
    type = extension;
  } else {
    throw std::runtime_error("Could not auto-detect file type to read/write point cloud from " + filename);
  }

  // Check if this is one of the filetypes we're aware of
  if (std::find(std::begin(supportedPointCloudTypes), std::end(supportedPointCloudTypes), type) ==
      std::end(supportedPointCloudTypes)) {
    throw std::runtime_error("Detected file type " + type + " to read/write point cloud from " + filename +
                             ". This is not a supported file type.");
  }

  return type;
}

} // namespace

// === Readers ===

// Particular per-filetype readers called by the general versions below
namespace {

std::tuple<std::unique_ptr<PointCloud>, std::unique_ptr<PointPositionGeometry>> readPointCloud_obj(std::istream& in) {
  // Bootstrap off the mesh reader
  surface::SimplePolygonMesh mesh(in, "obj");

  std::unique_ptr<PointCloud> cloud(new PointCloud(mesh.nVertices()));
  std::unique_ptr<PointPositionGeometry> geom(new PointPositionGeometry(*cloud));
  for (size_t i = 0; i < mesh.nVertices(); i++) {
    geom->positions[i] = mesh.vertexCoordinates[i];
  }

  return std::make_tuple(std::move(cloud), std::move(geom));
}

std::tuple<std::unique_ptr<PointCloud>, std::unique_ptr<PointPositionGeometry>> readPointCloud_ply(std::istream& in) {
  happly::PLYData plyIn(in);

  std::vector<std::array<double, 3>> vPos = plyIn.getVertexPositions();
  size_t N = vPos.size();

  std::unique_ptr<PointCloud> cloud(new PointCloud(N));
  std::unique_ptr<PointPositionGeometry> geom(new PointPositionGeometry(*cloud));
  for (size_t i = 0; i < N; i++) {
    for (int j = 0; j < 3; j++) {
      geom->positions[i][j] = vPos[i][j];
    }
  }

  return std::make_tuple(std::move(cloud), std::move(geom));
}

} // namespace

std::tuple<std::unique_ptr<PointCloud>, std::unique_ptr<PointPositionGeometry>> readPointCloud(std::string filename,
                                                                                               std::string type) {

  // Attempt to detect filename
  bool typeGiven = type != "";
  if (!typeGiven) {
    type = typeFromFilename(filename);
  }

  // == Open the file and load it
  // NOTE: Intentionally always open the stream as binary, even though some of the subsequent formats are plaintext and
  // others are binary.  The only real difference is that non-binary mode performs automatic translation of line ending
  // characters (e.g. \r\n --> \n from DOS). However, this behavior is platform-dependent and having platform-dependent
  // behavior seems more confusing then just handling the newlines properly in the parsers.
  std::ifstream inStream(filename, std::ios::binary);
  if (!inStream) throw std::runtime_error("couldn't open file " + filename);

  return readPointCloud(inStream, type);
}

std::tuple<std::unique_ptr<PointCloud>, std::unique_ptr<PointPositionGeometry>> readPointCloud(std::istream& in,
                                                                                               std::string type) {

  if (type == "obj") {
    return readPointCloud_obj(in);
  } else if (type == "ply") {
    return readPointCloud_ply(in);
  } else {
    throw std::runtime_error("Did not recognize point cloud file type " + type);
  }

  return std::tuple<std::unique_ptr<PointCloud>, std::unique_ptr<PointPositionGeometry>>{nullptr,
                                                                                         nullptr}; // unreachable
}

// === Writers ===

// Particular per-filetype writers called by the general versions below
namespace {

void writePointCloud_obj(PointCloud& cloud, PointPositionGeometry& geometry, std::ostream& out) {
  // Just bootstrap off the mesh version
  surface::SimplePolygonMesh mesh;
  mesh.vertexCoordinates.resize(cloud.nPoints());
  for (size_t i = 0; i < cloud.nPoints(); i++) {
    mesh.vertexCoordinates[i] = geometry.positions[i];
  }

  mesh.writeMesh(out, "obj");
}

void writePointCloud_ply(PointCloud& cloud, PointPositionGeometry& geometry, std::ostream& out) {
  // Copy to a std vector
  std::vector<std::array<double, 3>> vertexPos;
  vertexPos.resize(cloud.nPoints());
  for (size_t i = 0; i < cloud.nPoints(); i++) {
    for (int j = 0; j < 3; j++) {
      vertexPos[i][j] = geometry.positions[i][j];
    }
  }

  happly::PLYData plyOut;
  plyOut.addVertexPositions(vertexPos);
  plyOut.write(out, happly::DataFormat::Binary);
}


} // namespace

void writePointCloud(PointCloud& cloud, PointPositionGeometry& geometry, std::string filename, std::string type) {

  // Attempt to detect filename
  bool typeGiven = type != "";
  if (!typeGiven) {
    type = typeFromFilename(filename);
  }

  // == Open the file and write it
  // NOTE: Intentionally always open the stream as binary, even though some of the subsequent formats are plaintext and
  // others are binary.  The only real difference is that non-binary mode performs automatic translation of line ending
  // characters (e.g. \r\n --> \n from DOS). However, this behavior is platform-dependent and having platform-dependent
  // behavior seems more confusing then just handling the newlines properly in the parsers.
  std::ofstream outStream(filename, std::ios::binary);
  if (!outStream) throw std::runtime_error("couldn't open file " + filename);

  return writePointCloud(cloud, geometry, outStream, type);
}

void writePointCloud(PointCloud& cloud, PointPositionGeometry& geometry, std::ostream& out, std::string type) {

  if (type == "obj") {
    return writePointCloud_obj(cloud, geometry, out);
  } else if (type == "ply") {
    return writePointCloud_ply(cloud, geometry, out);
  } else {
    throw std::runtime_error("Did not recognize point cloud file type " + type);
  }
}

} // namespace pointcloud
} // namespace geometrycentral
