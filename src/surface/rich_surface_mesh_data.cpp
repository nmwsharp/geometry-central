#include "geometrycentral/surface/rich_surface_mesh_data.h"

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


RichSurfaceMeshData::RichSurfaceMeshData(SurfaceMesh& mesh_, std::string filename_)
    : plyData(filename_), mesh(&mesh_) {}

RichSurfaceMeshData::RichSurfaceMeshData(SurfaceMesh& mesh_, std::istream& in_) : plyData(in_), mesh(&mesh_) {}

RichSurfaceMeshData::RichSurfaceMeshData(SurfaceMesh& mesh_) : plyData(), mesh(&mesh_) {}

RichSurfaceMeshData::RichSurfaceMeshData(std::istream& in_) : plyData(in_) { loadMeshFromFile(); }
RichSurfaceMeshData::RichSurfaceMeshData(std::string filename_) : plyData(filename_) { loadMeshFromFile(); }

void RichSurfaceMeshData::loadMeshFromFile() {
  if (mesh != nullptr) throw std::runtime_error("cannot load mesh multiple times");

  if (!plyData.hasElement("gc_internal_halfedge")) {
    throw std::runtime_error(
        "Cannot load mesh from file, file was saved without connectiviy. Call addMeshConnectivity() before saving.");
  }

  // Annoyingly, ply doesn't allow uint64_t (aka size_t on most systems)... see note in addMeshConnectivity()
  auto fromSmallerVec = [](const std::vector<uint32_t>& vec) {
    std::vector<size_t> out(vec.size());
    for (size_t i = 0; i < vec.size(); i++) {
      size_t val = vec[i];
      if (val == std::numeric_limits<uint32_t>::max()) val = INVALID_IND;
      out[i] = static_cast<size_t>(val);
    }
    return out;
  };

  // clang-format off

  // Read the necessary arrays
  std::vector<size_t> heNextArr     = fromSmallerVec(plyData.getElement("gc_internal_halfedge").getProperty<uint32_t>("gc_internal_heNextArr"));
  std::vector<size_t> heVertexArr   = fromSmallerVec(plyData.getElement("gc_internal_halfedge").getProperty<uint32_t>("gc_internal_heVertexArr"));
  std::vector<size_t> heFaceArr     = fromSmallerVec(plyData.getElement("gc_internal_halfedge").getProperty<uint32_t>("gc_internal_heFaceArr"));
  std::vector<size_t> vHalfedgeArr  = fromSmallerVec(plyData.getElement("gc_internal_vertex").getProperty<uint32_t>("gc_internal_vHalfedgeArr"));
  std::vector<size_t> fHalfedgeArr  = fromSmallerVec(plyData.getElement("gc_internal_face").getProperty<uint32_t>("gc_internal_fHalfedgeArr"));
  std::vector<size_t> fHalfedgeArrB = fromSmallerVec(plyData.getElement("gc_internal_bl").getProperty<uint32_t>("gc_internal_blHalfedgeArr"));
  fHalfedgeArr.insert(fHalfedgeArr.end(), fHalfedgeArrB.begin(), fHalfedgeArrB.end());


  std::vector<size_t> heSiblingArr; 
  std::vector<size_t> heEdgeArr;   
  std::vector<char> heOrientArr;   
  std::vector<size_t> eHalfedgeArr;
  bool useImplicitTwin = !plyData.getElement("gc_internal_halfedge").hasProperty("gc_internal_heSiblingArr");
  if(!useImplicitTwin) {
    heSiblingArr    = fromSmallerVec(plyData.getElement("gc_internal_halfedge").getProperty<uint32_t>("gc_internal_heSiblingArr"));
    heEdgeArr       = fromSmallerVec(plyData.getElement("gc_internal_halfedge").getProperty<uint32_t>("gc_internal_heEdgeArr"));
    heOrientArr     = plyData.getElement("gc_internal_halfedge").getProperty<char>("gc_internal_heOrientArr");
    eHalfedgeArr    = fromSmallerVec(plyData.getElement("gc_internal_edge").getProperty<uint32_t>("gc_internal_eHalfedgeArr"));
  }

  // clang-format on

  // Build the actual mesh
  if (useImplicitTwin) {
    mesh = new ManifoldSurfaceMesh(heNextArr, heVertexArr, heFaceArr, vHalfedgeArr, fHalfedgeArr, fHalfedgeArrB.size());
  } else {
    mesh = new SurfaceMesh(heNextArr, heVertexArr, heFaceArr, vHalfedgeArr, fHalfedgeArr, heSiblingArr, heEdgeArr,
                           heOrientArr, eHalfedgeArr, fHalfedgeArrB.size());
  }
}

std::tuple<std::unique_ptr<SurfaceMesh>, std::unique_ptr<RichSurfaceMeshData>>
RichSurfaceMeshData::readMeshAndData(std::string filename) {
  std::unique_ptr<RichSurfaceMeshData> data(new RichSurfaceMeshData(filename));
  return std::make_tuple(std::unique_ptr<SurfaceMesh>(data->mesh), std::move(data));
}
std::tuple<std::unique_ptr<SurfaceMesh>, std::unique_ptr<RichSurfaceMeshData>>
RichSurfaceMeshData::readMeshAndData(std::istream& in) {
  std::unique_ptr<RichSurfaceMeshData> data(new RichSurfaceMeshData(in));
  return std::make_tuple(std::unique_ptr<SurfaceMesh>(data->mesh), std::move(data));
}
std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<RichSurfaceMeshData>>
RichSurfaceMeshData::readManifoldMeshAndData(std::string filename) {
  std::unique_ptr<RichSurfaceMeshData> data(new RichSurfaceMeshData(filename));
  ManifoldSurfaceMesh* manifMesh = dynamic_cast<ManifoldSurfaceMesh*>(data->mesh);
  if (manifMesh == nullptr) {
    throw std::runtime_error("tried to read ManifoldSurfaceMesh, but file was not written from a ManifoldSurfaceMesh");
  }
  return std::make_tuple(std::unique_ptr<ManifoldSurfaceMesh>(manifMesh), std::move(data));
}
std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<RichSurfaceMeshData>>
RichSurfaceMeshData::readManifoldMeshAndData(std::istream& in) {
  std::unique_ptr<RichSurfaceMeshData> data(new RichSurfaceMeshData(in));
  ManifoldSurfaceMesh* manifMesh = dynamic_cast<ManifoldSurfaceMesh*>(data->mesh);
  if (manifMesh == nullptr) {
    throw std::runtime_error("tried to read ManifoldSurfaceMesh, but file was not written from a ManifoldSurfaceMesh");
  }
  return std::make_tuple(std::unique_ptr<ManifoldSurfaceMesh>(manifMesh), std::move(data));
}

void RichSurfaceMeshData::addMeshConnectivity() {

  // I _think_ the code below does not require compressed
  // GC_SAFETY_ASSERT(mesh.isCompressed(), "Mesh must be compressed to use RichSurfaceMeshData. Call mesh.compress().");

  // == Write connectiviy as indices
  std::vector<std::vector<size_t>> faceIndices = mesh->getFaceVertexList();
  plyData.addFaceIndices(faceIndices);

  // == Write rich connectivity
  // These needed to be added on distinct elements, since in an uncompressed mesh the length of the arrays might be
  // different (also avoids name collisions)

  // Annoyingly, ply doesn't allow uint64_t (aka size_t on most systems).  It does allow uint32_t, so pack to one of
  // those... This will fail on sufficiently large data, and presumes INVALID_IND is the only really-big value used.
  auto toSmallerVec = [](std::vector<size_t>::iterator b, std::vector<size_t>::iterator e) {
    size_t count = std::distance(b, e);
    std::vector<uint32_t> out(count);
    for (size_t i = 0; i < count; i++) {
      size_t val = *b;
      if (val == INVALID_IND) val = std::numeric_limits<uint32_t>::max();
      out[i] = static_cast<uint32_t>(val);
      b++;
    }
    return out;
  };

  // clang-format off

  // Add the elements themselves
  plyData.addElement("gc_internal_vertex", mesh->nVerticesFillCount);
  plyData.addElement("gc_internal_halfedge", mesh->nHalfedgesFillCount);
  plyData.addElement("gc_internal_edge", mesh->nEdgesFillCount);
  plyData.addElement("gc_internal_face", mesh->nFacesFillCount);
  plyData.addElement("gc_internal_bl", mesh->nBoundaryLoopsFillCount);

  auto& elemVertex = plyData.getElement("gc_internal_vertex");
  auto& elemHalfedge = plyData.getElement("gc_internal_halfedge");
  auto& elemEdge = plyData.getElement("gc_internal_edge");
  auto& elemFace = plyData.getElement("gc_internal_face");
  auto& elemBl= plyData.getElement("gc_internal_bl");

  // Halfedge properties
  elemHalfedge.addProperty<uint32_t>("gc_internal_heNextArr", toSmallerVec(mesh->heNextArr.begin(), mesh->heNextArr.begin() + mesh->nHalfedgesFillCount));
  elemHalfedge.addProperty<uint32_t>("gc_internal_heFaceArr", toSmallerVec(mesh->heFaceArr.begin(), mesh->heFaceArr.begin() + mesh->nHalfedgesFillCount));
  elemHalfedge.addProperty<uint32_t>("gc_internal_heVertexArr", toSmallerVec(mesh->heVertexArr.begin(), mesh->heVertexArr.begin() + mesh->nHalfedgesFillCount));
  if(!mesh->usesImplicitTwin()) {
    elemHalfedge.addProperty<uint32_t>("gc_internal_heSiblingArr", toSmallerVec(mesh->heSiblingArr.begin(), mesh->heSiblingArr.begin() + mesh->nHalfedgesFillCount));
    elemHalfedge.addProperty<uint32_t>("gc_internal_heEdgeArr", toSmallerVec(mesh->heEdgeArr.begin(), mesh->heEdgeArr.begin() + mesh->nHalfedgesFillCount));
    elemHalfedge.addProperty<char>("gc_internal_heOrientArr", std::vector<char>(mesh->heOrientArr.begin(), mesh->heOrientArr.begin() + mesh->nHalfedgesFillCount));
  }

  // Vertex properties
  elemVertex.addProperty<uint32_t>("gc_internal_vHalfedgeArr", toSmallerVec(mesh->vHalfedgeArr.begin(), mesh->vHalfedgeArr.begin() + mesh->nVerticesFillCount));
 
  // Edge properties
  if(!mesh->usesImplicitTwin()) {
    elemEdge.addProperty<uint32_t>("gc_internal_eHalfedgeArr", toSmallerVec(mesh->eHalfedgeArr.begin(), mesh->eHalfedgeArr.begin() + mesh->nEdgesFillCount));
  }

  // Face properties
  elemFace.addProperty<uint32_t>("gc_internal_fHalfedgeArr", toSmallerVec(mesh->fHalfedgeArr.begin(), mesh->fHalfedgeArr.begin() + mesh->nFacesFillCount));
  
  // Boundary loop properties
  elemBl.addProperty<uint32_t>("gc_internal_blHalfedgeArr", toSmallerVec(mesh->fHalfedgeArr.end() - mesh->nBoundaryLoopsFillCount, mesh->fHalfedgeArr.end()));

  
  /*
    std::vector<size_t> heNextArr;                // he.next()
  std::vector<size_t> heVertexArr;  // he.vertex()
  std::vector<size_t> heFaceArr;    // he.face()
  std::vector<size_t> vHalfedgeArr; // v.halfedge()
  std::vector<size_t> fHalfedgeArr; // f.halfedge()
  // const bool useImplicitTwinFlag;
  std::vector<size_t> heSiblingArr; // he.sibling() and he.twin()
  std::vector<size_t> heEdgeArr;    // he.edge()
  std::vector<size_t> eHalfedgeArr; // e.halfedge()
  */

  // clang-format on
}

void RichSurfaceMeshData::addGeometry(EmbeddedGeometryInterface& geometry) {

  geometry.requireVertexPositions();

  // separate x/y/z coordinates
  VertexData<double> x(*mesh);
  VertexData<double> y(*mesh);
  VertexData<double> z(*mesh);

  for (Vertex v : mesh->vertices()) {
    Vector3 p = geometry.vertexPositions[v];
    x[v] = p.x;
    y[v] = p.y;
    z[v] = p.z;
  }

  addVertexProperty("x", x);
  addVertexProperty("y", y);
  addVertexProperty("z", z);

  geometry.unrequireVertexPositions();
}


void RichSurfaceMeshData::addIntrinsicGeometry(IntrinsicGeometryInterface& geometry) {
  geometry.requireEdgeLengths();
  addEdgeProperty("intrinsic_edge_lengths", geometry.edgeLengths);
  geometry.unrequireEdgeLengths();
}


std::unique_ptr<VertexPositionGeometry> RichSurfaceMeshData::getGeometry() {

  // Get x/y/z coordinates
  VertexData<double> x = getVertexProperty<double>("x");
  VertexData<double> y = getVertexProperty<double>("y");
  VertexData<double> z = getVertexProperty<double>("z");

  VertexData<Vector3> positions(*mesh);
  for (Vertex v : mesh->vertices()) {
    positions[v] = Vector3{x[v], y[v], z[v]};
  }

  // Return a new geometry object
  std::unique_ptr<VertexPositionGeometry> geom(new VertexPositionGeometry(*mesh, positions));
  return geom;
}

std::unique_ptr<EdgeLengthGeometry> RichSurfaceMeshData::getIntrinsicGeometry() {
  EdgeData<double> eLengths = getEdgeProperty<double>("intrinsic_edge_lengths");
  return std::unique_ptr<EdgeLengthGeometry>(new EdgeLengthGeometry(*mesh, eLengths));
}


VertexData<Vector3> RichSurfaceMeshData::getVertexColors() {

  VertexData<Vector3> color(*mesh);

  try {
    // Try uchar first
    VertexData<unsigned char> r = getVertexProperty<unsigned char>("red");
    VertexData<unsigned char> g = getVertexProperty<unsigned char>("green");
    VertexData<unsigned char> b = getVertexProperty<unsigned char>("blue");
    for (Vertex v : mesh->vertices()) {
      color[v][0] = r[v] / 255.0;
      color[v][1] = g[v] / 255.0;
      color[v][2] = b[v] / 255.0;
    }
    return color;

  } catch (std::runtime_error &orig_e) {

    // If that doesn't work, try float
    try {
      VertexData<double> r = getVertexProperty<double>("red");
      VertexData<double> g = getVertexProperty<double>("green");
      VertexData<double> b = getVertexProperty<double>("blue");
      for (Vertex v : mesh->vertices()) {
        color[v][0] = r[v];
        color[v][1] = g[v];
        color[v][2] = b[v];
      }
      return color;
    } catch (std::runtime_error &second_e) {
      throw std::runtime_error("Could not find vertex colors in PLY file, as uchar or float");
    }
  }
}

void RichSurfaceMeshData::write(std::string filename) { plyData.write(filename, outputFormat); }
void RichSurfaceMeshData::write(std::ostream& out) { plyData.write(out, outputFormat); }

} // namespace surface
} // namespace geometrycentral
