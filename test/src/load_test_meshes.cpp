#include "load_test_meshes.h"

#include "geometrycentral/surface/meshio.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

using std::cout;
using std::endl;

// Helpers
namespace {
std::string guessNiceNameFromPath(std::string fullname) {
  size_t startInd = 0;
  for (std::string sep : {"/", "\\"}) {
    size_t pos = fullname.rfind(sep);
    if (pos != std::string::npos) {
      startInd = std::max(startInd, pos + 1);
    }
  }

  return fullname.substr(startInd, fullname.size() - startInd);
};
} // namespace

MeshAsset::MeshAsset(std::string localPath) {
  name = guessNiceNameFromPath(localPath);
  std::string fullPath = std::string(GC_TEST_ASSETS_ABS_PATH) + "/" + name;
  cout << "  -- info: Loading mesh asset " << name << " from " << fullPath << endl;
  sourcePath = fullPath;
  std::tie(mesh, geometry) = loadMesh(fullPath);

  hasBoundary = mesh->hasBoundary();
  isTriangular = mesh->isTriangular();
}

MeshAsset MeshAsset::copy() const {
  MeshAsset newM;

  newM.name = name;
  newM.sourcePath = sourcePath;
  newM.mesh = mesh->copy();
  if (geometry) {
     newM.geometry = geometry->reinterpretTo(*newM.mesh);
  }

  newM.hasBoundary = hasBoundary;
  newM.isTriangular = isTriangular;
  newM.isPolygonalComplex = isPolygonalComplex;

  return newM;
}

void MeshAsset::printThyName() { cout << "  testing on mesh: " << name << endl; }

// Static storage for mesh assets
std::vector<MeshAsset> MeshAssetSuite::allMeshAssets;

void MeshAssetSuite::SetUpTestSuite() {
  // no need to set up more than once
  if(allMeshAssets.size() > 0) return;


  allMeshAssets.emplace_back("tet.obj");
  allMeshAssets.emplace_back("spot.ply");
  allMeshAssets.emplace_back("bob_small.ply");
  allMeshAssets.emplace_back("sphere_small.ply");
  allMeshAssets.emplace_back("lego.ply");
  allMeshAssets.emplace_back("dodecahedron_poly.obj");
  allMeshAssets.emplace_back("platonic_shelf.obj");
}


MeshAsset MeshAssetSuite::getAsset(std::string name) {
  for (MeshAsset& a : allMeshAssets) {
    if (a.name == name) {
      return a.copy();
    }
  }
  throw std::runtime_error("no mesh asset named " + name);
}

std::vector<MeshAsset> MeshAssetSuite::allMeshes(bool includeNoGeom) {
  std::vector<MeshAsset> result;
  for (MeshAsset& a : allMeshAssets) {
    if (includeNoGeom || a.geometry) result.push_back(a.copy());
  }
  return result;
}

std::vector<MeshAsset> MeshAssetSuite::closedMeshes(bool includeNoGeom) {
  std::vector<MeshAsset> result;
  for (MeshAsset& a : allMeshAssets) {
    if (!a.hasBoundary) {
      if (includeNoGeom || a.geometry) result.push_back(a.copy());
    }
  }
  return result;
}

std::vector<MeshAsset> MeshAssetSuite::boundaryMeshes(bool includeNoGeom) {
  std::vector<MeshAsset> result;
  for (MeshAsset& a : allMeshAssets) {
    if (a.hasBoundary) {
      if (includeNoGeom || a.geometry) result.push_back(a.copy());
    }
  }
  return result;
}

std::vector<MeshAsset> MeshAssetSuite::triangularMeshes(bool includeNoGeom) {
  std::vector<MeshAsset> result;
  for (MeshAsset& a : allMeshAssets) {
    if (a.isTriangular) {
      if (includeNoGeom || a.geometry) result.push_back(a.copy());
    }
  }
  return result;
}
