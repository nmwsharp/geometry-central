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

MeshAsset::MeshAsset(std::string localPath, bool loadManifold) {
  name = guessNiceNameFromPath(localPath);
  std::string fullPath = std::string(GC_TEST_ASSETS_ABS_PATH) + "/" + name;
  cout << "  -- info: Loading mesh asset " << name << " from " << fullPath << endl;
  sourcePath = fullPath;

  if (loadManifold) {
    std::unique_ptr<ManifoldSurfaceMesh> manifMesh;
    std::tie(manifMesh, geometry) = readManifoldSurfaceMesh(fullPath);
    manifoldMesh = manifMesh.release();
    mesh.reset(manifoldMesh);
    isSubclassManifoldSurfaceMesh = true;
  } else {
    std::tie(mesh, geometry) = readSurfaceMesh(fullPath);
    isSubclassManifoldSurfaceMesh = false;
  }

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

  newM.isSubclassManifoldSurfaceMesh = isSubclassManifoldSurfaceMesh;
  if (isSubclassManifoldSurfaceMesh) {
    newM.manifoldMesh = dynamic_cast<ManifoldSurfaceMesh*>(newM.mesh.get());
  }

  newM.hasBoundary = hasBoundary;
  newM.isTriangular = isTriangular;
  newM.isPolygonalComplex = isPolygonalComplex;

  return newM;
}

void MeshAsset::printThyName() const { cout << "  testing on mesh (" << (isSubclassManifoldSurfaceMesh ? "manifold" : "general") <<  "): " << name << endl; }

// Static storage for mesh assets
std::vector<MeshAsset> MeshAssetSuite::allMeshAssets;

void MeshAssetSuite::SetUpTestSuite() {
  // no need to set up more than once
  if (allMeshAssets.size() > 0) return;


  // Load manifold surface mesh variants
  allMeshAssets.emplace_back("tet.obj", true);
  allMeshAssets.emplace_back("spot.ply", true);
  allMeshAssets.emplace_back("bob_small.ply", true);
  allMeshAssets.emplace_back("sphere_small.ply", true);
  allMeshAssets.emplace_back("lego.ply", true);
  allMeshAssets.emplace_back("dodecahedron_poly.obj", true);
  allMeshAssets.emplace_back("platonic_shelf.obj", true);
  allMeshAssets.emplace_back("fox.ply", true);
  allMeshAssets.emplace_back("cat_head.obj", true);

  // Load general surface mesh variants
  allMeshAssets.emplace_back("tet.obj", false);
  allMeshAssets.emplace_back("spot.ply", false);
  allMeshAssets.emplace_back("bob_small.ply", false);
  allMeshAssets.emplace_back("sphere_small.ply", false);
  allMeshAssets.emplace_back("lego.ply", false);
  allMeshAssets.emplace_back("dodecahedron_poly.obj", false);
  allMeshAssets.emplace_back("platonic_shelf.obj", false);
  allMeshAssets.emplace_back("fox.ply", false);
  allMeshAssets.emplace_back("cat_head.obj", false);

  // Load nonmanifold models
  allMeshAssets.emplace_back("fan3.obj", false);
  allMeshAssets.emplace_back("hourglass_ico.obj", false);
  allMeshAssets.emplace_back("triple_vierbein.obj", false);
  allMeshAssets.emplace_back("moebius.obj", false);
}


MeshAsset MeshAssetSuite::getAsset(std::string name, bool loadManifold) {
  for (MeshAsset& a : allMeshAssets) {
    if (a.name == name && a.isSubclassManifoldSurfaceMesh == loadManifold) {
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

std::vector<MeshAsset> MeshAssetSuite::manifoldSurfaceMeshes(bool includeNoGeom) {
  std::vector<MeshAsset> result;
  for (MeshAsset& a : allMeshAssets) {
    if (a.isSubclassManifoldSurfaceMesh) {
      if (includeNoGeom || a.geometry) result.push_back(a.copy());
    }
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

std::vector<MeshAsset> MeshAssetSuite::polygonalComplexMeshes(bool includeNoGeom) {
  std::vector<MeshAsset> result;
  for (MeshAsset& a : allMeshAssets) {
    if (a.isPolygonalComplex) {
      if (includeNoGeom || a.geometry) result.push_back(a.copy());
    }
  }
  return result;
}

std::vector<MeshAsset> MeshAssetSuite::triangularMeshes(bool includeNoGeom, bool includeNonmanifold) {
  std::vector<MeshAsset> result;
  for (MeshAsset& a : allMeshAssets) {
    if (a.isTriangular) {
      if (!includeNoGeom && !a.geometry) continue;
      if (!includeNonmanifold && !a.isSubclassManifoldSurfaceMesh) continue;
      result.push_back(a.copy());
    }
  }
  return result;
}
