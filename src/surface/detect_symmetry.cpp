#include "geometrycentral/surface/detect_symmetry.h"

#include "nanoflann/KDTreeVectorOfVectorsAdaptor.h"
#include "nanoflann/nanoflann.hpp"

#include <array>
#include <vector>

// Interal implementations that hide the NN lookup while allowing it to be
// shared

using std::cout;
using std::endl;

namespace geometrycentral {
namespace surface {

namespace {

// Stupid nanoflann wrapper
typedef KDTreeVectorOfVectorsAdaptor<std::vector<std::vector<double>>, double> KdTree;
KdTree* buildKDTree(Geometry<Euclidean>* geom) {
  HalfedgeMesh* mesh = geom->getMesh();

  // Pack data in a vector of vectors
  std::vector<std::vector<double>>* pts = new std::vector<std::vector<double>>(mesh->nVertices());
  for (size_t i = 0; i < mesh->nVertices(); i++) {
    Vector3 p = geom->position(mesh->vertex(i));
    (*pts)[i] = {p.x, p.y, p.z};
  }

  KdTree* tree = new KdTree(3, *pts);
  tree->index->buildIndex();

  return tree;
}

bool findPoint(KdTree* tree, Vector3 target, double toleranceRadius, size_t& result) {
  std::vector<size_t> ret_indexes(1);
  std::vector<double> out_dists_sqr(1);
  nanoflann::KNNResultSet<double> resultSet(1);
  resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);

  std::vector<double> query = {target.x, target.y, target.z};

  bool success = tree->index->findNeighbors(resultSet, &query[0], nanoflann::SearchParams());

  // Nothing found
  if (!success) return false;

  double dist = std::sqrt(out_dists_sqr[0]);
  if (dist > toleranceRadius) return false; // too far

  // Point found
  result = ret_indexes[0];
  return true;
}

SymmetryResult detectSymmetryMirror(Geometry<Euclidean>* geom, Vector3 planeNormal, Vector3 planePoint, KdTree* tree) {
  HalfedgeMesh* mesh = geom->getMesh();
  planeNormal = unit(planeNormal);
  double toleranceRadius = geom->lengthScale() * 1e-5;

  SymmetryResult result;
  result.symmetryFound = false;
  result.symmetrySet = VertexData<std::vector<Vertex>>(*mesh, std::vector<Vertex>());

  for (Vertex v : mesh->vertices()) {
    // Compute the symmetric point
    Vector3 pos = geom->position(v);
    Vector3 vecToPlane = dot(planeNormal, planePoint - pos) * planeNormal;
    Vector3 mirrorPos = pos + 2 * vecToPlane;

    // If this point is on the positive side of the plane, it's canonical
    bool isCanonical = dot(planeNormal, pos - planePoint) > -toleranceRadius; // small tolerance for points on plane
    if (isCanonical) {
      result.canonicalVertices.push_back(v);
    }

    // If this point is its own pair, there's no mirror to look for (assumes no
    // duplicate verts)
    if (norm(pos - mirrorPos) < toleranceRadius) continue;

    // Search for the point
    size_t mirrorInd;
    bool success = findPoint(tree, mirrorPos, toleranceRadius, mirrorInd);
    if (!success) {
      return result;
    }

    // If found, add to lists
    if (isCanonical) {
      result.symmetrySet[v].push_back(mesh->vertex(mirrorInd));
    }
  }

  result.symmetryFound = true;
  return result;
}

SymmetryResult detectSymmetryRotation(Geometry<Euclidean>* geom, Vector3 rotAxis, Vector3 rotPoint, int nSym,
                                      KdTree* tree) {
  HalfedgeMesh* mesh = geom->getMesh();
  rotAxis = unit(rotAxis);
  double toleranceRadius = geom->lengthScale() * 1e-5;
  double deltaTheta = 2 * PI / nSym;

  // Any axis orthogonal to the rotation axis
  Vector3 castAxis = Vector3{0.12345623, -.883034723, 0.54457119}; // provably random
  Vector3 orthAxis = unit(cross(rotAxis, castAxis));

  SymmetryResult result;
  result.symmetryFound = false;
  result.symmetrySet = VertexData<std::vector<Vertex>>(*mesh, std::vector<Vertex>());

  for (Vertex v : mesh->vertices()) {
    // Compute the symmetric point
    Vector3 pos = geom->position(v);

    // Test if canonical
    Vector3 vPlane = (pos - rotPoint) - dot(rotAxis, pos - rotPoint) * rotAxis;
    double canonicalAngle = angleInPlane(orthAxis, vPlane, rotAxis);
    bool isCanonical = norm(vPlane) < toleranceRadius || (canonicalAngle >= 0 && canonicalAngle < deltaTheta);
    if (isCanonical) {
      result.canonicalVertices.push_back(v);
    }

    for (int iRot = 1; iRot < nSym; iRot++) {
      double theta = iRot * deltaTheta;
      Vector3 rotPos = (pos - rotPoint).rotate_around(rotAxis, theta) + rotPoint;

      // If this point is its own pair, there's no mirror to look for (assumes
      // no duplicate verts)
      if (norm(pos - rotPos) < toleranceRadius) continue;

      // Search for the point
      size_t rotInd;
      bool success = findPoint(tree, rotPos, toleranceRadius, rotInd);
      if (!success) {
        return result;
      }

      // If found, add to lists
      if (isCanonical) {
        result.symmetrySet[v].push_back(mesh->vertex(rotInd));
      }
    }
  }

  result.symmetryFound = true;
  return result;
}

SymmetryResult detectSymmetryDoubleMirror(Geometry<Euclidean>* geom, KdTree* tree) {
  HalfedgeMesh* mesh = geom->getMesh();
  double toleranceRadius = geom->lengthScale() * 1e-5;

  SymmetryResult result;
  result.symmetryFound = false;
  result.symmetrySet = VertexData<std::vector<Vertex>>(*mesh, std::vector<Vertex>());

  for (Vertex v : mesh->vertices()) {
    // Compute the symmetric point
    Vector3 pos = geom->position(v);

    // Test if canonical
    bool isCanonical = pos.y >= 0 && pos.z >= 0;
    if (isCanonical) {
      result.canonicalVertices.push_back(v);
    }

    for (int iS = 1; iS < 4; iS++) {
      // Compute positions flipped across axes
      Vector3 mirrorPos = pos;
      if (iS % 2 == 1) {
        mirrorPos.y *= -1;
      }
      if (iS >= 2) {
        mirrorPos.z *= -1;
      }

      // If this point is its own pair, there's no mirror to look for (assumes
      // no duplicate verts)
      if (norm(pos - mirrorPos) < toleranceRadius) continue;

      // Search for the point
      size_t symInd;
      bool success = findPoint(tree, mirrorPos, toleranceRadius, symInd);
      if (!success) {
        return result;
      }

      // If found, add to lists
      if (isCanonical) {
        result.symmetrySet[v].push_back(mesh->vertex(symInd));
      }
    }
  }

  result.symmetryFound = true;
  return result;
}

} // namespace

SymmetryResult detectSymmetryMirror(Geometry<Euclidean>* geom, Vector3 planeNormal, Vector3 planePoint) {
  KdTree* tree = buildKDTree(geom);
  SymmetryResult r = detectSymmetryMirror(geom, planeNormal, planePoint, tree);
  delete tree;
  return r;
}

SymmetryResult detectSymmetryRotation(Geometry<Euclidean>* geom, Vector3 rotAxis, Vector3 rotPoint, int nSym) {
  KdTree* tree = buildKDTree(geom);
  SymmetryResult r = detectSymmetryRotation(geom, rotAxis, rotPoint, nSym, tree);
  delete tree;
  return r;
}

SymmetryResult detectSymmetryAuto(Geometry<Euclidean>* geom) {
  std::cout << "Attempting to automatically detect symmetry..." << std::endl;

  KdTree* tree = buildKDTree(geom);

  Vector3 center = geom->center();

  // == Mirror symmetry across coordinate axes, about center
  {
    SymmetryResult res = detectSymmetryMirror(geom, Vector3{1.0, 0.0, 0.0}, center, tree);
    if (res.symmetryFound) {
      cout << "  ... symmetry found across x-axis!" << endl;
      delete tree;
      return res;
    }
  }
  {
    SymmetryResult res = detectSymmetryMirror(geom, Vector3{0.0, 1.0, 0.0}, center, tree);
    if (res.symmetryFound) {
      cout << "  ... symmetry found across y-axis!" << endl;
      delete tree;
      return res;
    }
  }
  {
    SymmetryResult res = detectSymmetryMirror(geom, Vector3{0.0, 0.0, 1.0}, center, tree);
    if (res.symmetryFound) {
      cout << "  ... symmetry found across z-axis!" << endl;
      delete tree;
      return res;
    }
  }

  //  == Rotational symmetry about coordinate axes at center
  // (higher order symmetries are cooler)
  for (int nSym = 8; nSym >= 2; nSym--) {
    {
      SymmetryResult res = detectSymmetryRotation(geom, Vector3{1.0, 0.0, 0.0}, center, nSym, tree);
      if (res.symmetryFound) {
        cout << "  ... rotational symmetry found about x-axis with index " << nSym << "!" << endl;
        delete tree;
        return res;
      }
    }

    {
      SymmetryResult res = detectSymmetryRotation(geom, Vector3{0.0, 1.0, 0.0}, center, nSym, tree);
      if (res.symmetryFound) {
        cout << "  ... rotational symmetry found about y-axis with index " << nSym << "!" << endl;
        delete tree;
        return res;
      }
    }

    {
      SymmetryResult res = detectSymmetryRotation(geom, Vector3{0.0, 0.0, 1.0}, center, nSym, tree);
      if (res.symmetryFound) {
        cout << "  ... rotational symmetry found about z-axis with index " << nSym << "!" << endl;
        delete tree;
        return res;
      }
    }
  }

  cout << "  ...no symmetry found." << endl;
  delete tree;
  SymmetryResult r;
  r.symmetryFound = false;
  return r;
}

SymmetryResult detectSymmetryAutoMirror(Geometry<Euclidean>* geom) {
  cout << "Attempting to automatically detect mirror symmetry..." << endl;

  KdTree* tree = buildKDTree(geom);

  Vector3 center = geom->center();

  // == Mirror symmetry across coordinate axes, about center
  {
    SymmetryResult res = detectSymmetryMirror(geom, Vector3{1.0, 0.0, 0.0}, center, tree);
    if (res.symmetryFound) {
      cout << "  ... symmetry found across x-axis!" << endl;
      delete tree;
      return res;
    }
  }
  {
    SymmetryResult res = detectSymmetryMirror(geom, Vector3{0.0, 1.0, 0.0}, center, tree);
    if (res.symmetryFound) {
      cout << "  ... symmetry found across y-axis!" << endl;
      delete tree;
      return res;
    }
  }
  {
    SymmetryResult res = detectSymmetryMirror(geom, Vector3{0.0, 0.0, 1.0}, center, tree);
    if (res.symmetryFound) {
      cout << "  ... symmetry found across z-axis!" << endl;
      delete tree;
      return res;
    }
  }

  cout << "  ...no symmetry found." << endl;
  delete tree;
  SymmetryResult r;
  r.symmetryFound = false;
  return r;
}

SymmetryResult detectSymmetryAutoRotation(Geometry<Euclidean>* geom) {
  cout << "Attempting to automatically detect rotational symmetry..." << endl;

  KdTree* tree = buildKDTree(geom);

  Vector3 center = geom->center();

  //  == Rotational symmetry about coordinate axes at center
  // (higher order symmetries are cooler)
  for (int nSym = 8; nSym >= 2; nSym--) {
    {
      SymmetryResult res = detectSymmetryRotation(geom, Vector3{1.0, 0.0, 0.0}, center, nSym, tree);
      if (res.symmetryFound) {
        cout << "  ... rotational symmetry found about x-axis with index " << nSym << "!" << endl;
        delete tree;
        return res;
      }
    }

    {
      SymmetryResult res = detectSymmetryRotation(geom, Vector3{0.0, 1.0, 0.0}, center, nSym, tree);
      if (res.symmetryFound) {
        cout << "  ... rotational symmetry found about y-axis with index " << nSym << "!" << endl;
        delete tree;
        return res;
      }
    }

    {
      SymmetryResult res = detectSymmetryRotation(geom, Vector3{0.0, 0.0, 1.0}, center, nSym, tree);
      if (res.symmetryFound) {
        cout << "  ... rotational symmetry found about z-axis with index " << nSym << "!" << endl;
        delete tree;
        return res;
      }
    }
  }

  cout << "  ...no symmetry found." << endl;
  delete tree;
  SymmetryResult r;
  r.symmetryFound = false;
  return r;
}

SymmetryResult detectSymmetryDoubleMirror(Geometry<Euclidean>* geom) {
  KdTree* tree = buildKDTree(geom);
  SymmetryResult r = detectSymmetryDoubleMirror(geom, tree);
  delete tree;
  return r;
}

} // namespace surface
} // namespace geometrycentral
