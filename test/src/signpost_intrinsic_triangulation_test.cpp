#include "geometrycentral/surface/signpost_intrinsic_triangulation.h"

#include "load_test_meshes.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

// helpers
namespace {

std::mt19937 mt(42);

double unitRand() {
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  return dist(mt);
}

void insertVerticesAtRandom(SignpostIntrinsicTriangulation& signpostTri) {
  for (Face f : signpostTri.intrinsicMesh->faces()) {
    if (unitRand() > 0.5) {
      signpostTri.insertVertex(SurfacePoint{f, {1/3., 1/3., 1/3.}});
    }
  }
}

void flipEdgesAtRandom(SignpostIntrinsicTriangulation& signpostTri) {
  for (Edge e : signpostTri.intrinsicMesh->edges()) {
    if (unitRand() > 0.5) {
      signpostTri.flipEdgeIfPossible(e);
    }
  }
}

void mutateIntrinsicMesh(SignpostIntrinsicTriangulation& signpostTri) {
  insertVerticesAtRandom(signpostTri);
  flipEdgesAtRandom(signpostTri);
}

double getBoundingBoxDiagonal(VertexPositionGeometry& geometry) {
  Vector3 minP = Vector3::constant(std::numeric_limits<double>::infinity());
  Vector3 maxP = Vector3::constant(-std::numeric_limits<double>::infinity());
  for (Vertex v : geometry.mesh.vertices()) {
    Vector3 p = geometry.vertexPositions[v];
    minP = componentwiseMin(minP, p);
    maxP = componentwiseMax(maxP, p);
  }
  return norm(minP - maxP);
}

const int N = 7;
const double EPS = 1e-3;

} // namespace

class SignpostIntrinsicTriangulationSuite : public MeshAssetSuite {};

// ============================================================
// =============== Equivalent point tests
// ============================================================

TEST_F(SignpostIntrinsicTriangulationSuite, EquivalentPointOnIntrinsic_FacePoint) {
  for (auto& asset : {getAsset("bob_small.ply", true), getAsset("lego.ply", true), getAsset("sphere_small.ply", true), getAsset("spot.ply", true)}) {
    asset.printThyName();
    SignpostIntrinsicTriangulation signpostTri(*asset.manifoldMesh, *asset.geometry);

    mutateIntrinsicMesh(signpostTri);

    double boundingBoxDiagonal = getBoundingBoxDiagonal(*asset.geometry);

    for (Face f : signpostTri.inputMesh.faces()) {
      Vector3 i;
      for (i[0] = 1; i[0] < N; ++i[0]) {
        for (i[1] = 1; i[1] < N - i[0]; ++i[1]) {
          i[2] = N - i[0] - i[1];
          Vector3 faceCoords = i / N;
          SurfacePoint pointOnInput_before{f, faceCoords};
          SurfacePoint pointOnIntrinsic = signpostTri.equivalentPointOnIntrinsic(pointOnInput_before);
          SurfacePoint pointOnInput_after = signpostTri.equivalentPointOnInput(pointOnIntrinsic);
          double error = (pointOnInput_before.interpolate(asset.geometry->inputVertexPositions) - pointOnInput_after.interpolate(asset.geometry->inputVertexPositions)).norm();
          EXPECT_LT(error, boundingBoxDiagonal * EPS);
        }
      }
    }
  }
}

TEST_F(SignpostIntrinsicTriangulationSuite, EquivalentPointOnInput_FacePoint) {
  for (auto& asset : {getAsset("bob_small.ply", true), getAsset("lego.ply", true), getAsset("sphere_small.ply", true), getAsset("spot.ply", true)}) {
    asset.printThyName();
    SignpostIntrinsicTriangulation signpostTri(*asset.manifoldMesh, *asset.geometry);

    mutateIntrinsicMesh(signpostTri);

    double boundingBoxDiagonal = getBoundingBoxDiagonal(*asset.geometry);

    for (Face f : signpostTri.intrinsicMesh->faces()) {
      Vector3 i;
      for (i[0] = 1; i[0] < N; ++i[0]) {
        for (i[1] = 1; i[1] < N - i[0]; ++i[1]) {
          i[2] = N - i[0] - i[1];
          Vector3 faceCoords = i / N;
          SurfacePoint pointOnIntrinsic_before{f, faceCoords};
          SurfacePoint pointOnInput_before = signpostTri.equivalentPointOnInput(pointOnIntrinsic_before);
          SurfacePoint pointOnIntrinsic_after = signpostTri.equivalentPointOnIntrinsic(pointOnInput_before);
          SurfacePoint pointOnInput_after = signpostTri.equivalentPointOnInput(pointOnIntrinsic_after);
          double error = (pointOnInput_before.interpolate(asset.geometry->inputVertexPositions) - pointOnInput_after.interpolate(asset.geometry->inputVertexPositions)).norm();
          EXPECT_LT(error, boundingBoxDiagonal * EPS);
        }
      }
    }
  }
}

TEST_F(SignpostIntrinsicTriangulationSuite, EquivalentPointOnIntrinsic_EdgePoint) {
  for (auto& asset : {getAsset("bob_small.ply", true), getAsset("lego.ply", true), getAsset("sphere_small.ply", true), getAsset("spot.ply", true)}) {
    asset.printThyName();
    SignpostIntrinsicTriangulation signpostTri(*asset.manifoldMesh, *asset.geometry);

    mutateIntrinsicMesh(signpostTri);

    double boundingBoxDiagonal = getBoundingBoxDiagonal(*asset.geometry);

    for (Edge e : signpostTri.inputMesh.edges()) {
      for (double i = 1; i < N; ++i) {
        SurfacePoint pointOnInput_before{e, i / N};
        SurfacePoint pointOnIntrinsic = signpostTri.equivalentPointOnIntrinsic(pointOnInput_before);
        SurfacePoint pointOnInput_after = signpostTri.equivalentPointOnInput(pointOnIntrinsic);
        double error = (pointOnInput_before.interpolate(asset.geometry->inputVertexPositions) - pointOnInput_after.interpolate(asset.geometry->inputVertexPositions)).norm();
        EXPECT_LT(error, boundingBoxDiagonal * EPS);
      }
    }
  }
}

TEST_F(SignpostIntrinsicTriangulationSuite, EquivalentPointOnInput_EdgePoint) {
  for (auto& asset : {getAsset("bob_small.ply", true), getAsset("lego.ply", true), getAsset("sphere_small.ply", true), getAsset("spot.ply", true)}) {
    asset.printThyName();
    SignpostIntrinsicTriangulation signpostTri(*asset.manifoldMesh, *asset.geometry);

    mutateIntrinsicMesh(signpostTri);

    double boundingBoxDiagonal = getBoundingBoxDiagonal(*asset.geometry);

    for (Edge e : signpostTri.intrinsicMesh->edges()) {
      for (double i = 1; i < N; ++i) {
        SurfacePoint pointOnIntrinsic_before{e, i / N};
        SurfacePoint pointOnInput_before = signpostTri.equivalentPointOnInput(pointOnIntrinsic_before);
        SurfacePoint pointOnIntrinsic_after = signpostTri.equivalentPointOnIntrinsic(pointOnInput_before);
        SurfacePoint pointOnInput_after = signpostTri.equivalentPointOnInput(pointOnIntrinsic_after);
        double error = (pointOnInput_before.interpolate(asset.geometry->inputVertexPositions) - pointOnInput_after.interpolate(asset.geometry->inputVertexPositions)).norm();
        EXPECT_LT(error, boundingBoxDiagonal * EPS);
      }
    }
  }
}

TEST_F(SignpostIntrinsicTriangulationSuite, EquivalentPointOnIntrinsic_VertexPoint) {
  for (auto& asset : {getAsset("bob_small.ply", true), getAsset("lego.ply", true), getAsset("sphere_small.ply", true), getAsset("spot.ply", true)}) {
    asset.printThyName();
    SignpostIntrinsicTriangulation signpostTri(*asset.manifoldMesh, *asset.geometry);

    mutateIntrinsicMesh(signpostTri);

    double boundingBoxDiagonal = getBoundingBoxDiagonal(*asset.geometry);

    for (Vertex v : signpostTri.inputMesh.vertices()) {
      SurfacePoint pointOnInput_before{v};
      SurfacePoint pointOnIntrinsic = signpostTri.equivalentPointOnIntrinsic(pointOnInput_before);
      SurfacePoint pointOnInput_after = signpostTri.equivalentPointOnInput(pointOnIntrinsic);
      double error = (pointOnInput_before.interpolate(asset.geometry->inputVertexPositions) - pointOnInput_after.interpolate(asset.geometry->inputVertexPositions)).norm();
      EXPECT_LT(error, boundingBoxDiagonal * EPS);
    }
  }
}

TEST_F(SignpostIntrinsicTriangulationSuite, EquivalentPointOnInput_VertexPoint) {
  for (auto& asset : {getAsset("bob_small.ply", true), getAsset("lego.ply", true), getAsset("sphere_small.ply", true), getAsset("spot.ply", true)}) {
    asset.printThyName();
    SignpostIntrinsicTriangulation signpostTri(*asset.manifoldMesh, *asset.geometry);

    mutateIntrinsicMesh(signpostTri);

    double boundingBoxDiagonal = getBoundingBoxDiagonal(*asset.geometry);

    for (Vertex v : signpostTri.intrinsicMesh->vertices()) {
      SurfacePoint pointOnIntrinsic_before{v};
      SurfacePoint pointOnInput_before = signpostTri.equivalentPointOnInput(pointOnIntrinsic_before);
      SurfacePoint pointOnIntrinsic_after = signpostTri.equivalentPointOnIntrinsic(pointOnInput_before);
      SurfacePoint pointOnInput_after = signpostTri.equivalentPointOnInput(pointOnIntrinsic_after);
      double error = (pointOnInput_before.interpolate(asset.geometry->inputVertexPositions) - pointOnInput_after.interpolate(asset.geometry->inputVertexPositions)).norm();
      EXPECT_LT(error, boundingBoxDiagonal * EPS);
    }
  }
}
