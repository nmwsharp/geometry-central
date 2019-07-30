
#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/meshio.h"

#include "geometrycentral/surface/base_geometry_interface.h"
#include "geometrycentral/surface/edge_length_geometry.h"
#include "geometrycentral/surface/embedded_geometry_interface.h"
#include "geometrycentral/surface/extrinsic_geometry_interface.h"
#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "load_test_meshes.h"

#include "gtest/gtest.h"

#include <iostream>
#include <string>
#include <unordered_set>


using namespace geometrycentral;
using namespace geometrycentral::surface;
using std::cout;
using std::endl;

class HalfedgeMutationSuite : public MeshAssetSuite {};

// Flip a few edges on a bunch of meshes
TEST_F(HalfedgeMutationSuite, EdgeFlipTest) {

  for (MeshAsset& a : allMeshes()) {
    a.printThyName();

    int count = 10;
    int indInc = static_cast<int>(std::ceil(a.mesh->nVertices() / static_cast<double>(count)));

    int flipInd = 0;
    for (int i = 0; i < count; i++) {

      // Flip an edge
      Edge eFlip = a.mesh->edge(flipInd);
      a.mesh->flip(eFlip);
      a.mesh->validateConnectivity();

      flipInd = (flipInd + indInc) % a.mesh->nVertices();
    }
  }
}

// Flip a lot of edges on one mesh without boundary
TEST_F(HalfedgeMutationSuite, EdgeFlipClosedManyTest) {

  auto asset = getAsset("sphere_small.ply");
  HalfedgeMesh& mesh = *asset.mesh;

  int count = 1000;
  int indInc = static_cast<int>(std::ceil(mesh.nVertices() / static_cast<double>(count)));

  int flipInd = 0;
  for (int i = 0; i < count; i++) {

    // Flip an edge
    Edge eFlip = mesh.edge(flipInd);
    mesh.flip(eFlip);
    mesh.validateConnectivity();

    flipInd = (flipInd + 1) % mesh.nVertices();
  }
}

// Split a few edges on a bunch of meshes
TEST_F(HalfedgeMutationSuite, InsertVertexAlongEdgeTest) {

  for (MeshAsset& a : allMeshes()) {
    a.printThyName();

    int count = 10;
    int indInc = static_cast<int>(std::ceil(a.mesh->nVertices() / static_cast<double>(count)));

    int ind = 0;
    for (int i = 0; i < count; i++) {

      // Insert along an edge
      Edge e = a.mesh->edge(ind);
      a.mesh->insertVertexAlongEdge(e);
      a.mesh->validateConnectivity();

      ind = (ind + indInc) % a.mesh->nVertices();
    }
  }
}


// Insert a vertex along every edge and triangulate (not-quite subdivision)
TEST_F(HalfedgeMutationSuite, InsertVertexAndTriangulateSubdivideTest) {

  for (MeshAsset& a : allMeshes()) {
    a.printThyName();

    // Split every edge
    std::vector<Edge> origEdges;
    for (Edge e : a.mesh->edges()) {
      origEdges.push_back(e);
    }
    for (Edge e : origEdges) {
      a.mesh->insertVertexAlongEdge(e);
    }

    a.mesh->validateConnectivity();

    // Triangulate
    // TODO this loops while modifying. Do we allow that?
    for (Face f : a.mesh->faces()) {
      a.mesh->triangulate(f);
    }

    a.mesh->validateConnectivity();
  }
}

// Split every edge and then flip (regular subdivision)
TEST_F(HalfedgeMutationSuite, SplitFlipSubdivide) {

  for (MeshAsset& a : triangularMeshes()) {
    a.printThyName();

    VertexData<char> isNewVertex(*a.mesh, false);
    for (Vertex v : a.mesh->vertices()) {
      isNewVertex[v] = true;
    }

    // Split every edge
    std::vector<Edge> origEdges;
    for (Edge e : a.mesh->edges()) {
      origEdges.push_back(e);
    }
    for (Edge e : origEdges) {
      a.mesh->splitEdgeTriangular(e);
    }
    a.mesh->validateConnectivity();

    // Flip edges between old and new
    for (Edge e : a.mesh->edges()) {
      if (isNewVertex[e.halfedge().vertex()] != isNewVertex[e.halfedge().twin().vertex()]) {
        a.mesh->flip(e);
      }
    }
    a.mesh->validateConnectivity();

    // Should yield subdivision
    for (Face f : a.mesh->faces()) {
      EXPECT_TRUE(f.isTriangle());
    }
  }
}

// Split a few edges on a bunch of meshes
TEST_F(HalfedgeMutationSuite, EdgeSplitTest) {

  for (MeshAsset& a : triangularMeshes()) {
    a.printThyName();

    int count = 10;
    int indInc = static_cast<int>(std::ceil(a.mesh->nVertices() / static_cast<double>(count)));

    int splitInd = 0;
    for (int i = 0; i < count; i++) {

      // Split an edge
      Edge eSplit = a.mesh->edge(splitInd);
      a.mesh->splitEdgeTriangular(eSplit);
      a.mesh->validateConnectivity();

      splitInd = (splitInd + indInc) % a.mesh->nVertices();
    }
  }
}

// =====================================================
// ========= Container tests
// =====================================================

TEST_F(HalfedgeMutationSuite, ContainerExpandTest) {

  auto asset = getAsset("lego.ply");
  HalfedgeMesh& mesh = *asset.mesh;
  VertexPositionGeometry& origGeometry = *asset.geometry;

  // Initial element counts
  size_t nVertexOrig = mesh.nVertices();
  size_t nHalfedgeOrig = mesh.nHalfedges();
  size_t nCornerOrig = mesh.nCorners();
  size_t nEdgeOrig = mesh.nEdges();
  size_t nFaceOrig = mesh.nFaces();

  // Some containers. Set a default value too.
  VertexData<int> vData(mesh, 42);
  HalfedgeData<int> heData(mesh, 42);
  CornerData<int> cData(mesh, 42);
  EdgeData<int> eData(mesh, 42);
  FaceData<int> fData(mesh, 42);

  // Set a different value for all existing element
  // NOTE: does not test boundary loops
  for (Vertex e : mesh.vertices()) vData[e] = 17;
  for (Halfedge e : mesh.halfedges()) heData[e] = 17;
  for (Corner e : mesh.corners()) cData[e] = 17;
  for (Edge e : mesh.edges()) eData[e] = 17;
  for (Face e : mesh.faces()) fData[e] = 17;

  // Do just one opertation, to trigger a single resize
  // (this adds at least one of each element type
  mesh.splitEdgeTriangular(mesh.edge(0));

  // Be sure the mesh actually got bigger
  EXPECT_LT(nVertexOrig, mesh.nVertices());
  EXPECT_LT(nHalfedgeOrig, mesh.nHalfedges());
  EXPECT_LT(nCornerOrig, mesh.nCorners());
  EXPECT_LT(nEdgeOrig, mesh.nEdges());
  EXPECT_LT(nFaceOrig, mesh.nFaces());

  // == Index all containers to make sure they grew. Also, make sure new elements got the default value, not the value
  // we set for existing elements.

  { // vertices
    size_t origValCount = 0;
    for (Vertex e : mesh.vertices()) {
      EXPECT_TRUE(vData[e] == 17 || vData[e] == 42);
      if (vData[e] == 17) origValCount++;
    }
    EXPECT_EQ(origValCount, nVertexOrig);
  }

  { // halfedges
    size_t origValCount = 0;
    for (Halfedge e : mesh.halfedges()) {
      EXPECT_TRUE(heData[e] == 17 || heData[e] == 42);
      if (heData[e] == 17) origValCount++;
    }
    EXPECT_EQ(origValCount, nHalfedgeOrig);
  }

  { // corners
    size_t origValCount = 0;
    for (Corner e : mesh.corners()) {
      EXPECT_TRUE(cData[e] == 17 || cData[e] == 42);
      if (cData[e] == 17) origValCount++;
    }
    EXPECT_EQ(origValCount, nCornerOrig);
  }

  { // edges
    size_t origValCount = 0;
    for (Edge e : mesh.edges()) {
      EXPECT_TRUE(eData[e] == 17 || eData[e] == 42);
      if (eData[e] == 17) origValCount++;
    }
    EXPECT_EQ(origValCount, nEdgeOrig);
  }

  { // faces
    size_t origValCount = 0;
    for (Face e : mesh.faces()) {
      EXPECT_TRUE(fData[e] == 17 || fData[e] == 42);
      if (fData[e] == 17) origValCount++;
    }
    EXPECT_EQ(origValCount, nFaceOrig);
  }


  // Do a whole bunch of mesh operations, which should trigger several resizes
  for (int i = 0; i < 2; i++) {
    std::vector<Edge> origEdges;
    for (Edge e : mesh.edges()) {
      origEdges.push_back(e);
    }
    for (Edge e : origEdges) {
      mesh.splitEdgeTriangular(e);
    }
    mesh.validateConnectivity();
    for (Face f : mesh.faces()) {
      mesh.triangulate(f);
    }
  }

  // Check the same expansion invariants as above again
  { // vertices
    size_t origValCount = 0;
    for (Vertex e : mesh.vertices()) {
      EXPECT_TRUE(vData[e] == 17 || vData[e] == 42);
      if (vData[e] == 17) origValCount++;
    }
    EXPECT_EQ(origValCount, nVertexOrig);
  }

  { // halfedges
    size_t origValCount = 0;
    for (Halfedge e : mesh.halfedges()) {
      EXPECT_TRUE(heData[e] == 17 || heData[e] == 42);
      if (heData[e] == 17) origValCount++;
    }
    EXPECT_EQ(origValCount, nHalfedgeOrig);
  }

  { // corners
    size_t origValCount = 0;
    for (Corner e : mesh.corners()) {
      EXPECT_TRUE(cData[e] == 17 || cData[e] == 42);
      if (cData[e] == 17) origValCount++;
    }
    EXPECT_EQ(origValCount, nCornerOrig);
  }

  { // edges
    size_t origValCount = 0;
    for (Edge e : mesh.edges()) {
      EXPECT_TRUE(eData[e] == 17 || eData[e] == 42);
      if (eData[e] == 17) origValCount++;
    }
    EXPECT_EQ(origValCount, nEdgeOrig);
  }

  { // faces
    size_t origValCount = 0;
    for (Face e : mesh.faces()) {
      EXPECT_TRUE(fData[e] == 17 || fData[e] == 42);
      if (fData[e] == 17) origValCount++;
    }
    EXPECT_EQ(origValCount, nFaceOrig);
  }
}


// =====================================================
// ========= Mutation helper tests
// =====================================================
