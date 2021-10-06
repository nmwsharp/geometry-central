
#include "geometrycentral/surface/manifold_surface_mesh.h"
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


// Split a few edges on a bunch of meshes
TEST_F(HalfedgeMutationSuite, InsertVertexAlongEdgeTest) {

  for (MeshAsset& a : manifoldSurfaceMeshes()) {
    a.printThyName();

    int count = 10;
    int indInc = static_cast<int>(std::ceil(a.mesh->nVertices() / static_cast<double>(count)));

    int ind = 0;
    for (int i = 0; i < count; i++) {

      // Insert along an edge
      Edge e = a.manifoldMesh->edge(ind);
      a.manifoldMesh->insertVertexAlongEdge(e);
      a.manifoldMesh->validateConnectivity();

      ind = (ind + indInc) % a.manifoldMesh->nVertices();
    }
  }
}


// Insert a vertex along every edge and triangulate (not-quite subdivision)
TEST_F(HalfedgeMutationSuite, InsertVertexAndTriangulateSubdivideTest) {

  for (MeshAsset& a : manifoldSurfaceMeshes()) {
    a.printThyName();

    // Split every edge
    std::vector<Edge> origEdges;
    for (Edge e : a.manifoldMesh->edges()) {
      origEdges.push_back(e);
    }
    for (Edge e : origEdges) {
      a.manifoldMesh->insertVertexAlongEdge(e);
    }

    a.manifoldMesh->validateConnectivity();

    // Triangulate
    // TODO this loops while modifying. Do we allow that?
    for (Face f : a.manifoldMesh->faces()) {
      a.manifoldMesh->triangulate(f);
    }

    a.manifoldMesh->validateConnectivity();
  }
}

// Split every edge and then flip (regular subdivision)
TEST_F(HalfedgeMutationSuite, SplitFlipSubdivide) {

  for (MeshAsset& a : triangularMeshes()) {
    a.printThyName();

    VertexData<char> isNewVertex(*a.manifoldMesh, false);
    for (Vertex v : a.manifoldMesh->vertices()) {
      isNewVertex[v] = true;
    }

    // Split every edge
    std::vector<Edge> origEdges;
    for (Edge e : a.manifoldMesh->edges()) {
      origEdges.push_back(e);
    }
    for (Edge e : origEdges) {
      a.manifoldMesh->splitEdgeTriangular(e);
    }
    a.manifoldMesh->validateConnectivity();

    // Flip edges between old and new
    for (Edge e : a.manifoldMesh->edges()) {
      if (isNewVertex[e.halfedge().vertex()] != isNewVertex[e.halfedge().twin().vertex()]) {
        a.manifoldMesh->flip(e);
      }
    }
    a.manifoldMesh->validateConnectivity();

    // Should yield subdivision
    for (Face f : a.manifoldMesh->faces()) {
      EXPECT_TRUE(f.isTriangle());
    }
  }
}

// Split a few edges on a bunch of meshes
TEST_F(HalfedgeMutationSuite, EdgeSplitTest) {

  for (MeshAsset& a : triangularMeshes()) {
    a.printThyName();

    int count = 10;
    int indInc = static_cast<int>(std::ceil(a.manifoldMesh->nVertices() / static_cast<double>(count)));

    int splitInd = 0;
    for (int i = 0; i < count; i++) {

      // Split an edge
      Edge eSplit = a.manifoldMesh->edge(splitInd);
      a.manifoldMesh->splitEdgeTriangular(eSplit);
      a.manifoldMesh->validateConnectivity();

      splitInd = (splitInd + indInc) % a.manifoldMesh->nVertices();
    }
  }
}


// Invert face orientation on a bunch of meshes
TEST_F(HalfedgeMutationSuite, InvertOrientationTest) {

  for (MeshAsset& a : allMeshes()) {
    if (a.isSubclassManifoldSurfaceMesh) continue;
    a.printThyName();

    int count = 10;
    int ind = 0;
    for (int i = 0; i < count; i++) {

      // Invert
      Face f = a.mesh->face(ind);
      a.mesh->invertOrientation(f);
      a.mesh->validateConnectivity();

      ind = (ind + 1) % a.mesh->nFaces();
    }
  }
}


TEST_F(HalfedgeMutationSuite, DuplicateFaceTest) {

  for (const MeshAsset& a : {getAsset("triple_vierbein.obj", false)}) {
    a.printThyName();

    int count = 10;
    int ind = 0;
    for (int i = 0; i < count; i++) {

      // Invert
      Face f = a.mesh->face(ind);
      Face newF = a.mesh->duplicateFace(f);
      a.mesh->validateConnectivity();

      ind = (ind + 1) % a.mesh->nFaces();
    }
  }
}

// == A few higher level tests which do many operations

TEST_F(HalfedgeMutationSuite, SeparateEdgesTest) {

  for (const MeshAsset& a : {getAsset("triple_vierbein.obj", false)}) {
    a.printThyName();
    a.mesh->separateNonmanifoldEdges();
    a.mesh->validateConnectivity();

    EXPECT_TRUE(a.mesh->isEdgeManifold());
  }
}

TEST_F(HalfedgeMutationSuite, SeparateEdgesAndVerticesTest) {

  for (const MeshAsset& a : {getAsset("triple_vierbein.obj", false)}) {
    a.printThyName();
    a.mesh->separateNonmanifoldEdges();
    a.mesh->separateNonmanifoldVertices();
    a.mesh->validateConnectivity();

    EXPECT_TRUE(a.mesh->isManifold());
  }
}


TEST_F(HalfedgeMutationSuite, GreedyOrientTest) {

  for (const MeshAsset& a : {getAsset("bob_small.ply", false)}) {
    a.printThyName();

    // Do a bunch of random inversion
    int count = 100;
    int ind = 0;
    for (int i = 0; i < count; i++) {

      // Invert
      Face f = a.mesh->face(ind);
      a.mesh->invertOrientation(f);

      ind = (ind + 13) % a.mesh->nFaces();
    }


    a.mesh->greedilyOrientFaces();
    a.mesh->validateConnectivity();

    EXPECT_TRUE(a.mesh->isOriented());
  }
}

TEST_F(HalfedgeMutationSuite, ToManifoldTest) {

  for (const MeshAsset& a : {
           getAsset("bob_small.ply", false),
           getAsset("hourglass_ico.obj", false),
           getAsset("lego.ply", false),
       }) {
    a.printThyName();

    a.mesh->separateNonmanifoldEdges();
    EXPECT_TRUE(a.mesh->isEdgeManifold());

    a.mesh->separateNonmanifoldVertices();
    EXPECT_TRUE(a.mesh->isManifold());

    a.mesh->greedilyOrientFaces();
    EXPECT_TRUE(a.mesh->isOriented());

    a.mesh->validateConnectivity();

    std::unique_ptr<ManifoldSurfaceMesh> manifMesh = a.mesh->toManifoldMesh();
    manifMesh->validateConnectivity();

    EXPECT_EQ(a.mesh->nVertices(), manifMesh->nVertices());
    EXPECT_EQ(a.mesh->nFaces(), manifMesh->nFaces());
    EXPECT_EQ(a.mesh->nEdges(), manifMesh->nEdges());
  }
}


// Flip a lot of edges on one mesh without boundary
TEST_F(HalfedgeMutationSuite, EdgeFlipClosedManyTest) {

  for (const MeshAsset& asset : {getAsset("sphere_small.ply", true)}) {
    SurfaceMesh& mesh = *asset.mesh;

    int count = 1000;
    int indInc = static_cast<int>(std::ceil(mesh.nVertices() / static_cast<double>(count)));

    int flipInd = 0;
    for (int i = 0; i < count; i++) {

      // Flip an edge
      Edge eFlip = mesh.edge(flipInd);
      bool didFlip = mesh.flip(eFlip);
      // mesh.validateConnectivity();

      flipInd = (flipInd + 1) % mesh.nVertices();
    }

    mesh.validateConnectivity();
  }
}

// Flip a lot of edges and orientations on one mesh
TEST_F(HalfedgeMutationSuite, EdgeFlipInvertOrientClosedManyTest) {

  for (const MeshAsset& asset : {getAsset("sphere_small.ply", false)}) {
    SurfaceMesh& mesh = *asset.mesh;

    int count = 1000;
    int flipInd = 0;
    int invertInd = mesh.nEdges() / 2;
    for (int i = 0; i < count; i++) {

      // Invert a faces
      Face fInvert = mesh.face(invertInd);
      mesh.invertOrientation(fInvert);
      // mesh.validateConnectivity();

      // Flip an edge
      Edge eFlip = mesh.edge(flipInd);
      bool didFlip = mesh.flip(eFlip);
      // mesh.validateConnectivity();

      flipInd = (flipInd + 1) % mesh.nEdges();
      invertInd = (invertInd + 1) % mesh.nFaces();
    }

    mesh.validateConnectivity();
  }
}


// =====================================================
// ========= Removal tests
// =====================================================

TEST_F(HalfedgeMutationSuite, RemoveVertex) {

  for (const MeshAsset& a : {getAsset("bob_small.ply", true)}) {
    a.printThyName();

    // Remove some vertices
    a.manifoldMesh->validateConnectivity();
    a.manifoldMesh->removeVertex(a.manifoldMesh->vertex(7));
    a.manifoldMesh->validateConnectivity();
    a.manifoldMesh->removeVertex(a.manifoldMesh->vertex(12));
    a.manifoldMesh->validateConnectivity();
    a.manifoldMesh->removeVertex(a.manifoldMesh->vertex(44));
    a.manifoldMesh->validateConnectivity();
  }
}

TEST_F(HalfedgeMutationSuite, RemoveVertexAndCompress) {

  for (const MeshAsset& a : {getAsset("bob_small.ply", true)}) {
    a.printThyName();

    a.manifoldMesh->validateConnectivity();

    // Remove a vertex
    a.manifoldMesh->removeVertex(a.manifoldMesh->vertex(7));
    a.manifoldMesh->validateConnectivity();

    a.manifoldMesh->compress();
    a.manifoldMesh->validateConnectivity();
  }
}


TEST_F(HalfedgeMutationSuite, CollapseEdge) {

  for (const MeshAsset& a : {getAsset("bob_small.ply", true)}) {
    a.printThyName();

    a.manifoldMesh->validateConnectivity();

    // Collapse an edge
    a.manifoldMesh->collapseEdgeTriangular(a.manifoldMesh->edge(7));
    a.manifoldMesh->validateConnectivity();

    a.manifoldMesh->compress();
    a.manifoldMesh->validateConnectivity();
  }
}


TEST_F(HalfedgeMutationSuite, CollapseEdgeBoundary) {

  for (const MeshAsset& a : {getAsset("lego.ply", true)}) {
    a.printThyName();

    a.manifoldMesh->validateConnectivity();

    Edge boundaryEdge = a.manifoldMesh->edge(163);
    ASSERT_TRUE(boundaryEdge.isBoundary());

    // Collapse an edge
    a.manifoldMesh->collapseEdgeTriangular(boundaryEdge);
    a.manifoldMesh->validateConnectivity();

    a.manifoldMesh->compress();
    a.manifoldMesh->validateConnectivity();
  }
}

// =====================================================
// ========= Container tests
// =====================================================

TEST_F(HalfedgeMutationSuite, ContainerExpandTest) {

  auto asset = getAsset("lego.ply", true);
  ManifoldSurfaceMesh& mesh = *asset.manifoldMesh;
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

TEST_F(HalfedgeMutationSuite, ContainerCompress) {

  for (const MeshAsset& a : {getAsset("bob_small.ply", true)}) {
    a.printThyName();

    // Create a container
    VertexData<int> values(*a.manifoldMesh, 7);

    // Remove a vertex
    a.manifoldMesh->removeVertex(a.manifoldMesh->vertex(12));

    // Iterate through and check values
    EXPECT_GT(values.size(), a.manifoldMesh->nVertices()); // buffer should be larger after an expansion
    for (Vertex v : a.manifoldMesh->vertices()) {
      EXPECT_EQ(values[v], 7);
    }

    // Compress
    a.manifoldMesh->compress();

    // Iterate through and check values
    EXPECT_EQ(values.size(), a.manifoldMesh->nVertices());
    for (Vertex v : a.manifoldMesh->vertices()) {
      EXPECT_EQ(values[v], 7);
    }
  }
}

// Test that edge containers get updated properly (it's a bit of a special case in the implicit twin implementation)
TEST_F(HalfedgeMutationSuite, ContainerCompressEdge) {

  for (const MeshAsset& a : {getAsset("bob_small.ply", true)}) {
    a.printThyName();

    // Create a container
    EdgeData<int> values(*a.manifoldMesh, 7);

    // Remove a vertex
    a.manifoldMesh->removeVertex(a.manifoldMesh->vertex(12));

    // Iterate through and check values
    EXPECT_GT(values.size(), a.manifoldMesh->nEdges());
    for (Edge e : a.manifoldMesh->edges()) {
      EXPECT_EQ(values[e], 7);
    }

    // Compress
    a.manifoldMesh->compress();

    // Iterate through and check values
    EXPECT_EQ(values.size(), a.manifoldMesh->nEdges());
    for (Edge e : a.manifoldMesh->edges()) {
      EXPECT_EQ(values[e], 7);
    }
  }
}

// do ALL the containers
TEST_F(HalfedgeMutationSuite, ContainerCompressAll) {

  for (const MeshAsset& a : {getAsset("bob_small.ply", true), getAsset("lego.ply", true)}) {
    a.printThyName();

    // Create a container
    VertexData<int> values_vertex(*a.manifoldMesh, 7);
    HalfedgeData<int> values_halfedge(*a.manifoldMesh, 7);
    EdgeData<int> values_edge(*a.manifoldMesh, 7);
    FaceData<int> values_face(*a.manifoldMesh, 7);
    BoundaryLoopData<int> values_bl(*a.manifoldMesh, 7);

    auto checkVertex = [&]() {
      //EXPECT_EQ(values_vertex.size(), a.manifoldMesh->nVertices());
      for (Vertex v : a.manifoldMesh->vertices()) {
        EXPECT_EQ(values_vertex[v], 7);
      }
    };
    auto checkHalfedge = [&]() {
      //EXPECT_EQ(values_halfedge.size(), a.manifoldMesh->nHalfedges());
      for (Halfedge he : a.manifoldMesh->halfedges()) {
        EXPECT_EQ(values_halfedge[he], 7);
      }
    };
    auto checkEdge = [&]() {
      //EXPECT_EQ(values_edge.size(), a.manifoldMesh->nEdges());
      for (Edge e : a.manifoldMesh->edges()) {
        EXPECT_EQ(values_edge[e], 7);
      }
    };
    auto checkFace = [&]() {
      //EXPECT_EQ(values_face.size(), a.manifoldMesh->nFaces());
      for (Face f : a.manifoldMesh->faces()) {
        EXPECT_EQ(values_face[f], 7);
      }
    };
    auto checkBoundaryLoop = [&]() {
      //EXPECT_EQ(values_bl.size(), a.manifoldMesh->nBoundaryLoops());
      for (BoundaryLoop bl : a.manifoldMesh->boundaryLoops()) {
        EXPECT_EQ(values_bl[bl], 7);
      }
    };

    // Remove a vertex
    a.manifoldMesh->removeVertex(a.manifoldMesh->vertex(144));

    checkVertex();
    checkHalfedge();
    checkEdge();
    checkFace();
    checkBoundaryLoop();


    // Compress
    a.manifoldMesh->compress();

    checkVertex();
    checkHalfedge();
    checkEdge();
    checkFace();
    checkBoundaryLoop();
  }
}

// =====================================================
// ========= Mutation helper tests
// =====================================================
