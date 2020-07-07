#include "geometrycentral/surface/tufted_laplacian.h"
#include "geometrycentral/surface/edge_length_geometry.h"
#include "geometrycentral/surface/intrinsic_mollification.h"
#include "geometrycentral/surface/simple_idt.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

namespace geometrycentral {
namespace surface {

std::tuple<SparseMatrix<double>, SparseMatrix<double>>
buildTuftedLaplacian(SurfaceMesh& mesh, EmbeddedGeometryInterface& geom, double relativeMollificationFactor) {

  // Create a copy of the mesh / geometry to operate on
  std::unique_ptr<SurfaceMesh> tuftedMesh = mesh.copyToSurfaceMesh();
  geom.requireVertexPositions();
  VertexPositionGeometry tuftedGeom(*tuftedMesh, geom.vertexPositions.reinterpretTo(*tuftedMesh));
  tuftedGeom.requireEdgeLengths();
  EdgeData<double> tuftedEdgeLengths = tuftedGeom.edgeLengths;

  // Mollify, if requested
  if (relativeMollificationFactor > 0) {
    mollifyIntrinsic(*tuftedMesh, tuftedEdgeLengths, relativeMollificationFactor);
  }

  // Build the cover
  buildIntrinsicTuftedCover(*tuftedMesh, tuftedEdgeLengths, &tuftedGeom);

  // Flip to delaunay
  size_t nFlips = flipToDelaunay(*tuftedMesh, tuftedEdgeLengths);

  // Build the matrices
  EdgeLengthGeometry tuftedIntrinsicGeom(*tuftedMesh, tuftedEdgeLengths);
  tuftedIntrinsicGeom.requireCotanLaplacian();
  tuftedIntrinsicGeom.requireVertexLumpedMassMatrix();

  return std::make_tuple(0.5 * tuftedIntrinsicGeom.cotanLaplacian, 0.5 * tuftedIntrinsicGeom.vertexLumpedMassMatrix);
}

void buildIntrinsicTuftedCover(SurfaceMesh& mesh, EdgeData<double>& edgeLengths, EmbeddedGeometryInterface* posGeom) {

  if (posGeom) {
    posGeom->requireVertexPositions();
    posGeom->requireFaceNormals();
  }

  // == Transform the connectivity of the input in to that of the tufted cover

  // Create two copies of each input face
  HalfedgeData<Halfedge> otherSheet(mesh);
  FaceData<bool> isFront(mesh, true); // original copy serves as front
  for (Face fFront : mesh.faces()) {
    if (!isFront[fFront]) continue; // only process original faces

    // create the new face, orient it in the opposite direction
    Face fBack = mesh.duplicateFace(fFront);

    // read off the correspondence bewteen the halfedges, before inverting
    Halfedge heF = fFront.halfedge();
    Halfedge heB = fBack.halfedge();
    do {
      otherSheet[heF] = heB;
      otherSheet[heB] = heF;
      heF = heF.next();
      heB = heB.next();
    } while (heF != fFront.halfedge());

    mesh.invertOrientation(fBack);

    isFront[fBack] = false;
  }


  // Around each edge, glue back faces to front along newly created edges
  EdgeData<bool> isOrigEdge(mesh, true);
  for (Edge e : mesh.edges()) {
    if (!isOrigEdge[e]) continue;

    // Gather the original faces incident on the edge (represented by the halfedge along the edge)
    std::vector<Halfedge> edgeFaces;
    for (Halfedge he : e.adjacentHalfedges()) {
      if (isFront[he.face()]) edgeFaces.push_back(he);
    }

    // If we have normals, use them for an angular sort.
    // Otherwise, just use the "natural" (arbitrary) ordering
    if (posGeom) {

      Vector3 pTail = posGeom->vertexPositions[e.halfedge().tailVertex()];
      Vector3 pTip = posGeom->vertexPositions[e.halfedge().tipVertex()];
      Vector3 edgeVec = pTip - pTail;
      std::array<Vector3, 2> eBasis = edgeVec.buildTangentBasis();

      auto edgeAngle = [&](const Halfedge& he) {
        Vector3 oppVert = posGeom->vertexPositions[he.next().tipVertex()];
        Vector3 outVec = unit(oppVert - pTail);
        return ::std::atan2(dot(std::get<1>(eBasis), outVec), dot(std::get<0>(eBasis), outVec));
      };
      auto sortFunc = [&](const Halfedge& heA, const Halfedge& heB) -> bool { return edgeAngle(heA) > edgeAngle(heB); };

      std::sort(edgeFaces.begin(), edgeFaces.end(), sortFunc);
    }


    // Sequentially connect the faces
    Halfedge currHe = edgeFaces.front();
    if (currHe.orientation()) currHe = otherSheet[currHe];
    for (size_t i = 0; i < edgeFaces.size(); i++) {
      Halfedge nextHe = edgeFaces[(i + 1) % edgeFaces.size()];
      if (currHe.orientation() == nextHe.orientation()) nextHe = otherSheet[nextHe];
      Edge newE = mesh.separateToNewEdge(currHe, nextHe);
      isOrigEdge[newE] = false;
      edgeLengths[newE] = edgeLengths[e];
      currHe = otherSheet[nextHe];
    }
  }

  // sanity checks
  // if (mesh.hasBoundary()) throw std::runtime_error("has boundary");
  // if (!mesh.isEdgeManifold()) throw std::runtime_error("not edge manifold");
  // if (!mesh.isOriented()) throw std::runtime_error("not oriented");
}

} // namespace surface
} // namespace geometrycentral
