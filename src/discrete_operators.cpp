#include "geometrycentral/discrete_operators.h"

using namespace Eigen;

namespace geometrycentral {


Eigen::DiagonalMatrix<double, Eigen::Dynamic> buildHodge0(Geometry<Euclidean>* geometry) {
  HalfedgeMesh* mesh = geometry->getMesh();
  VertexData<size_t> vInd = mesh->getVertexIndices();
  size_t nVerts = mesh->nVertices();

  GeometryCache<Euclidean>& gc = geometry->cache;
  gc.requireVertexDualAreas();

  // Reserve space in the sparse matrix
  Eigen::VectorXd hodge0(nVerts);

  for (VertexPtr v : mesh->vertices()) {
    double primalArea = 1.0;
    double dualArea = gc.vertexDualAreas[v];
    double ratio = dualArea / primalArea;
    size_t iV = vInd[v];
    hodge0[iV] = ratio;
  }

  return hodge0.asDiagonal();
}

Eigen::DiagonalMatrix<double, Eigen::Dynamic> buildHodge1(Geometry<Euclidean>* geometry) {
  HalfedgeMesh* mesh = geometry->getMesh();
  EdgeData<size_t> eInd = mesh->getEdgeIndices();
  size_t nEdges = mesh->nEdges();

  GeometryCache<Euclidean>& gc = geometry->cache;
  gc.requireEdgeCotanWeights();

  Eigen::VectorXd hodge1(nEdges);

  for (EdgePtr e : mesh->edges()) {
    double ratio = gc.edgeCotanWeights[e];
    size_t iE = eInd[e];
    hodge1[iE] = ratio;
  }

  return hodge1.asDiagonal();
}

Eigen::DiagonalMatrix<double, Eigen::Dynamic> buildHodge2(Geometry<Euclidean>* geometry) {
  HalfedgeMesh* mesh = geometry->getMesh();
  FaceData<size_t> fInd = mesh->getFaceIndices();
  size_t nFaces = mesh->nFaces();
  
  GeometryCache<Euclidean>& gc = geometry->cache;
  gc.requireFaceAreas();

  Eigen::VectorXd hodge2(nFaces);

  for (FacePtr f : mesh->faces()) {
    double primalArea = gc.faceAreas[f];
    double dualArea = 1.0;
    double ratio = dualArea / primalArea;

    size_t iF = fInd[f];
    hodge2[iF] = ratio;
  }

  return hodge2.asDiagonal();
}

Eigen::SparseMatrix<double> buildDerivative0(HalfedgeMesh* mesh) {
  VertexData<size_t> vInd = mesh->getVertexIndices();
  EdgeData<size_t> eInd = mesh->getEdgeIndices();
  size_t nVerts = mesh->nVertices();
  size_t nEdges = mesh->nEdges();

  Eigen::SparseMatrix<double> d0 = Eigen::SparseMatrix<double>(nEdges, nVerts);
  std::vector<Eigen::Triplet<double>> tripletList;

  for (EdgePtr e : mesh->edges()) {
    size_t iEdge = eInd[e];
    HalfedgePtr he = e.halfedge();
    VertexPtr vTail = he.vertex();
    VertexPtr vHead = he.twin().vertex();

    size_t iVHead = vInd[vHead];
    tripletList.emplace_back(iEdge, iVHead, 1.0);

    size_t iVTail = vInd[vTail];
    tripletList.emplace_back(iEdge, iVTail, -1.0);
  }

  d0.setFromTriplets(tripletList.begin(), tripletList.end());
  return d0;
}

Eigen::SparseMatrix<double> buildDerivative1(HalfedgeMesh* mesh) {
  EdgeData<size_t> eInd = mesh->getEdgeIndices();
  FaceData<size_t> fInd = mesh->getFaceIndices();
  size_t nEdges = mesh->nEdges();
  size_t nFaces = mesh->nFaces();

  Eigen::SparseMatrix<double> d1 = Eigen::SparseMatrix<double>(nFaces, nEdges);
  std::vector<Eigen::Triplet<double>> tripletList;

  for (FacePtr f : mesh->faces()) {
    size_t iFace = fInd[f];

    for (HalfedgePtr he : f.adjacentHalfedges()) {
      size_t iEdge = eInd[he.edge()];
      double sign = (he == he.edge().halfedge()) ? (1.0) : (-1.0);
      tripletList.emplace_back(iFace, iEdge, sign);
    }
  }

  d1.setFromTriplets(tripletList.begin(), tripletList.end());
  return d1;
}

} // namespace geometrycentral
