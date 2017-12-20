#include "geometrycentral/discrete_operators.h"

using namespace Eigen;

namespace geometrycentral {


Eigen::SparseMatrix<double> buildHodge0(Geometry<Euclidean>* geometry) {
  HalfedgeMesh* mesh = geometry->getMesh();
  VertexData<size_t> vInd = mesh->getVertexIndices();
  size_t nVerts = mesh->nVertices();

  // Reserve space in the sparse matrix
  Eigen::SparseMatrix<double> hodge0 =
      Eigen::SparseMatrix<double>(nVerts, nVerts);
  std::vector<Eigen::Triplet<double>> tripletList;

  for (VertexPtr v : mesh->vertices()) {
    double primalArea = 1.0;
    double dualArea = geometry->dualArea(v);
    double ratio = dualArea / primalArea;
    size_t iV = vInd[v];
    tripletList.emplace_back(iV, iV, ratio);
  }

  hodge0.setFromTriplets(tripletList.begin(), tripletList.end());

  return hodge0;
}

Eigen::SparseMatrix<double> buildHodge1(Geometry<Euclidean>* geometry) {
  HalfedgeMesh* mesh = geometry->getMesh();
  EdgeData<size_t> eInd = mesh->getEdgeIndices();
  size_t nEdges = mesh->nEdges();

  Eigen::SparseMatrix<double> hodge1 =
      Eigen::SparseMatrix<double>(nEdges, nEdges);
  std::vector<Eigen::Triplet<double>> tripletList;

  // Get the cotan weights all at once
  EdgeData<double> cotanWeights(mesh);
  geometry->getEdgeCotanWeights(cotanWeights);

  for (EdgePtr e : mesh->edges()) {
    double ratio = cotanWeights[e];
    size_t iE = eInd[e];
    tripletList.emplace_back(iE, iE, ratio);
  }

  hodge1.setFromTriplets(tripletList.begin(), tripletList.end());
  return hodge1;
}

Eigen::SparseMatrix<double> buildHodge2(Geometry<Euclidean>* geometry) {
  HalfedgeMesh* mesh = geometry->getMesh();
  FaceData<size_t> fInd = mesh->getFaceIndices();
  size_t nFaces = mesh->nFaces();

  Eigen::SparseMatrix<double> hodge2 =
      Eigen::SparseMatrix<double>(nFaces, nFaces);
  std::vector<Eigen::Triplet<double>> tripletList;

  for (FacePtr f : mesh->faces()) {
    double primalArea = geometry->area(f);
    double dualArea = 1.0;
    double ratio = dualArea / primalArea;

    size_t iF = fInd[f];
    tripletList.emplace_back(iF, iF, ratio);
  }

  hodge2.setFromTriplets(tripletList.begin(), tripletList.end());
  return hodge2;
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

}  // namespace geometrycentral
