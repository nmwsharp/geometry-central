#include "geometrycentral/discrete_operators.h"

namespace geometrycentral {

using namespace Eigen;

// Hodge star on 0-forms. Returns a (nVerts, nVerts) matrix.
template <typename T, typename D>
geometrycentral::SparseMatrix<T> buildHodge0(Geometry<D>* geometry) {
  HalfedgeMesh* mesh = geometry->getMesh();
  VertexData<size_t> vInd = mesh->getVertexIndices();
  size_t nVerts = mesh->nVertices();

  // Reserve space in the sparse matrix
  geometrycentral::SparseMatrix<T> hodge0 =
      geometrycentral::SparseMatrix<T>(nVerts, nVerts);

  for (VertexPtr v : mesh->vertices()) {
    double primalArea = 1.0;
    double dualArea = geometry->area(v.dual());
    double ratio = dualArea / primalArea;
    size_t iV = vInd[v];
    hodge0(iV, iV) = ratio;
  }

  return hodge0;
}
template geometrycentral::SparseMatrix<double> buildHodge0<double, Euclidean>(
    Geometry<Euclidean>* geometry);
template geometrycentral::SparseMatrix<Complex> buildHodge0<Complex, Euclidean>(
    Geometry<Euclidean>* geometry);

// Hodge star on 1-forms. Returns a (nEdges, nEdges) matrix.
template <typename T, typename D>
geometrycentral::SparseMatrix<T> buildHodge1(Geometry<D>* geometry) {
  HalfedgeMesh* mesh = geometry->getMesh();
  EdgeData<size_t> eInd = mesh->getEdgeIndices();
  size_t nEdges = mesh->nEdges();

  geometrycentral::SparseMatrix<T> hodge1 =
      geometrycentral::SparseMatrix<T>(nEdges, nEdges);

  // Get the cotan weights all at once
  EdgeData<double> cotanWeights(mesh);
  geometry->getEdgeCotanWeights(cotanWeights);

  for (EdgePtr e : mesh->edges()) {
    double ratio = cotanWeights[e];
    size_t iE = eInd[e];
    hodge1(iE, iE) = ratio;
  }

  return hodge1;
}
template geometrycentral::SparseMatrix<double> buildHodge1<double, Euclidean>(
    Geometry<Euclidean>* geometry);
template geometrycentral::SparseMatrix<Complex> buildHodge1<Complex, Euclidean>(
    Geometry<Euclidean>* geometry);

// Hodge star on 2-forms. Returns a (nFaces, nFaces) matrix.
template <typename T, typename D>
geometrycentral::SparseMatrix<T> buildHodge2(Geometry<D>* geometry) {
  HalfedgeMesh* mesh = geometry->getMesh();
  FaceData<size_t> fInd = mesh->getFaceIndices();
  size_t nFaces = mesh->nFaces();

  // Reserve space in the sparse matrix
  geometrycentral::SparseMatrix<T> hodge2 =
      geometrycentral::SparseMatrix<T>(nFaces, nFaces);

  for (FacePtr f : mesh->faces()) {
    double primalArea = geometry->area(f);
    double dualArea = 1.0;
    double ratio = dualArea / primalArea;

    size_t iF = fInd[f];
    hodge2(iF, iF) = ratio;
  }

  return hodge2;
}
template geometrycentral::SparseMatrix<double> buildHodge2<double, Euclidean>(
    Geometry<Euclidean>* geometry);
template geometrycentral::SparseMatrix<Complex> buildHodge2<Complex, Euclidean>(
    Geometry<Euclidean>* geometry);

// Derivative on 0-forms. Returns a (nEdges, nVerts) matrix
template <typename T>
geometrycentral::SparseMatrix<T> buildDerivative0(HalfedgeMesh* mesh) {
  VertexData<size_t> vInd = mesh->getVertexIndices();
  EdgeData<size_t> eInd = mesh->getEdgeIndices();
  size_t nVerts = mesh->nVertices();
  size_t nEdges = mesh->nEdges();

  geometrycentral::SparseMatrix<T> d0 =
      geometrycentral::SparseMatrix<T>(nEdges, nVerts);

  for (EdgePtr e : mesh->edges()) {
    size_t iEdge = eInd[e];
    HalfedgePtr he = e.halfedge();
    VertexPtr vTail = he.vertex();
    VertexPtr vHead = he.twin().vertex();

    size_t iVHead = vInd[vHead];
    d0(iEdge, iVHead) = 1.0;

    size_t iVTail = vInd[vTail];
    d0(iEdge, iVTail) = -1.0;
  }

  return d0;
}
template geometrycentral::SparseMatrix<double> buildDerivative0<double>(
    HalfedgeMesh* mesh);
template geometrycentral::SparseMatrix<Complex> buildDerivative0<Complex>(
    HalfedgeMesh* mesh);

// Derivative on 1-forms. Returns a (nFaces, nEdges) matrix
template <typename T>
geometrycentral::SparseMatrix<T> buildDerivative1(HalfedgeMesh* mesh) {
  EdgeData<size_t> eInd = mesh->getEdgeIndices();
  FaceData<size_t> fInd = mesh->getFaceIndices();
  size_t nEdges = mesh->nEdges();
  size_t nFaces = mesh->nFaces();

  geometrycentral::SparseMatrix<T> d1 =
      geometrycentral::SparseMatrix<T>(nFaces, nEdges);

  for (FacePtr f : mesh->faces()) {
    size_t iFace = fInd[f];

    for (HalfedgePtr he : f.adjacentHalfedges()) {
      size_t iEdge = eInd[he.edge()];
      double sign = (he == he.edge().halfedge()) ? (1.0) : (-1.0);
      d1(iFace, iEdge) = sign;
    }
  }

  return d1;
}
template geometrycentral::SparseMatrix<double> buildDerivative1<double>(
    HalfedgeMesh* mesh);
template geometrycentral::SparseMatrix<Complex> buildDerivative1<Complex>(
    HalfedgeMesh* mesh);

Eigen::SparseMatrix<double> buildHodge0Eigen(Geometry<Euclidean>* geometry) {
  HalfedgeMesh* mesh = geometry->getMesh();
  VertexData<size_t> vInd = mesh->getVertexIndices();
  size_t nVerts = mesh->nVertices();

  // Reserve space in the sparse matrix
  Eigen::SparseMatrix<double> hodge0 =
      Eigen::SparseMatrix<double>(nVerts, nVerts);
  std::vector<Eigen::Triplet<double>> tripletList;

  for (VertexPtr v : mesh->vertices()) {
    double primalArea = 1.0;
    double dualArea = geometry->area(v.dual());
    double ratio = dualArea / primalArea;
    size_t iV = vInd[v];
    tripletList.emplace_back(iV, iV, ratio);
  }

  hodge0.setFromTriplets(tripletList.begin(), tripletList.end());

  return hodge0;
}

Eigen::SparseMatrix<double> buildHodge1Eigen(Geometry<Euclidean>* geometry) {
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

Eigen::SparseMatrix<double> buildHodge2Eigen(Geometry<Euclidean>* geometry) {
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

Eigen::SparseMatrix<double> buildDerivative0Eigen(HalfedgeMesh* mesh) {
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

Eigen::SparseMatrix<double> buildDerivative1Eigen(HalfedgeMesh* mesh) {
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
