#pragma once

#include "geometrycentral/halfedge_mesh.h"
#include "geometrycentral/unit_vector3.h"
#include "geometrycentral/vector2.h"
#include "geometrycentral/vector3.h"
#include "geometrycentral/dependent_quantity.h"

#include <Eigen/SparseCore>

#include <iostream>

namespace geometrycentral {


class IntrinsicGeometry {

public:

  // Constructor (doesn't do much)
  IntrinsicGeometry(HalfedgeMesh* mesh_);

  // == Members
  HalfedgeMesh* mesh = nullptr;

  // == Basic geometric quantities

  // Face areas
  inline void requireFaceAreas() { faceAreasQ.require(); }
  FaceData<double> faceAreas;

  // Vertex dual areas
  inline void requireVertexDualAreas() { vertexDualAreasQ.require(); }
  VertexData<double> vertexDualAreas;

  // Edge lengths
  inline void requireEdgeLengths() { edgeLengthsQ.require(); }
  EdgeData<double> edgeLengths;
  
  // Halfedge cotan weights
  inline void requireHalfedgeCotanWeights() { halfedgeCotanWeightsQ.require(); }
  HalfedgeData<double> halfedgeCotanWeights;

  // Edge cotan weights
  inline void requireEdgeCotanWeights() { edgeCotanWeightsQ.require(); }
  EdgeData<double> edgeCotanWeights;

  // Angle defect at vertices
  inline void requireVertexAngleDefects() { vertexAngleDefectsQ.require(); }
  VertexData<double> vertexAngleDefects;


  // == Vector fields, angles, and transport

  // The coordinate of each halfedge in the basis of he.face()
  // NOTE: These HAVE magnitude, unlike the vertex version (confusingly)
  inline void requireHalfedgeFaceCoords() { halfedgeFaceCoordsQ.require(); }
  HalfedgeData<Complex> halfedgeFaceCoords;

  // Transport an intrinsic vector field in he.face() to he.twin().face() by multiplying
  // by complex z = e^(theta I) given here
  inline void requireFaceTransportCoefs() { faceTransportCoefsQ.require(); }
  HalfedgeData<Complex> faceTransportCoefs;

  // Halfedge opposite angles
  inline void requireHalfedgeOppositeAngles() { halfedgeOppositeAnglesQ.require(); }
  HalfedgeData<double> halfedgeOppositeAngles;

  // Halfedge opposite angles (scaled by the angle defect to sum to 2 PI at each vertex)
  inline void requireHalfedgeRescaledOppositeAngles() { halfedgeRescaledOppositeAnglesQ.require(); }
  HalfedgeData<double> halfedgeRescaledOppositeAngles;

  // The coordinate of each halfedge in the basis of he.vertex(), rescaled so the sum around each vertex is 2*PI
  inline void requireHalfedgeVertexCoords() { halfedgeVertexCoordsQ.require(); }
  HalfedgeData<Complex> halfedgeVertexCoords;

  // Transport an intrinsic vector field in he.vertex() to he.twin().vertex() by multiplying
  // by complex z = e^(theta I) given here
  inline void requireVertexTransportCoefs() { vertexTransportCoefsQ.require(); }
  HalfedgeData<Complex> vertexTransportCoefs;


  // == Indices

  inline void requireVertexIndices() { vertexIndicesQ.require(); }
  VertexData<size_t> vertexIndices;

  inline void requireInteriorVertexIndices() { interiorVertexIndicesQ.require(); }
  VertexData<size_t> interiorVertexIndices;

  inline void requireFaceIndices() { faceIndicesQ.require(); }
  FaceData<size_t> faceIndices;

  inline void requireEdgeIndices() { edgeIndicesQ.require(); }
  EdgeData<size_t> edgeIndices;

  inline void requireHalfedgeIndices() { halfedgeIndicesQ.require(); }
  HalfedgeData<size_t> halfedgeIndices;


  // == Operators
  // Note: These don't quite follow the usual naming scheme, for the sake of grouping common operators
  // TODO factorizations?

  // All of the basic DEC operators and their inverses
  inline void requireBasicDECOperators() { basicDECOperatorsQ.require(); }
  Eigen::SparseMatrix<double> d0, d1, hodge0, hodge1, hodge2;
  Eigen::SparseMatrix<double> hodge0Inv, hodge1Inv, hodge2Inv;

  // Cotan-laplace operator
  // Remember, this DOES NOT include the mass matrix (hodge0)
  inline void requireZeroFormWeakLaplacian() { zeroFormWeakLaplacianQ.require(); }
  Eigen::SparseMatrix<double> zeroFormWeakLaplacian;

protected:
  std::vector<DependentQuantity*> allQuantities;

  // === Internal interface for all quantities
  virtual void buildDependencies();

  // == Basic geometric quantities

  DependentQuantity faceAreasQ;
  virtual void computeFaceAreas();

  DependentQuantity vertexDualAreasQ;
  virtual void computeVertexDualAreas();

  DependentQuantity edgeLengthsQ;
  virtual void computeEdgeLengths() = 0;

  DependentQuantity halfedgeCotanWeightsQ;
  virtual void computeHalfedgeCotanWeights();

  DependentQuantity edgeCotanWeightsQ;
  virtual void computeEdgeCotanWeights();

  DependentQuantity vertexAngleDefectsQ;
  virtual void computeVertexAngleDefects();

  // == Vector fields, angles, and transport

  DependentQuantity halfedgeFaceCoordsQ;
  virtual void computeHalfedgeFaceCoords();

  DependentQuantity faceTransportCoefsQ;
  virtual void computeFaceTransportCoefs();

  DependentQuantity halfedgeOppositeAnglesQ;
  virtual void computeHalfedgeOppositeAngles();

  DependentQuantity halfedgeRescaledOppositeAnglesQ;
  virtual void computeHalfedgeRescaledOppositeAngles();

  DependentQuantity halfedgeVertexCoordsQ;
  virtual void computeHalfedgeVertexCoords();

  DependentQuantity vertexTransportCoefsQ;
  virtual void computeVertexTransportCoefs();


  // == Indices

  DependentQuantity vertexIndicesQ;
  virtual void computeVertexIndices();

  DependentQuantity interiorVertexIndicesQ;
  virtual void computeInteriorVertexIndices();

  DependentQuantity faceIndicesQ;
  virtual void computeFaceIndices();

  DependentQuantity edgeIndicesQ;
  virtual void computeEdgeIndices();

  DependentQuantity halfedgeIndicesQ;
  virtual void computeHalfedgeIndices();

  // == Operators

  DependentQuantity basicDECOperatorsQ;
  virtual void computeBasicDECOperators();

  DependentQuantity zeroFormWeakLaplacianQ;
  virtual void computeZeroFormWeakLaplacian();

  // == Helpers

  // Throws an error if mesh is not triangular
  void verifyTriangular(HalfedgeMesh* m);

};

} // namespace geometrycentral
