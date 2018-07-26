#pragma once

#include "geometrycentral/dependent_quantity.h"
#include "geometrycentral/halfedge_mesh.h"
#include "geometrycentral/unit_vector3.h"
#include "geometrycentral/vector2.h"
#include "geometrycentral/vector3.h"

#include <functional>
#include <vector>

#include <Eigen/SparseCore>

namespace geometrycentral {


template <typename G>
class GeometryCache {

public:
  // Create a new cache
  GeometryCache(Geometry<G>* geometry_);

  // Clear out the cache and immediately recompute all required quantities. Should be called after modifying mesh or
  // geometry.
  void repopulate();

  // Get the mesh
  HalfedgeMesh* getMesh();

  // Get the geometry
  Geometry<G>* getGeometry();

  // Hide copy and move constructors, they are unlikely to be what is actually wanted and shouldn't be used accidentally
  GeometryCache(const GeometryCache& other) = delete;
  GeometryCache& operator=(const GeometryCache& other) = delete;
  GeometryCache(GeometryCache&& other) = delete;
  GeometryCache& operator=(GeometryCache&& other) = delete;

  // === Cached quantities
  // For each quantity XYZ, there is a method gc.requireXYZ() (in camel case). The data can then be accessed via
  // gc.XYZ[v]; Typical usage then looks like:
  //
  //     gc.requireVertexNormals();
  //     ...
  //     ...
  //     Vector3 n = gc.vertexNormals[v];
  //

  // To implement a new quantity, follow the pattern below. For a hypothetical quantity called 'XX', the necessary
  // quantities are, for example
  //   - void computeXX(); (private)
  //   - VertexData<double> XX; (public)
  //   - void requireXX(); (public)
  //   - DependentQuantity XXQ{); (private)
  //     |--> initialize in GeometryCache constructor
  //
  // Note: the ordering in which the quantities are listed is used to implicitly encode the DAG amongst the
  // dependencies... a quantity may only depend on those which appear above it!

  // (note, more public members below)
private:
  // The mesh for which this cache applies
  HalfedgeMesh* mesh;
  Geometry<G>* geometry;

  std::vector<DependentQuantity*> allQuantities;

  // === Internal interface for all quantities

  // == Basic geometric quantities

  DependentQuantity faceAreaNormalsQ;
  void computeFaceAreaNormals();

  DependentQuantity faceAreasQ;
  void computeFaceAreas();

  DependentQuantity faceNormalsQ;
  void computeFaceNormals();

  DependentQuantity vertexNormalsQ;
  void computeVertexNormals();

  DependentQuantity vertexDualAreasQ;
  void computeVertexDualAreas();

  DependentQuantity halfedgeVectorsQ;
  void computeHalfedgeVectors();

  DependentQuantity edgeLengthsQ;
  void computeEdgeLengths();

  DependentQuantity dihedralAnglesQ;
  void computeDihedralAngles();

  DependentQuantity halfedgeCotanWeightsQ;
  void computeHalfedgeCotanWeights();

  DependentQuantity edgeCotanWeightsQ;
  void computeEdgeCotanWeights();

  DependentQuantity vertexAngleDefectsQ;
  void computeVertexAngleDefects();

  // == Vector fields, angles, and transport

  DependentQuantity faceBasesQ;
  void computeFaceBases();

  DependentQuantity vertexBasesQ;
  void computeVertexBases();

  DependentQuantity halfedgeFaceCoordsQ;
  void computeHalfedgeFaceCoords();

  DependentQuantity faceTransportCoefsQ;
  void computeFaceTransportCoefs();

  DependentQuantity halfedgeOppositeAnglesQ;
  void computeHalfedgeOppositeAngles();

  DependentQuantity halfedgeRescaledOppositeAnglesQ;
  void computeHalfedgeRescaledOppositeAngles();

  DependentQuantity halfedgeVertexCoordsQ;
  void computeHalfedgeVertexCoords();

  DependentQuantity vertexTransportCoefsQ;
  void computeVertexTransportCoefs();

  DependentQuantity vertexFaceTransportCoefsQ;
  void computeVertexFaceTransportCoefs();

  DependentQuantity principalDirectionsQ;
  void computePrincipalDirections();


  // == Indices

  DependentQuantity vertexIndicesQ;
  void computeVertexIndices();

  DependentQuantity interiorVertexIndicesQ;
  void computeInteriorVertexIndices();

  DependentQuantity faceIndicesQ;
  void computeFaceIndices();

  DependentQuantity edgeIndicesQ;
  void computeEdgeIndices();

  DependentQuantity halfedgeIndicesQ;
  void computeHalfedgeIndices();

  // == Operators

  DependentQuantity basicDECOperatorsQ;
  void computeBasicDECOperators();

  DependentQuantity modifiedDECOperatorsQ;
  void computeModifiedDECOperators();

  DependentQuantity zeroFormWeakLaplacianQ;
  void computeZeroFormWeakLaplacian();


public:
  // == Basic geometric quantities

  // Face area normals
  // vector which points in the normal direction and has magnitude equal to area of face
  inline void requireFaceAreaNormals() { faceAreaNormalsQ.require(); }
  FaceData<Vector3> faceAreaNormals;

  // Face areas
  inline void requireFaceAreas() { faceAreasQ.require(); }
  FaceData<double> faceAreas;

  // Face normals
  inline void requireFaceNormals() { faceNormalsQ.require(); }
  FaceData<Vector3> faceNormals;

  // Vertex normals
  inline void requireVertexNormals() { vertexNormalsQ.require(); }
  VertexData<Vector3> vertexNormals;

  // Vertex dual areas
  inline void requireVertexDualAreas() { vertexDualAreasQ.require(); }
  VertexData<double> vertexDualAreas;

  // Halfedge cotans
  inline void requireHalfedgeVectors() { halfedgeVectorsQ.require(); }
  HalfedgeData<Vector3> halfedgeVectors;

  // Edge lengths
  inline void requireEdgeLengths() { edgeLengthsQ.require(); }
  EdgeData<double> edgeLengths;

  // Dihedral angle
  inline void requireDihedralAngles() { dihedralAnglesQ.require(); }
  EdgeData<double> dihedralAngles;

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

  // Extrinsic basis vector pair in each face
  inline void requireFaceBases() { faceBasesQ.require(); }
  FaceData<std::array<Vector3, 2>> faceBases;

  // Extrinsic basis vector pair in each vertex's tangent plane
  inline void requireVertexBases() { vertexBasesQ.require(); }
  VertexData<std::array<Vector3, 2>> vertexBases;

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

  // Transport an intrinsic vector field in he.vertex() to he.face() by multiplying
  // by complex z = e^(theta I) given here
  inline void requireVertexFaceTransportCoefs() { vertexFaceTransportCoefsQ.require(); }
  HalfedgeData<Complex> vertexFaceTransportCoefs;

  // The two-symmetric vector field encoding the principal directions and their strength
  inline void requirePrincipalDirections() { principalDirectionsQ.require(); }
  VertexData<Complex> principalDirections;


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

  // All of the basic DEC operators
  inline void requireBasicDECOperators() { basicDECOperatorsQ.require(); }
  Eigen::SparseMatrix<double> d0, d1, hodge0, hodge1, hodge2;

  // Includes inverses
  inline void requireModifiedDECOperators() { modifiedDECOperatorsQ.require(); }
  Eigen::SparseMatrix<double> hodge0Inv, hodge1Inv, hodge2Inv;

  // Cotan-laplace operator
  // Remember, this DOES NOT include the mass matrix (hodge0)
  inline void requireZeroFormWeakLaplacian() { zeroFormWeakLaplacianQ.require(); }
  Eigen::SparseMatrix<double> zeroFormWeakLaplacian;
};

}; // namespace geometrycentral
