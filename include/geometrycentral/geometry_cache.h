#pragma once

#include "geometrycentral/halfedge_mesh.h"
#include "geometrycentral/unit_vector3.h"
#include "geometrycentral/vector2.h"
#include "geometrycentral/vector3.h"

#include <functional>
#include <vector>

#include <Eigen/SparseCore>

namespace geometrycentral {

// Helper class which manages a dependency graph of quantities
class DependentQuantity {

public:
  DependentQuantity(){};

  DependentQuantity(std::vector<DependentQuantity*> dependencies_, std::function<void()> evaluateFunc_)
      : dependencies(dependencies_), evaluateFunc(evaluateFunc_) {}

  std::vector<DependentQuantity*> dependencies;

  // Compute the quantity, if we don't have it already
  void ensureHave();

  // Compute the quantity if we need it and don't have it already
  void ensureHaveIfRequired();

  // Note that something will reqiure this quantity (increments a count of such requirements),
  // and ensure that we have this quantity
  void require();

  // Decrement the count of requirements of this quantity
  void unrequire();

  bool computed = false;
  int requireCount = 0;

  std::function<void()> evaluateFunc;
};

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

  // == Internal interface for all quantities

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

  DependentQuantity faceBasesQ;
  void computeFaceBases();

  DependentQuantity faceTransportCoefsQ;
  void computeFaceTransportCoefs();

  DependentQuantity dihedralAnglesQ;
  void computeDihedralAngles();

  DependentQuantity halfedgeCotanWeightsQ;
  void computeHalfedgeCotanWeights();
  
  DependentQuantity edgeCotanWeightsQ;
  void computeEdgeCotanWeights();


  // Indices

  DependentQuantity vertexIndicesQ;
  void computeVertexIndices();

  DependentQuantity faceIndicesQ;
  void computeFaceIndices();

  DependentQuantity edgeIndicesQ;
  void computeEdgeIndices();

  DependentQuantity halfedgeIndicesQ;
  void computeHalfedgeIndices();

  // Operators

  DependentQuantity basicDECOperatorsQ;
  void computeBasicDECOperators();

  DependentQuantity modifiedDECOperatorsQ;
  void computeModifiedDECOperators();

  DependentQuantity zeroFormWeakLaplacianQ;
  void computeZeroFormWeakLaplacian();


public:
  // face area normals
  // vector which points in the normal direction and has magnitude equal to area of face
  inline void requireFaceAreaNormals() { faceAreaNormalsQ.require(); }
  FaceData<Vector3> faceAreaNormals;

  // face areas
  inline void requireFaceAreas() { faceAreasQ.require(); }
  FaceData<double> faceAreas;

  // face normals
  inline void requireFaceNormals() { faceNormalsQ.require(); }
  FaceData<Vector3> faceNormals;

  // vertex normals
  inline void requireVertexNormals() { vertexNormalsQ.require(); }
  VertexData<Vector3> vertexNormals;

  // vertex dual areas
  inline void requireVertexDualAreas() { vertexDualAreasQ.require(); }
  VertexData<double> vertexDualAreas;
  
  // halfedge cotans 
  inline void requireHalfedgeVectors() { halfedgeVectorsQ.require(); }
  HalfedgeData<Vector3> halfedgeVectors;

  // edge lengths
  inline void requireEdgeLengths() { edgeLengthsQ.require(); }
  EdgeData<double> edgeLengths;

  // extrinsic basis vector pair in each face
  inline void requireFaceBases() { faceBasesQ.require(); }
  FaceData<std::array<Vector3, 2>> faceBases;

  // transport angle in he.face() to he.twin().face() by multiplying
  // complex e^(theta I)
  inline void requireFaceTransportCoefs() { faceTransportCoefsQ.require(); }
  HalfedgeData<Complex> faceTransportCoefs;

  // dihedral angle
  inline void requireDihedralAngles() { dihedralAnglesQ.require(); }
  EdgeData<double> dihedralAngles;
  
  // halfedge cotans 
  inline void requireHalfedgeCotanWeights() { halfedgeCotanWeightsQ.require(); }
  HalfedgeData<double> halfedgeCotanWeights;

  // edge cotan weights
  inline void requireEdgeCotanWeights() { edgeCotanWeightsQ.require(); }
  EdgeData<double> edgeCotanWeights;

  // === Indices

  inline void requireVertexIndices() { vertexIndicesQ.require(); }
  VertexData<size_t> vertexIndices;

  inline void requireFaceIndices() { faceIndicesQ.require(); }
  FaceData<size_t> faceIndices;

  inline void requireEdgeIndices() { edgeIndicesQ.require(); }
  EdgeData<size_t> edgeIndices;

  inline void requireHalfedgeIndices() { halfedgeIndicesQ.require(); }
  HalfedgeData<size_t> halfedgeIndices;


  // === Operators
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
