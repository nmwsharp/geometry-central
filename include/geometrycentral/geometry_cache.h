#pragma once

#include "geometrycentral/halfedge_mesh.h"
#include "geometrycentral/unit_vector3.h"
#include "geometrycentral/vector2.h"
#include "geometrycentral/vector3.h"

#include <functional>
#include <vector>

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
  // Convention is that the dependent quantity management object is postfixed by 'Q'.
  // Typical usage then looks like:
  //
  //     gc.vertexNormalsQ.require();
  //     ...
  //     ...
  //     Vector3 n = gc.vertexNormals[v];
  //

  // To implement a new quantity, follow the pattern below. For a hypothetical quantity called 'XX', the necessary
  // quantities are, for example
  //   - void computeXX(); (private)
  //   - const VertexData<double> XX = XXRaw; (public)
  //   - DependentQuantity XXQ{); (public)
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
  void computeFaceAreaNormals();
  void computeFaceAreas();
  void computeFaceNormals();
  void computeVertexNormals();
  void computeVertexDualAreas();
  void computeEdgeLengths();
  void computeFaceBasis();
  void computeFaceTransportCoefs();
  void computeDihedralAngle();

public:
  // face area normals
  // vector which points in the normal direction and has magnitude equal to area of face
  FaceData<Vector3> faceAreaNormals;
  DependentQuantity faceAreaNormalsQ;

  // face areas
  FaceData<double> faceAreas;
  DependentQuantity faceAreasQ;

  // face normals
  FaceData<Vector3> faceNormals;
  DependentQuantity faceNormalsQ;

  // vertex normals
  VertexData<Vector3> vertexNormals;
  DependentQuantity vertexNormalsQ;

  // vertex dual areas
  VertexData<double> vertexDualAreas;
  DependentQuantity vertexDualAreasQ;

  // edge lengths
  EdgeData<double> edgeLengths;
  DependentQuantity edgeLengthsQ;
  
  // extrinsic basis vector pair in each face 
  FaceData<std::array<Vector3,2>> faceBasis;
  DependentQuantity faceBasisQ;
  
  // edge lengths
  HalfedgeData<Complex> faceTransportCoefs;
  DependentQuantity faceTransportCoefsQ;
  
  // dihedral angle 
  EdgeData<double> dihedralAngle;
  DependentQuantity dihedralAngleQ;
};

}; // namespace geometrycentral