#pragma once

#include "geometrycentral/geometry.h"

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


class GeometryCache {

public:
  // Create a new cache
  GeometryCache(Geometry<Euclidean>* geometry_);

  // Clear out the cache and immediately recompute all required quantities. Should be called after modifying mesh or
  // geometry.
  void repopulate();

  // Get the mesh
  HalfedgeMesh* getMesh();

  // Get the geometry
  Geometry<Euclidean>* getGeometry();

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
  //   - VertexData<double> XXRaw; (private)
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
  Geometry<Euclidean>* geometry;

  std::vector<DependentQuantity*> allQuantities;

  // == Internal interface for all quantities

  FaceData<Vector3> faceAreaNormalsRaw;
  void computeFaceAreaNormals();

  FaceData<double> faceAreasRaw;
  void computeFaceAreas();

  FaceData<Vector3> faceNormalsRaw;
  void computeFaceNormals();

  VertexData<Vector3> vertexNormalsRaw;
  void computeVertexNormals();

  VertexData<double> vertexDualAreasRaw;
  void computeVertexDualAreas();

public:
  // face area normals
  // vector which points in the normal direction and has magnitude equal to area of face
  const FaceData<Vector3> faceAreaNormals = faceAreaNormalsRaw;
  DependentQuantity faceAreaNormalsQ;

  // face areas
  const FaceData<double> faceAreas = faceAreasRaw;
  DependentQuantity faceAreasQ;

  // face normals
  const FaceData<Vector3> faceNormals = faceNormalsRaw;
  DependentQuantity faceNormalsQ;

  // vertex normals
  const VertexData<Vector3> vertexNormals = vertexNormalsRaw;
  DependentQuantity vertexNormalsQ;

  // vertex dual areas
  const VertexData<double> vertexDualAreas = vertexDualAreasRaw;
  DependentQuantity vertexDualAreasQ;
};

}; // namespace geometrycentral