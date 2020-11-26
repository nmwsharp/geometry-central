#pragma once

#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/utilities/vector2.h"

#include <Eigen/SparseCore>

namespace geometrycentral {
namespace surface {


class ExtrinsicGeometryInterface : public IntrinsicGeometryInterface {

protected:
  // Constructor is protected, because this class is an interface which is not meant to be instantiated directly.
  // Instantiate it via some realization which encapsulates input data, like VertexPositionGeometry.
  ExtrinsicGeometryInterface(SurfaceMesh& mesh_);

public:
  virtual ~ExtrinsicGeometryInterface() {}

  // Edge dihedral angle
  EdgeData<double> edgeDihedralAngles;
  void requireEdgeDihedralAngles();
  void unrequireEdgeDihedralAngles();

  // Vertex principal curvature
  VertexData<Vector2> vertexPrincipalCurvatureDirections;
  void requireVertexPrincipalCurvatureDirections();
  void unrequireVertexPrincipalCurvatureDirections();

  // Face principal curvature
  FaceData<Vector2> facePrincipalCurvatureDirections;
  void requireFacePrincipalCurvatureDirections();
  void unrequireFacePrincipalCurvatureDirections();

protected:
  // Edge dihedral angle
  DependentQuantityD<EdgeData<double>> edgeDihedralAnglesQ;
  virtual void computeEdgeDihedralAngles() = 0;

  // Vertex principal curvature
  DependentQuantityD<VertexData<Vector2>> vertexPrincipalCurvatureDirectionsQ;
  virtual void computeVertexPrincipalCurvatureDirections();

  // Face principal curvature
  DependentQuantityD<FaceData<Vector2>> facePrincipalCurvatureDirectionsQ;
  virtual void computeFacePrincipalCurvatureDirections();
};

} // namespace surface
} // namespace geometrycentral
