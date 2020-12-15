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
  
  // Vertex mean curvature
  VertexData<double> vertexMeanCurvatures;
  void requireVertexMeanCurvatures();
  void unrequireVertexMeanCurvatures();

  // Vertex min principal curvature
  VertexData<double> vertexMinPrincipalCurvatures;
  void requireVertexMinPrincipalCurvatures();
  void unrequireVertexMinPrincipalCurvatures();

  // Vertex max principal curvature
  VertexData<double> vertexMaxPrincipalCurvatures;
  void requireVertexMaxPrincipalCurvatures();
  void unrequireVertexMaxPrincipalCurvatures();

  // Vertex principal curvature direction
  VertexData<Vector2> vertexPrincipalCurvatureDirections;
  void requireVertexPrincipalCurvatureDirections();
  void unrequireVertexPrincipalCurvatureDirections();

  // Face principal curvature direction
  FaceData<Vector2> facePrincipalCurvatureDirections;
  void requireFacePrincipalCurvatureDirections();
  void unrequireFacePrincipalCurvatureDirections();

protected:
  // Edge dihedral angle
  DependentQuantityD<EdgeData<double>> edgeDihedralAnglesQ;
  virtual void computeEdgeDihedralAngles() = 0;

  // Vertex mean curvature
  DependentQuantityD<VertexData<double>> vertexMeanCurvaturesQ;
  virtual void computeVertexMeanCurvatures();

  // Vertex min principal curvature
  DependentQuantityD<VertexData<double>> vertexMinPrincipalCurvaturesQ;
  virtual void computeVertexMinPrincipalCurvatures();

  // Vertex max principal curvature
  DependentQuantityD<VertexData<double>> vertexMaxPrincipalCurvaturesQ;
  virtual void computeVertexMaxPrincipalCurvatures();

  virtual void computePrincipalCurvatures( int whichCurvature, VertexData<double>& kappa );

  // Vertex principal curvature direction
  DependentQuantityD<VertexData<Vector2>> vertexPrincipalCurvatureDirectionsQ;
  virtual void computeVertexPrincipalCurvatureDirections();

  // Face principal curvature direction
  DependentQuantityD<FaceData<Vector2>> facePrincipalCurvatureDirectionsQ;
  virtual void computeFacePrincipalCurvatureDirections();
};

} // namespace surface
} // namespace geometrycentral
