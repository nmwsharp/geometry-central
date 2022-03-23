#pragma once

#include "geometrycentral/surface/embedded_geometry_interface.h"
#include "geometrycentral/surface/surface_mesh.h"

#include <Eigen/SparseCore>

namespace geometrycentral {
namespace surface {


class VertexPositionGeometry : public EmbeddedGeometryInterface {

public:
  // Construct empty -- all positions initially set to the origin
  VertexPositionGeometry(SurfaceMesh& mesh_);

  // Construct from positions
  VertexPositionGeometry(SurfaceMesh& mesh_, const VertexData<Vector3>& inputVertexPositions);

  // Construct from positions (stored in an Eigen matrix)
  template <typename T>
  VertexPositionGeometry(SurfaceMesh& mesh_, const Eigen::MatrixBase<T>& vertexPositions);

  // Boring destructor
  virtual ~VertexPositionGeometry() {}

  // Construct a new geometry which is exactly the same as this one, on the same mesh.
  // This is a deep copy, no quantites are shared, etc. Require counts/computed quantities are not copied.
  std::unique_ptr<VertexPositionGeometry> copy();

  // Construct a new geometry which is exactly the same as this one, on another mesh.
  // This is a deep copy, no quantites are shared, etc. Require counts/computed quantities are not copied.
  // The meshes must be in correspondence (have the same connectivity).
  std::unique_ptr<VertexPositionGeometry> reinterpretTo(SurfaceMesh& targetMesh);


  // == Members

  // The actual input data which defines the geometry
  // In a previous version of the library, this was a distinct field which got copied in to `vertexPositions`. However,
  // now they are simply aliases for the same buffer.
  VertexData<Vector3>& inputVertexPositions;

  // == Immediates
  double edgeLength(Edge e) const;
  double faceArea(Face f) const;
  double vertexDualArea(Vertex v) const;
  double cornerAngle(Corner c) const;
  double halfedgeCotanWeight(Halfedge he) const;
  double edgeCotanWeight(Edge e) const;
  Vector3 faceNormal(Face f) const;
  Vector3 halfedgeVector(Halfedge he) const;
  double edgeDihedralAngle(Edge e) const;
  double vertexMeanCurvature(Vertex v) const;
  double vertexGaussianCurvature(Vertex v) const;
  double vertexMinPrincipalCurvature(Vertex v) const;
  double vertexMaxPrincipalCurvature(Vertex v) const;
  Vector3 vertexDualMeanCurvatureNormal(Vertex v) const;

protected:
  // Override the compute vertex positions method for embedded geometry
  virtual void computeVertexPositions() override;

  double vertexPrincipalCurvature(int whichCurvature, Vertex v) const;

private:
};

} // namespace surface
} // namespace geometrycentral

#include "geometrycentral/surface/vertex_position_geometry.ipp"
