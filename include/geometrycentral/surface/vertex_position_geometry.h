#pragma once

#include "geometrycentral/surface/embedded_geometry_interface.h"
#include "geometrycentral/surface/halfedge_mesh.h"

#include <Eigen/SparseCore>

namespace geometrycentral {
namespace surface {


class VertexPositionGeometry : public EmbeddedGeometryInterface {

public:
  // Construct empty -- all positions initially set to the origin
  VertexPositionGeometry(HalfedgeMesh& mesh_);

  // Construct from positions
  VertexPositionGeometry(HalfedgeMesh& mesh_, VertexData<Vector3>& inputVertexPositions);

  // Boring destructor
  virtual ~VertexPositionGeometry() {}

  // Construct a new geometry which is exactly the same as this one, on the same mesh.
  // This is a deep copy, no quantites are shared, etc. Require counts/computed quantities are not copied.
  std::unique_ptr<VertexPositionGeometry> copy();

  // Construct a new geometry which is exactly the same as this one, on another mesh.
  // This is a deep copy, no quantites are shared, etc. Require counts/computed quantities are not copied.
  // The meshes must be in correspondence (have the same connectivity).
  std::unique_ptr<VertexPositionGeometry> reinterpretTo(HalfedgeMesh& targetMesh);


  // == Members

  // The actual input data which defines the geometry
  VertexData<Vector3> inputVertexPositions;

  // == Immediates
  double edgeLength(Edge e) const;
  double faceArea(Face f) const;
  double cornerAngle(Corner c) const;
  double halfedgeCotanWeight(Halfedge he) const;
  double edgeCotanWeight(Edge e) const;
  Vector3 faceNormal(Face f) const;


protected:
  // Override the compute vertex positions method for embedded geometry
  virtual void computeVertexPositions() override;


private:
};

} // namespace surface
} // namespace geometrycentral

#include "geometrycentral/surface/vertex_position_geometry.ipp"
