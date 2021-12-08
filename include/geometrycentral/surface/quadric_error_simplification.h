#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/mutation_manager.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "geometrycentral/utilities/elementary_geometry.h"

#include <queue>

namespace geometrycentral {
namespace surface {
class Quadric {
public:
  Quadric();
  Quadric(const Eigen::Matrix3d& A_, const Eigen::Vector3d& b_, double c_);
  Quadric(const Quadric& Q1, const Quadric& Q2);
  double cost(const Eigen::Vector3d& v);
  Eigen::Vector3d optimalPoint();

  Quadric operator+=(const Quadric& Q);

protected:
  Eigen::Matrix3d A;
  Eigen::Vector3d b;
  double c;

  friend Quadric operator+(const Quadric& Q1, const Quadric& Q2);
};

Quadric operator+(const Quadric& Q1, const Quadric& Q2);

void quadricErrorSimplify(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geo, double tol = 0.05);
void quadricErrorSimplify(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geo, double tol, MutationManager& mm);

} // namespace surface
} // namespace geometrycentral
