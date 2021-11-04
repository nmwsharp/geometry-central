#include "geometrycentral/surface/quadric_error_simplification.h"

namespace geometrycentral {
namespace surface {

Quadric::Quadric() {
  A = Eigen::Matrix3d::Zero();
  b = Eigen::Vector3d::Zero();
  c = 0;
}

Quadric::Quadric(const Quadric& Q1, const Quadric& Q2) {
  A = Q1.A + Q2.A;
  b = Q1.b + Q2.b;
  c = Q1.c + Q2.c;
}

Quadric::Quadric(const Eigen::Matrix3d& A_, const Eigen::Vector3d& b_, double c_) : A(A_), b(b_), c(c_) {}

double Quadric::cost(const Eigen::Vector3d& v) { return v.dot(A * v) + 2 * b.dot(v) + c; }

Eigen::Vector3d Quadric::optimalPoint() { return -A.inverse() * b; }

Quadric Quadric::operator+=(const Quadric& Q) {
  A += Q.A;
  b += Q.b;
  c += Q.c;

  return *this;
}

Quadric operator+(const Quadric& Q1, const Quadric& Q2) { return Quadric(Q1.A + Q2.A, Q1.b + Q2.b, Q1.c + Q2.c); }

void quadricErrorSimplify(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geo, double tol) {
  MutationManager mm(mesh, geo);
  quadricErrorSimplify(mesh, geo, tol, mm);
}

void quadricErrorSimplify(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geo, double tol, MutationManager& mm) {
  auto toEigen = [](Vector3 v) -> Eigen::Vector3d {
    Eigen::Vector3d ret;
    ret << v.x, v.y, v.z;
    return ret;
  };
  auto fromEigen = [](Eigen::Vector3d v) -> Vector3 { return Vector3{v(0), v(1), v(2)}; };

  VertexData<Quadric> Q(mesh, Quadric());

  geo.requireFaceNormals();

  for (Face f : mesh.faces()) {
    Eigen::Vector3d n = toEigen(geo.faceNormals[f]);
    Eigen::Matrix3d M = n * n.transpose();
    for (Vertex v : f.adjacentVertices()) {
      Eigen::Vector3d q = toEigen(geo.inputVertexPositions[v]);
      double d = -n.dot(q);

      Q[v] += Quadric(M, d * n, d * d);
    }
  }

  using PotentialEdge = std::tuple<double, Edge>;

  auto cmp = [](const PotentialEdge& a, const PotentialEdge& b) -> bool { return std::get<0>(a) > std::get<0>(b); };

  std::priority_queue<PotentialEdge, std::vector<PotentialEdge>, decltype(cmp)> edgesToCheck(cmp);

  for (Edge e : mesh.edges()) {
    Quadric Qe = Q[e.halfedge().tailVertex()] + Q[e.halfedge().tipVertex()];
    Eigen::Vector3d q = Qe.optimalPoint();
    double cost = Qe.cost(q);
    edgesToCheck.push(std::make_tuple(cost, e));
  }

  while (!edgesToCheck.empty()) {
    PotentialEdge best = edgesToCheck.top();
    edgesToCheck.pop();

    // Stop when collapse becomes too expensive
    double cost = std::get<0>(best);
    if (cost > tol) break;

    Edge e = std::get<1>(best);
    if (!e.isDead()) {

      Vertex v1 = e.halfedge().tailVertex();
      Vertex v2 = e.halfedge().tipVertex();

      // Get edge quadric
      Quadric Qe(Q[v1], Q[v2]);
      Eigen::Vector3d q = Qe.optimalPoint();

      // If either vertex has been collapsed since the edge was pushed
      // onto the queue, the old cost was wrong. In that case, give up
      if (abs(cost - Qe.cost(q)) > 1e-8) continue;

      Vertex v = mm.collapseEdge(e, fromEigen(q));
      if (v == Vertex()) continue;
      Q[v] = Qe;

      for (Edge f : v.adjacentEdges()) {
        Quadric Qf(Q[f.halfedge().tailVertex()], Q[f.halfedge().tipVertex()]);
        Eigen::Vector3d q = Qf.optimalPoint();
        double cost = Qf.cost(q);
        edgesToCheck.push(std::make_tuple(cost, f));
      }
    }
  }

  mesh.compress();
  return;
}
} // namespace surface
} // namespace geometrycentral
