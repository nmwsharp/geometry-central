#include "quadric_error_simplification.h"

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
  MutationManager mm(mesh);
  quadricErrorSimplify(mesh, geo, tol, mm);
}

void quadricErrorSimplify(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geo, double tol, MutationManager& mm) {
  VertexData<QES::Quadric> Q(mesh, QES::Quadric());

  geo.requireFaceNormals();

  for (Face f : mesh.faces()) {
    Eigen::Vector3d n = toEigen(geo.faceNormals[f]);
    Eigen::Matrix3d M = n * n.transpose();
    for (Vertex v : f.adjacentVertices()) {
      Eigen::Vector3d q = toEigen(geo.inputVertexPositions[v]);
      double d = -n.dot(q);

      Q[v] += QES::Quadric(M, d * n, d * d);
    }
  }

  using PotentialEdge = std::tuple<double, Edge>;

  auto cmp = [](const PotentialEdge& a, const PotentialEdge& b) -> bool { return std::get<0>(a) > std::get<0>(b); };

  std::priority_queue<PotentialEdge, std::vector<PotentialEdge>, decltype(cmp)> edgesToCheck(cmp);

  for (Edge e : mesh.edges()) {
    QES::Quadric Qe = Q[e.src()] + Q[e.dst()];
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

      // Get edge quadric
      QES::Quadric Qe(Q[e.src()], Q[e.dst()]);
      Eigen::Vector3d q = Qe.optimalPoint();

      // If either vertex has been collapsed since the edge was pushed
      // onto the queue, the old cost was wrong. In that case, give up
      if (abs(cost - Qe.cost(q)) > 1e-8) continue;

      Vertex v;
      try {
        v = mm.collapseEdge(e);
      } catch (std::runtime_error& e) {
      }

      // v == Vertex() if collapseEdge threw an exception, or if
      // collapseEdge failed and returned Vertex(). In either case, we
      // give up
      if (v == Vertex()) continue;
      geo.vertexPositions[v] = fromEigen(q);
      Q[v] = Qe;

      for (Edge f : v.adjacentEdges()) {
        QES::Quadric Qf(Q[f.src()], Q[f.dst()]);
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
