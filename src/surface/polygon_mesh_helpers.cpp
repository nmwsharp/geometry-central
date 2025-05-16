#include "geometrycentral/surface/polygon_mesh_helpers.h"

namespace geometrycentral {
namespace surface {

// Helper functions for Bunge et al. "Polygon Laplacian Made Simple" (2020).
// Copyright (C) 2020 Astrid Bunge, Philipp Herholz, Misha Kazhdan, Mario Botsch, MIT license
// (Modified to work in geometry-central. Original code can be found here:
// https://github.com/mbotsch/polygon-laplacian)

Eigen::VectorXd simplePolygonVirtualVertex(const Eigen::MatrixXd& poly) {

  // Given a polygon face, computes the affine weights that determine the position of the virtual vertex that minimizes
  // the sum of the squared areas of the triangles in the induced triangle fan. While the location of this vertex (the
  // minimizer) is unique, its expression as an affine combination of the polygon verties may not be -- regularize by
  // picking the weights with minimum L_2 norm, which encourages the weights to be as uniform as possible.

  int n = poly.rows();
  Eigen::VectorXd weights(n);
  Eigen::MatrixXd J(n, n);
  Eigen::VectorXd b(n);
  for (int i = 0; i < n; i++) {
    Eigen::Vector3d pk = poly.row(i);

    double Bk1_d2 = 0.0;
    double Bk1_d1 = 0.0;

    double Bk2_d0 = 0.0;
    double Bk2_d2 = 0.0;

    double Bk3_d0 = 0.0;
    double Bk3_d1 = 0.0;

    double CBk = 0.0;
    Eigen::Vector3d d = Eigen::MatrixXd::Zero(3, 1);

    for (int j = 0; j < n; j++) {
      Eigen::Vector3d pi = poly.row(j);
      Eigen::Vector3d pj = poly.row((j + 1) % n);
      d = pi - pj;

      double Bik1 = d(1) * pk(2) - d(2) * pk(1);
      double Bik2 = d(2) * pk(0) - d(0) * pk(2);
      double Bik3 = d(0) * pk(1) - d(1) * pk(0);

      double Ci1 = d(1) * pi(2) - d(2) * pi(1);
      double Ci2 = d(2) * pi(0) - d(0) * pi(2);
      double Ci3 = d(0) * pi(1) - d(1) * pi(0);

      Bk1_d1 += d(1) * Bik1;
      Bk1_d2 += d(2) * Bik1;

      Bk2_d0 += d(0) * Bik2;
      Bk2_d2 += d(2) * Bik2;

      Bk3_d0 += d(0) * Bik3;
      Bk3_d1 += d(1) * Bik3;

      CBk += Ci1 * Bik1 + Ci2 * Bik2 + Ci3 * Bik3;
    }
    for (int k = 0; k < n; k++) {
      Eigen::Vector3d xj = poly.row(k);
      J(i, k) =
          0.5 * (xj(2) * Bk1_d1 - xj(1) * Bk1_d2 + xj(0) * Bk2_d2 - xj(2) * Bk2_d0 + xj(1) * Bk3_d0 - xj(0) * Bk3_d1);
    }
    b(i) = 0.5 * CBk;
  }

  Eigen::MatrixXd M(n + 1, n);
  M.block(0, 0, n, n) = 4 * J;
  M.block(n, 0, 1, n).setOnes();

  Eigen::VectorXd b_(n + 1);
  b_.block(0, 0, n, 1) = 4 * b;

  b_(n) = 1.;
  weights = M.completeOrthogonalDecomposition().solve(b_).topRows(n);

  return weights;
}

// Helper functions for de Goes et al. "Discrete Differential Operators on Polygonal Meshes" (2020).
// Use of this source code is governed by a LGPL-3.0 license.
// (Modified to work in geometry-central. Original code can be found here:
// https://github.com/DGtal-team/DGtal/blob/master/src/DGtal/dec/PolygonalCalculus.h)

Eigen::MatrixXd polygonAveragingMatrix(const Face& f) {
  size_t d = f.degree();
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(d, d);
  for (size_t i = 0; i < d; i++) {
    A(i, (i + 1) % d) = 0.5;
    A(i, i) = 0.5;
  }
  return A;
}

Eigen::MatrixXd polygonDerivativeMatrix(const Face& f) {
  size_t d = f.degree();
  Eigen::MatrixXd D = Eigen::MatrixXd::Zero(d, d);
  for (size_t i = 0; i < d; i++) {
    D(i, (i + 1) % d) = 1.;
    D(i, i) = -1.;
  }
  return D;
}

Eigen::Matrix3d bracket(const Eigen::Vector3d& n) {
  Eigen::Matrix3d B;
  B << 0., -n[2], n[1], n[2], 0., -n[0], -n[1], n[0], 0.;
  return B;
}

Eigen::Vector3d project(const Eigen::Vector3d& u, const Eigen::Vector3d& n) {
  // Project u on the orthgonal of n.
  return u - (u.dot(n) / n.squaredNorm()) * n;
}

Eigen::MatrixXd kroneckerWithI2(const Eigen::MatrixXd& M) {
  size_t h = M.rows();
  size_t w = M.cols();
  Eigen::MatrixXd MK = Eigen::MatrixXd::Zero(h * 2, w * 2);
  for (size_t j = 0; j < h; j++)
    for (size_t i = 0; i < w; i++) {
      MK(2 * j, 2 * i) = M(j, i);
      MK(2 * j + 1, 2 * i + 1) = M(j, i);
    }
  return MK;
}

} // namespace surface
} // namespace geometrycentral