#include "geometrycentral/surface/polygon_mesh_helpers.h"

namespace geometrycentral {
namespace surface {

// Helper functions for Bunge et al. "Polygon Laplacian Made Simple" (2020).
// Copyright (C) 2020 Astrid Bunge, Philipp Herholz, Misha Kazhdan, Mario Botsch, MIT license
// (Modified to work in geometry-central. Original code can be found here:
// https://github.com/mbotsch/polygon-laplacian)

Eigen::MatrixXd simplePolygonMassMatrix(EmbeddedGeometryInterface& geometry,
                                        const FaceData<Eigen::VectorXd>& virtualRefinementAreaWeights, const Face& f) {

  size_t n = f.degree();
  Eigen::MatrixXd poly = polygonPositionMatrix(geometry, f);
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(n, n);
  const Eigen::VectorXd& weights = virtualRefinementAreaWeights[f];
  Eigen::Vector3d virtualVertex = poly.transpose() * weights;
  Eigen::VectorXd ln = Eigen::VectorXd::Zero(n + 1);
  double l[3], l2[3]; // lengths, lengths squared
  // Build triangle fan mass and cotan matrices
  for (size_t i = 0; i < n; i++) {
    const size_t i1 = (i + 1) % n;
    l2[2] = (poly.row(i) - poly.row(i1)).squaredNorm();
    l2[0] = (poly.row(i1) - virtualVertex.transpose()).squaredNorm();
    l2[1] = (poly.row(i) - virtualVertex.transpose()).squaredNorm();
    l[0] = std::sqrt(l2[0]);
    l[1] = std::sqrt(l2[1]);
    l[2] = std::sqrt(l2[2]);
    const double arg =
        (l[0] + (l[1] + l[2])) * (l[2] - (l[0] - l[1])) * (l[2] + (l[0] - l[1])) * (l[0] + (l[1] - l[2]));
    const double area = 0.25 * std::sqrt(arg);
    l[0] = 1.0 / 6.0 * area;
    l[1] = 1.0 / 12.0 * area;
    M(i1, i1) += 1.0 / 6.0 * area;
    M(i, i) += 1.0 / 6.0 * area;
    M(i1, i) += 1.0 / 12.0 * area;
    M(i, i1) += 1.0 / 12.0 * area;
    ln(i1) += l[1];
    ln(i) += l[1];
    ln(n) += l[0];
  }
  // Apply prolongation
  for (size_t j = 0; j < n; ++j)
    for (size_t i = 0; i < n; ++i) M(i, j) += weights(i) * ln(j) + weights(j) * ln(i) + weights(i) * weights(j) * ln(n);

  return M;
}

Eigen::MatrixXd simplePolygonStiffnessMatrix(EmbeddedGeometryInterface& geometry,
                                             const FaceData<Eigen::VectorXd>& virtualRefinementAreaWeights,
                                             const Face& f) {

  size_t n = f.degree();
  Eigen::MatrixXd poly = polygonPositionMatrix(geometry, f);
  Eigen::MatrixXd S = Eigen::MatrixXd::Zero(n, n);
  const Eigen::VectorXd& weights = virtualRefinementAreaWeights[f];
  Eigen::Vector3d virtualVertex = poly.transpose() * weights;
  Eigen::VectorXd ln = Eigen::VectorXd::Zero(n + 1);
  double l[3], l2[3]; // lengths, lengths squared
  // Build triangle fan mass and cotan matrices
  for (size_t i = 0; i < n; i++) {
    const size_t i1 = (i + 1) % n;
    l2[2] = (poly.row(i) - poly.row(i1)).squaredNorm();
    l2[0] = (poly.row(i1) - virtualVertex.transpose()).squaredNorm();
    l2[1] = (poly.row(i) - virtualVertex.transpose()).squaredNorm();
    l[0] = std::sqrt(l2[0]);
    l[1] = std::sqrt(l2[1]);
    l[2] = std::sqrt(l2[2]);
    const double arg =
        (l[0] + (l[1] + l[2])) * (l[2] - (l[0] - l[1])) * (l[2] + (l[0] - l[1])) * (l[0] + (l[1] - l[2]));
    const double area = 0.5 * std::sqrt(arg);
    if (area > 1e-7) {
      l[0] = 0.25 * (l2[1] + l2[2] - l2[0]) / area;
      l[1] = 0.25 * (l2[2] + l2[0] - l2[1]) / area;
      l[2] = 0.25 * (l2[0] + l2[1] - l2[2]) / area;

      S(i1, i1) += l[0];
      S(i, i) += l[1];
      S(i1, i) -= l[2];
      S(i, i1) -= l[2];
      S(i, i) += l[2];
      S(i1, i1) += l[2];

      ln(i1) -= l[0];
      ln(i) -= l[1];
      ln(n) += l[0] + l[1];
    }
  }
  // Apply prolongation
  for (size_t j = 0; j < n; ++j)
    for (size_t i = 0; i < n; ++i) S(i, j) += weights(i) * ln(j) + weights(j) * ln(i) + weights(i) * weights(j) * ln(n);

  return S;
}

Eigen::MatrixXd polygonPositionMatrix(EmbeddedGeometryInterface& geometry, const Face& f) {

  Eigen::MatrixXd poly(f.degree(), 3);
  int i = 0;
  for (Vertex v : f.adjacentVertices()) {
    for (int j = 0; j < 3; j++) {
      poly(i, j) = geometry.vertexPositions[v][j];
    }
    i++;
  }
  return poly;
}

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

Eigen::MatrixXd polygonPerFaceLaplacian(EmbeddedGeometryInterface& geometry, const Face& f,
                                        const double polygonLambda) {
  Eigen::MatrixXd Df = polygonDerivativeMatrix(f);
  return Df.transpose() * polygonPerFaceInnerProductMatrix(geometry, f, polygonLambda) * Df; // build positive-definite
}

Eigen::MatrixXd polygonPerFaceInnerProductMatrix(EmbeddedGeometryInterface& geometry, const Face& f,
                                                 const double polygonLambda) {
  // faceAreas is required in computePolygonLaplacian(), which loops over all faces and calls
  // polygonPerFaceLaplacian(), which calls polygonPerFaceInnerProductMatrix().
  Eigen::MatrixXd Uf = polygonSharp(geometry, f);
  Eigen::MatrixXd Pf = polygonProjectionMatrix(geometry, f);
  double A = geometry.faceAreas[f];
  return A * Uf.transpose() * Uf + polygonLambda * Pf.transpose() * Pf;
}

Eigen::MatrixXd polygonPerFaceConnectionLaplacian(EmbeddedGeometryInterface& geometry, const Face& f,
                                                  const double polygonLambda) {
  // faceAreas is required in computePolygonVertexConnectionLaplacian(), which loops over all
  // faces and calls polygonPerFaceConnectionLaplacian().
  Eigen::MatrixXd G = polygonCovariantGradient(geometry, f);
  Eigen::MatrixXd P = polygonCovariantProjection(geometry, f);
  double A = geometry.faceAreas[f];
  Eigen::MatrixXd L = A * G.transpose() * G + polygonLambda * P.transpose() * P;
  return L;
}

Eigen::MatrixXd polygonBlockConnection(EmbeddedGeometryInterface& geometry, const Face& f) {
  size_t d = f.degree();
  Eigen::MatrixXd R = Eigen::MatrixXd::Zero(2 * d, 2 * d);
  size_t cpt = 0;
  for (Vertex v : f.adjacentVertices()) {
    Eigen::Matrix2d Rv = Rvf(geometry, v, f);
    R.block<2, 2>(2 * cpt, 2 * cpt) = Rv;
    cpt++;
  }
  return R;
}

Eigen::MatrixXd polygonCovariantGradient(EmbeddedGeometryInterface& geometry, const Face& f) {
  return kroneckerWithI2(Tf(geometry, f).transpose() * polygonPerFaceGradientMatrix(geometry, f)) *
         polygonBlockConnection(geometry, f);
}

Eigen::MatrixXd polygonCovariantProjection(EmbeddedGeometryInterface& geometry, const Face& f) {
  return kroneckerWithI2(polygonProjectionMatrix(geometry, f) * polygonDerivativeMatrix(f)) *
         polygonBlockConnection(geometry, f);
}

Eigen::MatrixXd polygonProjectionMatrix(EmbeddedGeometryInterface& geometry, const Face& f) {
  size_t d = f.degree();
  Eigen::MatrixXd P = Eigen::MatrixXd::Identity(d, d) - polygonFlat(geometry, f) * polygonSharp(geometry, f);
  return P;
}

Eigen::MatrixXd polygonPerFaceGradientMatrix(EmbeddedGeometryInterface& geometry, const Face& f) {
  // faceNormals, faceAreas are required in computePolygonDECOperators(), etc. which
  // is ultimately where this function is used (repeatedly in a loop over all faces).

  // equivalent to applying the (per-face) sharp operator to the (per-face) exterior derivative
  double A = geometry.faceAreas[f];
  Vector3 n = geometry.faceNormals[f];
  Eigen::Vector3d N = {n[0], n[1], n[2]};
  return 1. / A * bracket(N) * polygonCoGradientMatrix(geometry, f);
}

Eigen::MatrixXd polygonCoGradientMatrix(EmbeddedGeometryInterface& geometry, const Face& f) {
  return polygonEdgeVectorMatrix(geometry, f).transpose() * polygonAveragingMatrix(f);
}

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

Eigen::MatrixXd polygonEdgeVectorMatrix(EmbeddedGeometryInterface& geometry, const Face& f) {
  return polygonDerivativeMatrix(f) * polygonPositionMatrix(geometry, f);
}

Eigen::MatrixXd polygonEdgeMidpointMatrix(EmbeddedGeometryInterface& geometry, const Face& f) {
  return polygonAveragingMatrix(f) * polygonPositionMatrix(geometry, f);
}

Eigen::MatrixXd polygonFlat(EmbeddedGeometryInterface& geometry, const Face& f) {
  // faceNormals are required in computePolygonDECOperators(), which
  // is ultimately where this function is used (repeatedly in a loop over all faces).
  Vector3 n = geometry.faceNormals[f];
  Eigen::Vector3d N = {n[0], n[1], n[2]};
  return polygonEdgeVectorMatrix(geometry, f) * (Eigen::MatrixXd::Identity(3, 3) - N * N.transpose());
}

Eigen::MatrixXd polygonSharp(EmbeddedGeometryInterface& geometry, const Face& f) {
  // faceAreas, faceNormals are required in computePolygonDECOperators(), which
  // is ultimately where this function is used (repeatedly in a loop over all faces).
  size_t d = f.degree();
  double A = geometry.faceAreas[f];
  Vector3 n = geometry.faceNormals[f];
  Eigen::Vector3d N = {n[0], n[1], n[2]};
  Eigen::Vector3d c = polygonCentroid(geometry, f);
  return 1. / A * bracket(N) *
         (polygonEdgeMidpointMatrix(geometry, f).transpose() - c * Eigen::VectorXd::Ones(d).transpose());
}

Eigen::Vector3d polygonCentroid(EmbeddedGeometryInterface& geometry, const Face& f) {
  // vertexPositions are required in computePolygonVertexConnectionLaplacian(), which
  // is ultimately where this function is used (repeatedly in a loop over all faces).
  Vector3 c = {0, 0, 0};
  for (Vertex v : f.adjacentVertices()) {
    c += geometry.vertexPositions[v];
  }
  c /= f.degree();
  return Eigen::Vector3d(c[0], c[1], c[2]);
}

Eigen::MatrixXd Tv(EmbeddedGeometryInterface& geometry, const Vertex& v) {
  // vertexTangentBasis is required in computePolygonVertexConnectionLaplacian(), which is
  // ultimately where this function is used (repeatedly in a loop over all faces).

  // Return 3 x 2 matrix defining the tangent space at vertex v, with basis vectors in columns.
  Vector3 xVec = geometry.vertexTangentBasis[v][0];
  Vector3 yVec = geometry.vertexTangentBasis[v][1];
  Eigen::Vector3d uu = {xVec[0], xVec[1], xVec[2]};
  Eigen::Vector3d vv = {yVec[0], yVec[1], yVec[2]};
  Eigen::MatrixXd B(3, 2);
  B.col(0) = uu;
  B.col(1) = vv;
  return B;
}

Eigen::MatrixXd Tf(EmbeddedGeometryInterface& geometry, const Face& f) {
  // faceTangentBasis is required in computePolygonVertexConnectionLaplacian(), which is
  // ultimately where this function is used (repeatedly in a loop over all faces).

  // Return 3 x 2 matrix defining the tangent space at face f, with basis vectors in columns.
  Vector3 xVec = geometry.faceTangentBasis[f][0];
  Vector3 yVec = geometry.faceTangentBasis[f][1];
  Eigen::Vector3d uu = {xVec[0], xVec[1], xVec[2]};
  Eigen::Vector3d vv = {yVec[0], yVec[1], yVec[2]};
  Eigen::MatrixXd B(3, 2);
  B.col(0) = uu;
  B.col(1) = vv;
  return B;
}

Eigen::Matrix2d Rvf(EmbeddedGeometryInterface& geometry, const Vertex& v, const Face& f) {
  return Tf(geometry, f).transpose() * Qvf(geometry, v, f) * Tv(geometry, v);
}

Eigen::Matrix3d Qvf(EmbeddedGeometryInterface& geometry, const Vertex& v, const Face& f) {
  // faceNormals, polygonVertexNormals are required in
  // computePolygonVertexConnectionLaplacian(), which is ultimately where this function is
  // used (repeatedly in a loop over all faces).

  // Return 3 x 3 rotation matrix to align n_v to n_f.
  Vector3 n = geometry.faceNormals[f];
  Eigen::Vector3d nf = {n[0], n[1], n[2]};
  Eigen::Vector3d nv = geometry.polygonVertexNormals[v];
  double c = nv.dot(nf);

  // Special case for opposite nv and nf vectors.
  if (std::abs(c + 1.0) < 1e-5) return -Eigen::Matrix3d::Identity();

  Eigen::Vector3d vv = nv.cross(nf);
  Eigen::Matrix3d skew = bracket(vv);
  return Eigen::Matrix3d::Identity() + skew + 1.0 / (1.0 + c) * skew * skew;
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