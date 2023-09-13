#include "geometrycentral/surface/stripe_patterns.h"

#include "geometrycentral/numerical/linear_solvers.h"
#include "geometrycentral/surface/direction_fields.h"

#include <Eigen/Sparse>


using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace geometrycentral {
namespace surface {
namespace {

// Compute the 1-form \omega_{ij} such as defined in eq.7 of [Knoppel et al. 2015]
double computeOmega(IntrinsicGeometryInterface& geometry, const VertexData<Vector2>& directionField,
                    const VertexData<double>& frequencies, const Edge& e, bool* crossesSheets = nullptr) {

  geometry.requireEdgeLengths();
  geometry.requireVertexIndices();
  geometry.requireHalfedgeVectorsInVertex();
  geometry.requireTransportVectorsAlongHalfedge();

  // compute roots of direction field (power representation)
  Vector2 Xi = Vector2::fromAngle(directionField[e.firstVertex()].arg() / 2);
  Vector2 Xj = Vector2::fromAngle(directionField[e.secondVertex()].arg() / 2);

  // check if the directions point the same direction
  Vector2 rij = geometry.transportVectorsAlongHalfedge[e.halfedge()];

  double s = dot(rij * Xi, Xj) > 0 ? 1 : -1;
  if (crossesSheets) *crossesSheets = (s < 0);

  // compute the 1-form value along edge ij
  double lij = geometry.edgeLengths[e];
  double phiI = Xi.arg();
  double phiJ = (s * Xj).arg();

  // get the angle of the edge w.r.t. the endpoints' bases
  double thetaI = geometry.halfedgeVectorsInVertex[e.halfedge()].arg();
  double thetaJ = thetaI + rij.arg();

  double omegaIJ = (lij / 2.) * (frequencies[e.firstVertex()] * cos(phiI - thetaI) +
                                 frequencies[e.secondVertex()] * cos(phiJ - thetaJ));

  return omegaIJ;
}

// Build a Laplace-like matrix with double entries (necessary to represent complex conjugation)
SparseMatrix<double> buildVertexEnergyMatrix(IntrinsicGeometryInterface& geometry,
                                             const VertexData<Vector2>& directionField,
                                             const FaceData<int>& branchIndices,
                                             const VertexData<double>& frequencies) {
  SurfaceMesh& mesh = geometry.mesh;

  geometry.requireVertexIndices();
  geometry.requireHalfedgeCotanWeights();

  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(12 * mesh.nEdges());
  for (Edge e : mesh.edges()) {
    // compute the discrete 1-form
    bool crossesSheet;
    double omegaIJ = computeOmega(geometry, directionField, frequencies, e, &crossesSheet);

    // compute the cotan weight
    double w = 0;
    if (branchIndices[e.halfedge().face()] == 0) {
      w += geometry.halfedgeCotanWeights[e.halfedge()];
    }
    if (!e.isBoundary() && branchIndices[e.halfedge().twin().face()] == 0) {
      w += geometry.halfedgeCotanWeights[e.halfedge().twin()];
    }

    int i = 2 * geometry.vertexIndices[e.halfedge().vertex()];
    int j = 2 * geometry.vertexIndices[e.halfedge().twin().vertex()];

    // add the diagonal terms
    triplets.emplace_back(i + 0, i + 0, w);
    triplets.emplace_back(i + 1, i + 1, w);

    triplets.emplace_back(j + 0, j + 0, w);
    triplets.emplace_back(j + 1, j + 1, w);

    // compute the new transport coefficient
    Vector2 rij = w * Vector2::fromAngle(omegaIJ);

    // these terms are the same in both cases
    triplets.emplace_back(i + 0, j + 0, -rij.x);
    triplets.emplace_back(i + 1, j + 0, rij.y);

    triplets.emplace_back(j + 0, i + 0, -rij.x);
    triplets.emplace_back(j + 0, i + 1, rij.y);

    // if both vectors don't point the same direction, the block represents conjugation as well as multiplication
    if (crossesSheet) rij *= -1;

    triplets.emplace_back(i + 0, j + 1, -rij.y);
    triplets.emplace_back(i + 1, j + 1, -rij.x);

    triplets.emplace_back(j + 1, i + 0, -rij.y);
    triplets.emplace_back(j + 1, i + 1, -rij.x);
  }

  // assemble matrix from triplets
  SparseMatrix<double> vertexEnergyMatrix(2 * mesh.nVertices(), 2 * mesh.nVertices());
  vertexEnergyMatrix.setFromTriplets(triplets.begin(), triplets.end());

  // Shift to avoid singularity
  SparseMatrix<double> eye(2 * mesh.nVertices(), 2 * mesh.nVertices());
  eye.setIdentity();
  vertexEnergyMatrix += 1e-4 * eye;

  return vertexEnergyMatrix;
}

// Build a lumped mass matrix with double entries
SparseMatrix<double> computeRealVertexMassMatrix(IntrinsicGeometryInterface& geometry) {

  SurfaceMesh& mesh = geometry.mesh;

  geometry.requireVertexDualAreas();

  std::vector<Eigen::Triplet<double>> triplets(2 * mesh.nVertices());
  for (size_t i = 0; i < mesh.nVertices(); ++i) {
    double area = geometry.vertexDualAreas[i];
    triplets[2 * i] = Eigen::Triplet<double>(2 * i, 2 * i, area);
    triplets[2 * i + 1] = Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, area);
  }

  // assemble matrix from triplets
  SparseMatrix<double> realVertexMassMatrix(2 * mesh.nVertices(), 2 * mesh.nVertices());
  realVertexMassMatrix.setFromTriplets(triplets.begin(), triplets.end());

  return realVertexMassMatrix;
}

// Solve the generalized eigenvalue problem in equation 9 [Knoppel et al. 2015]
VertexData<Vector2> computeParameterization(IntrinsicGeometryInterface& geometry,
                                            const VertexData<Vector2>& directionField,
                                            const FaceData<int>& branchIndices, const VertexData<double>& frequencies) {

  SurfaceMesh& mesh = geometry.mesh;

  // Compute vertex energy matrix A and mass matrix B
  SparseMatrix<double> energyMatrix = buildVertexEnergyMatrix(geometry, directionField, branchIndices, frequencies);
  SparseMatrix<double> massMatrix = computeRealVertexMassMatrix(geometry);

  // Find the smallest eigenvector
  Vector<double> solution = smallestEigenvectorPositiveDefinite(energyMatrix, massMatrix, 20);

  // Copy the result to a VertexData vector
  VertexData<Vector2> toReturn(mesh);
  for (size_t i = 0; i < mesh.nVertices(); ++i) {
    toReturn[i].x = solution(2 * i);
    toReturn[i].y = solution(2 * i + 1);
    toReturn[i] = toReturn[i].normalize();
  }
  return toReturn;
}

// extract the final texture coordinates from the parameterization
std::tuple<CornerData<double>, FaceData<int>> computeTextureCoordinates(IntrinsicGeometryInterface& geometry,
                                                                        const VertexData<Vector2>& directionField,
                                                                        const VertexData<double>& frequencies,
                                                                        const VertexData<Vector2>& parameterization) {

  SurfaceMesh& mesh = geometry.mesh;

  CornerData<double> textureCoordinates(mesh);
  FaceData<int> paramIndices(mesh);

  geometry.requireTransportVectorsAlongHalfedge();

  for (Face f : mesh.faces()) {
    textureCoordinates[f.halfedge().corner()] = parameterization[f.halfedge().vertex()].arg();

    for (Halfedge he : f.adjacentHalfedges()) {
      // grab the parameter values at vertices
      Vector2 psiI = parameterization[he.vertex()];
      Vector2 psiJ = parameterization[he.next().vertex()];

      // is each halfedge canonical?
      double cIJ = (he.edge().halfedge() != he ? -1 : 1);

      // grab the connection coeffients
      double omegaIJ = cIJ * computeOmega(geometry, directionField, frequencies, he.edge());

      auto crossesSheets = [&](Halfedge hTarget) {
        Halfedge h = f.halfedge();
        Vector2 Xi = Vector2::fromAngle(directionField[h.vertex()].arg() / 2);
        Vector2 Xj = Vector2::fromAngle(directionField[hTarget.vertex()].arg() / 2);
        while (h != hTarget) {
          Vector2 r = geometry.transportVectorsAlongHalfedge[h];
          Xi = r * Xi;
          h = h.next();
        }
        return dot(Xi, Xj) <= 0;
      };

      if (crossesSheets(he)) {
        psiI = psiI.conj();
        omegaIJ *= -cIJ;
      }
      if (crossesSheets(he.next())) {
        psiJ = psiJ.conj();
        omegaIJ *= cIJ;
      }

      // construct complex transport coefficients
      Vector2 rij = Vector2::fromAngle(omegaIJ);

      if (he.next() != f.halfedge()) {
        textureCoordinates[he.next().corner()] = textureCoordinates[he.corner()] + omegaIJ - (rij * psiI / psiJ).arg();
      } else {
        double alpha = textureCoordinates[he.corner()] + omegaIJ - (rij * psiI / psiJ).arg();
        paramIndices[f] = std::round((alpha - textureCoordinates[he.next().corner()]) / (2 * PI));
      }
    }
  }

  return std::tie(textureCoordinates, paramIndices);
}

} // namespace


std::tuple<CornerData<double>, FaceData<int>, FaceData<int>>
computeStripePattern(IntrinsicGeometryInterface& geometry, const VertexData<double>& frequencies,
                     const VertexData<Vector2>& directionField) {
  // find singularities of the direction field
  FaceData<int> branchIndices = computeFaceIndex(geometry, directionField, 2);

  // solve the eigenvalue problem (multiply by 2pi to get the right frequencies)
  VertexData<Vector2> parameterization =
      computeParameterization(geometry, directionField, branchIndices, 2 * PI * frequencies);

  // compute the final corner-based values, along with singularities of the stripe pattern
  CornerData<double> textureCoordinates;
  FaceData<int> zeroIndices;
  std::tie(textureCoordinates, zeroIndices) =
      computeTextureCoordinates(geometry, directionField, 2 * PI * frequencies, parameterization);

  return std::tie(textureCoordinates, zeroIndices, branchIndices);
}

namespace {

// lists all the k's \in \mathbb{Z} such that val1 < 2 k \pi < val2 or val2 < 2 k \pi < val1
std::vector<std::pair<double, int>> crossingsModulo2Pi(double val1, double val2) {
  std::vector<std::pair<double, int>> barys;
  if (val1 == val2) return barys;

  int maxCrossings = std::ceil(std::abs(val1 - val2) / (2 * PI));

  for (int i = 0; i < maxCrossings; ++i) {
    int k = std::ceil(std::min(val1, val2) / (2 * PI)) + i;
    double isoval = 2 * PI * k;

    if (std::max(val1, val2) > isoval) {
      barys.emplace_back((isoval - val2) / (val1 - val2), k);
    }
  }
  return barys;
}

std::tuple<std::vector<Vector3>, std::vector<std::array<int, 2>>, EdgeData<std::vector<int>>>
extractCrossingsFromStripePattern(EmbeddedGeometryInterface& geometry, const CornerData<double>& stripeValues,
                                  const FaceData<int>& stripeIndices, const FaceData<int>& fieldIndices) {
  SurfaceMesh& mesh = geometry.mesh;

  std::vector<Vector3> vertices;
  std::vector<std::array<int, 2>> edges;
  // list of per-edge indices pointing to the vertices list
  EdgeData<std::vector<int>> polylineIndices(mesh);

  for (Face f : mesh.faces()) {
    if (stripeIndices[f] != 0 || fieldIndices[f] != 0) continue; // singularities are ignored in this function

    std::vector<int> faceIsovalues;
    std::vector<int> edgeIndices;
    for (Halfedge he : f.adjacentHalfedges()) {
      // list all the crossings along halfedge he
      std::vector<std::pair<double, int>> isoPoints =
          crossingsModulo2Pi(stripeValues[he.corner()], stripeValues[he.next().corner()]);

      // add them to vertices and isovalues lists
      for (const auto& pair : isoPoints) {
        faceIsovalues.push_back(pair.second);
      }
      // only add new isoline vertices if he.edge() has not yet been visited
      if (polylineIndices[he.edge()].size() == 0) {
        // add polyline vertex indices to the corresponding edge
        for (const auto& pair : isoPoints) {
          double bary = pair.first;
          polylineIndices[he.edge()].push_back(vertices.size());
          vertices.push_back(bary * geometry.vertexPositions[he.tailVertex()] +
                             (1 - bary) * geometry.vertexPositions[he.tipVertex()]);
        }
        edgeIndices.insert(edgeIndices.end(), polylineIndices[he.edge()].begin(), polylineIndices[he.edge()].end());
      }
      // otherwise reuse the existing polylineIndices to build the edgeIndices
      else {
        // need to make sure the polylineIndices are enumerated in the right order, otherwise we flip them
        if (stripeValues[he.corner()] < stripeValues[he.next().corner()] &&
                stripeValues[he.twin().corner()] < stripeValues[he.twin().next().corner()] ||
            stripeValues[he.corner()] > stripeValues[he.next().corner()] &&
                stripeValues[he.twin().corner()] > stripeValues[he.twin().next().corner()])
          edgeIndices.insert(edgeIndices.end(), polylineIndices[he.edge()].rbegin(), polylineIndices[he.edge()].rend());
        else
          edgeIndices.insert(edgeIndices.end(), polylineIndices[he.edge()].begin(), polylineIndices[he.edge()].end());
      }
    }

    // build up list of edges
    std::vector<bool> visited(faceIsovalues.size(), false);
    std::vector<std::array<int, 2>> faceEdges;
    // match vertices by their isovalues (if point hasn't been visited yet, find the one which has the same isovalue)
    for (int i = 0; i < faceIsovalues.size(); ++i) {
      if (visited[i]) continue;

      visited[i] = true;
      // find the (unique) point with index j which has the same isovalue as point i
      for (int j = i; j < faceIsovalues.size(); ++j) {
        if (visited[j]) continue;

        if (faceIsovalues[j] == faceIsovalues[i]) {
          visited[j] = true;
          faceEdges.push_back({edgeIndices[i], edgeIndices[j]});
        }
      }
    }
    edges.insert(edges.end(), faceEdges.begin(), faceEdges.end());
  }
  return std::make_tuple(vertices, edges, polylineIndices);
}


/**
 * This part tries to connect the isolines on singularities
 * tries to connect isolines that are (a) close together and
 * (b) whose crossing is as aligned with the input direction field as possible
 */
void connectIsolinesOnSingularities(EmbeddedGeometryInterface& geometry, const CornerData<double>& stripeValues,
                                    const FaceData<int>& stripeIndices, const FaceData<int>& branchIndices,
                                    EdgeData<std::vector<int>>& polylineIndices, const FaceData<Vector3>& vectorField,
                                    std::vector<Vector3>& points, std::vector<std::array<int, 2>>& edges) {
  SurfaceMesh& mesh = geometry.mesh;
  for (Face f : mesh.faces()) {
    if (branchIndices[f] != 0) {
      // do nothing for now
    } else if (stripeIndices[f] != 0) {
      int pointsCount = 0;
      std::vector<std::pair<std::vector<double>, Halfedge>> crossings;
      for (Halfedge he : f.adjacentHalfedges()) {
        std::vector<std::pair<double, int>> isoPoints;
        if (he.next() == f.halfedge())
          isoPoints = crossingsModulo2Pi(stripeValues[he.corner()],
                                         stripeValues[he.next().corner()] + 2 * stripeIndices[f] * PI);
        else
          isoPoints = crossingsModulo2Pi(stripeValues[he.corner()], stripeValues[he.next().corner()]);
        pointsCount += isoPoints.size();

        std::vector<double> barys;
        for (const auto& pair : isoPoints) {
          barys.push_back(pair.first);
        }
        crossings.emplace_back(barys, he);
      }

      std::array<std::vector<Vector3>, 3> facePoints;
      std::array<std::vector<int>, 3> faceIndices;
      int idx = 0;
      bool ignoreFace = false;
      std::sort(crossings.begin(), crossings.end(),
                [](std::pair<std::vector<double>, Halfedge> a, std::pair<std::vector<double>, Halfedge> b) {
                  return a.first.size() > b.first.size();
                });

      for (auto crossing : crossings) {
        if (crossing.first.size() == 0) continue;
        Halfedge he = crossing.second;

        for (double bary : crossing.first) {
          facePoints[idx].push_back(bary * geometry.vertexPositions[he.tailVertex()] +
                                    (1 - bary) * geometry.vertexPositions[he.tipVertex()]);
        }

        if (polylineIndices[he.edge()].size() > 0) { // use the existing indices
          // need to make sure the polylineIndices are enumerated in the right order, otherwise we flip them
          if (norm(facePoints[idx][0] - points[polylineIndices[he.edge()][0]]) > 1e-6) {
            faceIndices[idx].insert(faceIndices[idx].end(), polylineIndices[he.edge()].rbegin(),
                                    polylineIndices[he.edge()].rend());
          } else {
            faceIndices[idx] = polylineIndices[he.edge()];
          }
        } else { // add the intersections to the points lists and initialize the corresponding list of indices
          for (double bary : crossing.first) {
            points.push_back(bary * geometry.vertexPositions[he.tailVertex()] +
                             (1 - bary) * geometry.vertexPositions[he.tipVertex()]);
            polylineIndices[he.edge()].push_back(points.size() - 1);
          }
          faceIndices[idx] = polylineIndices[he.edge()];
        }
        ++idx;
      }

      if (ignoreFace) continue;

      std::vector<std::vector<bool>> matched(3);
      for (int i = 0; i < 3; ++i) matched[i] = std::vector<bool>(facePoints[i].size(), false);

      // match points v and w which minimize abs(dot(v - w, vectorField[f]))
      while (pointsCount > std::abs(stripeIndices[f])) {
        std::array<int, 4> minIdx;
        double minVal = std::numeric_limits<double>::max();

        // if triangular inequality doesn't hold, don't examine the second edge
        int nI;
        if (faceIndices[0].size() > faceIndices[1].size() + faceIndices[2].size())
          nI = 1;
        else
          nI = 2;

        for (int i = 0; i < nI; ++i) {
          for (int k = 0; k < facePoints[i].size(); ++k) {
            if (matched[i][k]) continue;

            for (int j = i + 1; j < 3; ++j) {
              for (int l = 0; l < facePoints[j].size(); ++l) {
                if (matched[j][l]) continue;

                Vector3 v = facePoints[i][k];
                Vector3 w = facePoints[j][l];
                if (abs(dot(v - w, vectorField[f])) < minVal) {
                  minVal = abs(dot(v - w, vectorField[f]));
                  minIdx = {i, k, j, l};
                }
              }
            }
          }
        }

        matched[minIdx[0]][minIdx[1]] = true;
        matched[minIdx[2]][minIdx[3]] = true;
        edges.push_back({faceIndices[minIdx[0]][minIdx[1]], faceIndices[minIdx[2]][minIdx[3]]});
        pointsCount -= 2;
      }
    }
  }
}

} // namespace

std::tuple<std::vector<Vector3>, std::vector<std::array<int, 2>>>
extractPolylinesFromStripePattern(EmbeddedGeometryInterface& geometry, const CornerData<double>& values,
                                  const FaceData<int>& stripeIndices, const FaceData<int>& fieldIndices,
                                  const VertexData<Vector2>& directionField, bool connectOnSingularities) {
  SurfaceMesh& mesh = geometry.mesh;

  std::vector<Vector3> points;
  std::vector<std::array<int, 2>> edges;
  EdgeData<std::vector<int>> indicesPerEdge;
  std::tie(points, edges, indicesPerEdge) =
      extractCrossingsFromStripePattern(geometry, values, stripeIndices, fieldIndices);

  if (connectOnSingularities) {
    // convert vertex-based direction field to a face-based representation
    geometry.requireVertexTangentBasis();

    FaceData<Vector3> vectorField(mesh);
    for (Face f : mesh.faces()) {
      Vector3 avgDir = {0, 0, 0};
      for (Vertex v : f.adjacentVertices()) {
        Vector2 vector = directionField[v].pow(0.5);
        avgDir +=
            (geometry.vertexTangentBasis[v][0] * vector.x + geometry.vertexTangentBasis[v][1] * vector.y).normalize();
      }

      vectorField[f] = avgDir.normalize();
    }

    connectIsolinesOnSingularities(geometry, values, stripeIndices, fieldIndices, indicesPerEdge, vectorField, points,
                                   edges);
  }

  return std::make_tuple(points, edges);
}


std::tuple<std::vector<Vector3>, std::vector<std::array<int, 2>>>
computeStripePatternPolylines(EmbeddedGeometryInterface& geometry, const VertexData<double>& frequencies,
                              const VertexData<Vector2>& directionField, bool connectIsolines) {
  CornerData<double> stripeValues;
  FaceData<int> stripeIndices, fieldIndices;
  std::tie(stripeValues, stripeIndices, fieldIndices) = computeStripePattern(geometry, frequencies, directionField);

  return extractPolylinesFromStripePattern(geometry, stripeValues, stripeIndices, fieldIndices, directionField,
                                           connectIsolines);
}

} // namespace surface
} // namespace geometrycentral