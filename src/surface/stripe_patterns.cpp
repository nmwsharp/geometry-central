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
std::vector<double> crossingsModulo2Pi(double val1, double val2) {
  std::vector<double> barys;
  if (val1 == val2) return barys;

  int maxCrossings = std::ceil(std::abs(val1 - val2) / (2 * PI));

  if (val1 < val2) {
    for (int i = 0; i < maxCrossings; ++i) {
      int k = std::ceil(val1 / (2 * PI)) + i;
      double isoval = 2 * PI * k;

      if (isoval < val2) {
        barys.push_back((isoval - val2) / (val1 - val2));
      }
    }
  } else {
    for (int i = maxCrossings - 1; i >= 0; --i) {
      int k = std::ceil(val2 / (2 * PI)) + i;
      double isoval = 2 * PI * k;

      if (isoval < val1) {
        barys.push_back((isoval - val2) / (val1 - val2));
      }
    }
  }

  return barys;
}

// Matches crossings based on a strategy proposed in "Navigating intrinsic triangulations" [Sharp et al. 2019].
// See https://github.com/nmwsharp/geometry-central/pull/89#issuecomment-936150222 for more details
std::vector<std::array<size_t, 2>> matchCrossings(const std::vector<std::vector<size_t>>& crossings) {
  assert(crossings.size() == 3);

  int idxIJ = 2;
  if (crossings[0].size() >= crossings[1].size() && crossings[0].size() >= crossings[2].size()) {
    idxIJ = 0;
  } else if (crossings[1].size() >= crossings[2].size() && crossings[1].size() >= crossings[0].size()) {
    idxIJ = 1;
  }
  int idxJK = (idxIJ + 1) % 3;
  int idxKI = (idxIJ + 2) % 3;

  const std::vector<size_t>& IJ = crossings[idxIJ];
  const std::vector<size_t>& JK = crossings[idxJK];
  const std::vector<size_t>& KI = crossings[idxKI];

  int nIJ = IJ.size();
  int nJK = JK.size();
  int nKI = KI.size();

  assert(nIJ >= nJK && nIJ >= nKI);
  assert(nIJ <= nJK + nKI);
  assert((nIJ + nJK + nKI) % 2 == 0);

  std::vector<std::array<size_t, 2>> matchings;
  if (nIJ == nJK + nKI) { // Case 1: all edges intersecting ijk cross a common edge ij
    // match IJ with IK
    for (int m = 0; m < nKI; ++m) {
      matchings.push_back({IJ[m], KI[nKI - m - 1]});
    }
    // match IJ with KJ
    for (int m = 0; m < nJK; ++m) {
      matchings.push_back({IJ[nKI + m], JK[nJK - m - 1]});
    }
  } else { // Case 2: there is no common edge
    int nRemainingCrossings = (nIJ + nJK + nKI) / 2;
    int m = 0;
    while (nRemainingCrossings > nJK) {
      matchings.push_back({IJ[m], KI[nKI - m - 1]});
      ++m;
      nRemainingCrossings -= 1;
    }

    int l = 0;
    while (nRemainingCrossings > nKI - m) {
      matchings.push_back({IJ[nIJ - 1 - l], JK[l]});
      nRemainingCrossings -= 1;
      ++l;
    }

    int p = 0;
    while (nRemainingCrossings > 0) {
      matchings.push_back({JK[nJK - 1 - p], KI[p]});
      ++p;
      nRemainingCrossings -= 1;
    }
  }

  return matchings;
}

std::tuple<std::vector<Vector3>, std::vector<std::array<size_t, 2>>, EdgeData<std::vector<size_t>>>
extractCrossingsFromStripePattern(EmbeddedGeometryInterface& geometry, const CornerData<double>& stripeValues,
                                  const FaceData<int>& stripeIndices, const FaceData<int>& fieldIndices) {
  SurfaceMesh& mesh = geometry.mesh;

  std::vector<Vector3> vertices;
  std::vector<std::array<size_t, 2>> edges;
  // list of per-edge indices pointing to the vertices list
  EdgeData<std::vector<size_t>> polylineIndices(mesh);

  for (Face f : mesh.faces()) {
    if (stripeIndices[f] != 0 || fieldIndices[f] != 0) continue; // singularities are ignored in this function

    std::vector<std::vector<size_t>> edgeIndices(3);
    int idx = 0;
    for (Halfedge he : f.adjacentHalfedges()) {
      // list all the crossings along halfedge he
      std::vector<double> isoPoints = crossingsModulo2Pi(stripeValues[he.corner()], stripeValues[he.next().corner()]);

      // only add new isoline vertices if he.edge() has not yet been visited
      if (polylineIndices[he.edge()].size() == 0) {
        // add polyline vertex indices to the corresponding edge
        for (double bary : isoPoints) {
          polylineIndices[he.edge()].push_back(vertices.size());
          vertices.push_back(bary * geometry.vertexPositions[he.tailVertex()] +
                             (1 - bary) * geometry.vertexPositions[he.tipVertex()]);
        }
        edgeIndices[idx] = polylineIndices[he.edge()];
      } else {
        edgeIndices[idx].insert(edgeIndices[idx].end(), polylineIndices[he.edge()].rbegin(),
                                polylineIndices[he.edge()].rend());
      }
      ++idx;
    }

    std::vector<std::array<size_t, 2>> matchings = matchCrossings(edgeIndices);
    edges.insert(edges.end(), matchings.begin(), matchings.end());
  }
  return std::make_tuple(vertices, edges, polylineIndices);
}

// Returns positions of maximum values, skipping some values to focus on IJ if the triangular inequality needs to be
// preserved. M = nIJ - nJK - nKI, should be at most 0
std::vector<std::vector<size_t>> removeCrossings(std::vector<std::pair<double, std::array<size_t, 2>>>& values, int N,
                                                 std::vector<std::vector<size_t>>& faceIndices) {
  using elemType = std::pair<double, std::array<size_t, 2>>;
  std::sort(values.begin(), values.end(), [](elemType& a, elemType& b) { return a.first > b.first; });

  std::array<std::vector<bool>, 3> removed;
  for (int i = 0; i < 3; ++i) removed[i] = std::vector<bool>(faceIndices[i].size(), false);

  assert(faceIndices.size() == 3);
  size_t idxIJ = 0;
  if (faceIndices[1].size() >= faceIndices[2].size() && faceIndices[1].size() >= faceIndices[0].size()) {
    idxIJ = 1;
  } else if (faceIndices[2].size() >= faceIndices[1].size() && faceIndices[2].size() >= faceIndices[0].size()) {
    idxIJ = 2;
  }

  int M = faceIndices[idxIJ].size() - faceIndices[(idxIJ + 1) % 3].size() - faceIndices[(idxIJ + 2) % 3].size();
  int nbRemoved = 0; // total number of removed points
  int i = 0;
  while (nbRemoved < N) {
    std::array<size_t, 2> indices = values[i].second;
    if (indices[0] == idxIJ) {
      --M;
    } else {
      if (N - nbRemoved <= M) {
        ++i;
        continue;
      }
      ++M;
    }

    removed[indices[0]][indices[1]] = true;
    ++nbRemoved;
    ++i;
  }
  assert(M <= 0);

  std::vector<std::vector<size_t>> crossingIndices(3);
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < faceIndices[i].size(); ++j)
      if (!removed[i][j]) crossingIndices[i].push_back(faceIndices[i][j]);

  return crossingIndices;
}

// finds the minimum of abs(dot(v - w, vector)) for all points w except those on ignoredEdge
double findMin(const std::vector<std::vector<Vector3>>& points, const Vector3& v, const Vector3& vector,
               size_t ignoredEdge) {
  double min = std::numeric_limits<double>::max();
  for (size_t i = 0; i < points.size(); ++i) {
    if (i == ignoredEdge) continue;
    for (size_t j = 0; j < points[i].size(); ++j) {
      Vector3 w = points[i][j];
      if (abs(dot(v - w, vector)) < min) {
        min = abs(dot(v - w, vector));
      }
    }
  }
  return min;
}

/**
 * This part tries to connect the isolines on singularities
 * Tries to connect isolines that are (a) close together and
 * (b) whose crossing is as aligned with the input direction field as possible
 */
void connectIsolinesOnSingularities(EmbeddedGeometryInterface& geometry, const CornerData<double>& stripeValues,
                                    const FaceData<int>& stripeIndices, const FaceData<int>& branchIndices,
                                    EdgeData<std::vector<size_t>>& polylineIndices,
                                    const FaceData<Vector3>& vectorField, std::vector<Vector3>& vertices,
                                    std::vector<std::array<size_t, 2>>& edges) {
  SurfaceMesh& mesh = geometry.mesh;
  geometry.requireFaceIndices();

  for (Face f : mesh.faces()) {
    if (stripeIndices[f] != 0) {
      std::vector<std::vector<Vector3>> facePoints;
      std::vector<std::vector<size_t>> faceIndices;
      for (Halfedge he : f.adjacentHalfedges()) {
        std::vector<double> isoPoints;
        if (he.next() == f.halfedge())
          isoPoints = crossingsModulo2Pi(stripeValues[he.corner()],
                                         stripeValues[he.next().corner()] + 2 * stripeIndices[f] * PI);
        else
          isoPoints = crossingsModulo2Pi(stripeValues[he.corner()], stripeValues[he.next().corner()]);

        std::vector<Vector3> points;
        for (double bary : isoPoints) {
          points.push_back(bary * geometry.vertexPositions[he.tailVertex()] +
                           (1 - bary) * geometry.vertexPositions[he.tipVertex()]);
        }
        if (polylineIndices[he.edge()].size() > 0) { // use the existing indices in reverse order
          std::vector<size_t> reverseIndices(polylineIndices[he.edge()].rbegin(), polylineIndices[he.edge()].rend());
          faceIndices.push_back(reverseIndices);
        } else { // add the intersections to the vertices lists and initialize the corresponding list of indices
          for (Vector3 point : points) {
            polylineIndices[he.edge()].push_back(vertices.size());
            vertices.push_back(point);
          }
          faceIndices.push_back(polylineIndices[he.edge()]);
        }
        facePoints.push_back(points);
      }

      // compute minimum of abs(dot(v - w, vectorField[f])) for all pairs v, w
      std::vector<std::pair<double, std::array<size_t, 2>>> minValues;
      for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < facePoints[i].size(); ++j) {
          Vector3 v = facePoints[i][j];
          double min = findMin(facePoints, v, vectorField[f], i);
          minValues.push_back({min, {i, j}});
        }
      }

      // remove the indices corresponding to the largest elements in minValues
      int N = abs(stripeIndices[f]);
      std::vector<std::vector<size_t>> crossingIndices = removeCrossings(minValues, N, faceIndices);

      std::vector<std::array<size_t, 2>> matchings = matchCrossings(crossingIndices);
      edges.insert(edges.end(), matchings.begin(), matchings.end());
    }
  }
}

} // namespace

std::tuple<std::vector<Vector3>, std::vector<std::array<size_t, 2>>>
extractPolylinesFromStripePattern(EmbeddedGeometryInterface& geometry, const CornerData<double>& values,
                                  const FaceData<int>& stripeIndices, const FaceData<int>& fieldIndices,
                                  const VertexData<Vector2>& directionField, bool connectOnSingularities) {
  SurfaceMesh& mesh = geometry.mesh;

  std::vector<Vector3> points;
  std::vector<std::array<size_t, 2>> edges;
  EdgeData<std::vector<size_t>> indicesPerEdge;
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


std::tuple<std::vector<Vector3>, std::vector<std::array<size_t, 2>>>
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