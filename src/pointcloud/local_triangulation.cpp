#include "geometrycentral/pointcloud/local_triangulation.h"


#include "geometrycentral/utilities/elementary_geometry.h"


namespace geometrycentral {
namespace pointcloud {

PointData<std::vector<std::array<Point, 3>>> buildLocalTriangulations(PointCloud& cloud, PointPositionGeometry& geom,
                                                                      bool withDegeneracyHeuristic) {


  // NOTE: This is not robust if the entire neighbohood is coincident (or very nearly coincident) with the centerpoint.
  // Though in that case, the generating normals will probably also have issues.

  // An innocent numerical parameter used for the degeneracy heuristic
  const double DEGENERATE_THRESH = 1e-7; // in units of relative length

  geom.requireNeighbors();
  geom.requireTangentCoordinates();

  PointData<std::vector<std::array<Point, 3>>> result(cloud);

  for (Point p : cloud.points()) {
    size_t nNeigh = geom.neighbors->neighbors[p].size();


    double lenScale = 0;
    { // Compute a lengthscale for the neighborhood as the radius of the most distant point
      double lenScale2 = 0;
      for (size_t iNeigh = 0; iNeigh < nNeigh; iNeigh++) {
        Vector2 neighPt = geom.tangentCoordinates[p][iNeigh];
        double dist2 = norm2(neighPt);
        lenScale2 = std::fmax(lenScale2, dist2);
      }
      lenScale = std::sqrt(lenScale2);
    }

    // Something is hopelessly degenerate, don't even bother trying. No triangles for this point.
    if (!std::isfinite(lenScale) || lenScale <= 0) {
      std::cerr << "skipping degenerate neighborhood" << std::endl;
      continue;
    }

    // Local copies of points
    std::vector<Vector2> perturbPoints = geom.tangentCoordinates[p];


    if (withDegeneracyHeuristic) {

      { // Perturb points which are extremely close to the source
        for (size_t iNeigh = 0; iNeigh < nNeigh; iNeigh++) {
          Vector2& neighPt = perturbPoints[iNeigh];
          double dist = norm(neighPt);
          if (dist < lenScale * DEGENERATE_THRESH) { // need to perturb
            Vector2 dir = normalize(neighPt);
            if (!isfinite(dir)) { // even direction is degenerate :(
              // pick a direction from index
              double thetaDir = (2. * M_PI * iNeigh) / nNeigh;
              dir = Vector2::fromAngle(thetaDir);
            }

            // Set the distance from the origin for the pertubed point. Including the index avoids creating many
            // co-circular points; no need to stress the Delaunay triangulation unnessecarily.
            double len = (1. + static_cast<double>(iNeigh) / nNeigh) * lenScale * DEGENERATE_THRESH * 10;

            neighPt = len * dir; // update the point
          }
        }
      }
    }

    std::vector<size_t> sortInds;
    { // = Angularly sort the points CCW, such that the closest point comes first

      // sentinel value for below
      double BAD_ANGLE = -777;

      std::vector<double> pointAngles;
      for (size_t i = 0; i < nNeigh; i++) {
        double angle = arg(unit(perturbPoints[i]));
        if (!std::isfinite(angle)) {
          angle = BAD_ANGLE;
        }
        sortInds.push_back(i);
        pointAngles.push_back(angle);
      }

      // Angular sort
      std::sort(sortInds.begin(), sortInds.end(),
                [&](const size_t& a, const size_t& b) -> bool { return pointAngles[a] < pointAngles[b]; });

      // Immediately skip any invalid indices in the search below, by detecting sentinels from above
      for (size_t i = 0; i < nNeigh; i++) {
        if (pointAngles[sortInds[i]] == BAD_ANGLE) {
          sortInds[i] = INVALID_IND;
        }
      }
    }


    // == Find the local Delaunay triangulation
    // Strategy: we start with an angularly-sorted list of points around the center point. The output we seek is a
    // subset of this list, which corresponds to the 1-ring of the center vertex in the local Delaunay triangulation. We
    // will repeatedly remove points from the sorted list (marking them as removed with INVALID_IND) if the diamond they
    // form is not Delaunay. Note that is is also necessary to handle the case where all points lie in some half-space,
    // there is an absent triangle between some pair of points.
    // NOTE: Doing this with while(changed) {} loop as below is vulnerable to N^2 behavior on e.g. a spiral
    // configuration. But that seems sufficiently unlikely in practice. Could be remedied with a queue.

    auto isBoundary = [&](size_t localIndA, size_t localIndB) {
      Vector2 pA = perturbPoints[localIndA];
      Vector2 pB = perturbPoints[localIndB];

      return cross(pA, pB) <= 0.;
    };

    // The main loop, discarding (aka flipping) until we have the Delaunay triagulation
    bool anyChanged = true;
    size_t N = sortInds.size();

    while (anyChanged) {
      anyChanged = false;
      for (size_t iMiddle = 0; iMiddle < N; iMiddle++) {
        if (sortInds[iMiddle] == INVALID_IND) continue; // skip unused indices

        // Find the previous and next indices
        size_t iPrev = iMiddle;
        do {
          iPrev = (iPrev + N - 1) % N;
        } while (sortInds[iPrev] == INVALID_IND);
        size_t iNext = iMiddle;
        do {
          iNext = (iNext + 1) % N;
        } while (sortInds[iNext] == INVALID_IND);

        // Indices in to the local neighbor list
        size_t prev = sortInds[iPrev];
        size_t curr = sortInds[iMiddle];
        size_t next = sortInds[iNext];

        // Degenerate cases
        if (curr == prev || curr == next || prev == next) continue;

        // For any colllinear points, keep only the closest
        if (withDegeneracyHeuristic) {
          double lenPrev = norm(perturbPoints[prev]);
          double lenCurr = norm(perturbPoints[curr]);
          double lenNext = norm(perturbPoints[next]);

          bool collinearPrev =
              std::abs(cross(perturbPoints[curr], perturbPoints[prev])) < (lenPrev * lenCurr) * DEGENERATE_THRESH &&
              dot(perturbPoints[curr], perturbPoints[prev]) > 0;
          bool collinearNext =
              std::abs(cross(perturbPoints[curr], perturbPoints[next])) < (lenNext * lenCurr) * DEGENERATE_THRESH &&
              dot(perturbPoints[curr], perturbPoints[next]) > 0;

          if ((collinearNext && lenCurr > lenNext) || (collinearPrev && lenCurr > lenPrev)) {
            sortInds[iMiddle] = INVALID_IND;
            anyChanged = true;
            continue;
          }
        }

        // If either of the triangles is empty (aka actually the boundary), skip this
        if (isBoundary(prev, curr) || isBoundary(curr, next)) continue;

        // Test if the triangles should be merged
        if (!inCircleTest(Vector2{0., 0.}, perturbPoints[prev], perturbPoints[next], perturbPoints[curr])) {
          sortInds[iMiddle] = INVALID_IND;
          anyChanged = true;
        }
      }
    }

    // Emit the actual triangles
    std::vector<std::array<Point, 3>>& thisPointTriangles = result[p];
    for (size_t iPrev = 0; iPrev < N; iPrev++) {
      if (sortInds[iPrev] == INVALID_IND) continue; // skip unused indices

      size_t iNext = iPrev;
      do {
        iNext = (iNext + 1) % N;
      } while (sortInds[iNext] == INVALID_IND);

      if (iPrev == iNext) continue;

      // Indices in to the local neighbor list
      size_t prev = sortInds[iPrev];
      size_t next = sortInds[iNext];

      if (!isBoundary(prev, next)) {
        thisPointTriangles.push_back(
            std::array<Point, 3>{p, geom.neighbors->neighbors[p][prev], geom.neighbors->neighbors[p][next]});
      }
    }
  }

  geom.unrequireNeighbors();
  geom.unrequireTangentCoordinates();

  return result;
}

PointData<std::vector<std::array<size_t, 3>>>
handleToInds(PointCloud& cloud, const PointData<std::vector<std::array<Point, 3>>>& localResult) {

  GC_SAFETY_ASSERT(cloud.isCompressed(), "cloud must be compressed");

  PointData<std::vector<std::array<size_t, 3>>> result(cloud);
  for (Point p : cloud.points()) {
    size_t nTri = localResult[p].size();
    result[p].resize(nTri);
    for (size_t i = 0; i < nTri; i++) {
      for (size_t j = 0; j < 3; j++) {
        result[p][i][j] = localResult[p][i][j].getIndex();
      }
    }
  }

  return result;
}

std::vector<std::vector<size_t>> handleToFlatInds(PointCloud& cloud,
                                                    const PointData<std::vector<std::array<Point, 3>>>& localResult) {

  GC_SAFETY_ASSERT(cloud.isCompressed(), "cloud must be compressed");

  std::vector<std::vector<size_t>> result;
  for (Point p : cloud.points()) {
    size_t nTri = localResult[p].size();
    for (size_t i = 0; i < nTri; i++) {
      std::vector<size_t> tri(3);
      for (size_t j = 0; j < 3; j++) {
        tri[j] = localResult[p][i][j].getIndex();
      }
      result.push_back(tri);
    }
  }

  return result;
}

} // namespace pointcloud
} // namespace geometrycentral
