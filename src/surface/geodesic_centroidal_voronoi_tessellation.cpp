#include "geometrycentral/surface/geodesic_centroidal_voronoi_tessellation.h"

namespace geometrycentral {
namespace surface {

namespace {
const bool VORONOI_PRINT = false;
}

// The default trace options
const VoronoiOptions defaultVoronoiOptions;

VoronoiResult computeGeodesicCentroidalVoronoiTessellation(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom,
                                                           VoronoiOptions options) {

  if (options.useDelaunay) {
    std::unique_ptr<SignpostIntrinsicTriangulation> intTri(new SignpostIntrinsicTriangulation(mesh, geom));
    intTri->flipToDelaunay();

    // Translate initial sites to intrinsic triangulation if any were given
    for (size_t iS = 0; iS < options.initialSites.size(); iS++) {
      options.initialSites[iS] = intTri->equivalentPointOnIntrinsic(options.initialSites[iS]);
    }

    // Get solutions on intrinsic triangulation
    options.useDelaunay = false;
    VoronoiResult result = computeGeodesicCentroidalVoronoiTessellation(*intTri->intrinsicMesh, *intTri, options);

    // Translate solutions back to original triangulation
    for (size_t iS = 0; iS < result.siteLocations.size(); iS++) {
      result.siteLocations[iS] = intTri->equivalentPointOnInput(result.siteLocations[iS]);
    }

    return result;
  }

  VectorHeatMethodSolver vSolver(geom, options.tCoef);

  auto computeRHS = [&](const std::vector<SurfacePoint>& points) -> VertexData<double> {
    VertexData<double> rhs(mesh, 0);
    for (const SurfacePoint& point : points) {
      SurfacePoint facePoint = point.inSomeFace();
      Halfedge he = facePoint.face.halfedge();

      rhs[he.vertex()] += facePoint.faceCoords.x;
      rhs[he.next().vertex()] += facePoint.faceCoords.y;
      rhs[he.next().next().vertex()] += facePoint.faceCoords.z;
    }
    return rhs;
  };

  // Set points to start
  std::vector<SurfacePoint> siteLocations = options.initialSites;
  if (siteLocations.empty()) {
    for (size_t i = 0; i < options.nSites; i++) {
      Face startF = mesh.face(randomIndex(mesh.nFaces()));
      double u = unitRand();
      Vector3 bCoord{u, 0.5 * (1.0 - u), 0.5 * (1.0 - u)};
      SurfacePoint fp{startF, bCoord};
      siteLocations.push_back(fp);
    }
  }
  size_t nSites = siteLocations.size();

  geom.requireVertexDualAreas();

  // == Iterations
  for (size_t iIter = 0; iIter < options.iterations; iIter++) {

    // Compute the normalizer distribution
    VertexData<double> normRHS = computeRHS(siteLocations);
    VertexData<double> normD = vSolver.scalarDiffuse(normRHS);

    double energy = 0;
    std::vector<SurfacePoint> newSiteLocations;

    for (size_t iSite = 0; iSite < nSites; iSite++) {
      SurfacePoint site = siteLocations[iSite];

      // === Compute the nearest distribution
      VertexData<double> unitRHS = computeRHS({site});
      VertexData<double> thisFracD = vSolver.scalarDiffuse(unitRHS);
      for (Vertex v : mesh.vertices()) {
        thisFracD[v] /= normD[v];
      }

      for (size_t iSubIter = 0; iSubIter < options.nSubIterations; iSubIter++) {

        // === Compute the log map
        VertexData<Vector2> logmap = vSolver.computeLogMap(site);

        // Evaluate energy and gradient contribution
        Vector2 updateSum{0, 0};
        double updateWSum = 0.0;

        for (Vertex v : mesh.vertices()) {

          double weight = thisFracD[v] * geom.vertexDualAreas[v];
          Vector2 logVal = logmap[v];
          double dist = logVal.norm();

          updateSum += weight * logVal;
          updateWSum += weight;

          energy += dist * dist * weight;
        }
        updateSum /= updateWSum;
        Vector2 update = updateSum / updateWSum;

        // Take a step
        TraceGeodesicResult traceResult = traceGeodesic(geom, site, options.stepSize * update);
        site = traceResult.endPoint;
      }

      newSiteLocations.push_back(site);
    }

    siteLocations = newSiteLocations;
    if (VORONOI_PRINT) std::cout << "Finished iteration " << iIter << "  energy " << energy << std::endl;
  }

  geom.unrequireVertexDualAreas();

  VoronoiResult result;
  result.siteLocations = siteLocations;

  if (options.computeDistributions) {
    result.hasDistributions = true;

    // Compute the normalizer distribution
    VertexData<double> normRHS = computeRHS(siteLocations);
    VertexData<double> normD = vSolver.scalarDiffuse(normRHS);

    for (size_t iSite = 0; iSite < nSites; iSite++) {
      SurfacePoint site = siteLocations[iSite];

      // === Compute the nearest distribution
      VertexData<double> unitRHS = computeRHS({site});
      VertexData<double> thisFracD = vSolver.scalarDiffuse(unitRHS);
      for (Vertex v : mesh.vertices()) thisFracD[v] /= normD[v];
      result.siteDistributions.push_back(thisFracD);
    }
  }

  return result;
}
} // namespace surface
} // namespace geometrycentral
