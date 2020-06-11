#include "geometrycentral/surface/surface_centers.h"

#include "geometrycentral/surface/trace_geodesic.h"
#include "geometrycentral/surface/vector_heat_method.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "geometrycentral/utilities/utilities.h"

namespace geometrycentral {
namespace surface {


SurfacePoint findCenter(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom, const std::vector<Vertex>& vertexPts, int p) {
  VertexData<double> dist(geom.mesh, 0.);
  for (Vertex v : vertexPts) {
    dist[v] += 1.;
  }

  // Forward to the general version
  return findCenter(mesh, geom, dist, p);
}

SurfacePoint findCenter(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom, VectorHeatMethodSolver& solver,
                        const std::vector<Vertex>& vertexPts, int p) {
  VertexData<double> dist(geom.mesh, 0.);
  for (Vertex v : vertexPts) {
    dist[v] += 1.;
  }

  // Forward to the general version
  return findCenter(mesh, geom, solver, dist, p);
}

SurfacePoint findCenter(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom, const VertexData<double>& distribution, int p) {
  VectorHeatMethodSolver solver(geom);

  // Forward to the general version
  return findCenter(mesh, geom, solver, distribution, p);
}


SurfacePoint findCenter(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom, VectorHeatMethodSolver& solver,
                        const VertexData<double>& distribution, int p) {

  if (p != 1 && p != 2) {
    throw std::logic_error("only p=1 or p=2 is supported");
  }

  geom.requireFaceAreas();
  geom.requireHalfedgeVectorsInVertex();
  geom.requireHalfedgeVectorsInFace();

  SurfacePoint initialGuess;

  // Reasonable initial guess: the vertex with the most stuff at it
  // double biggestVal = -std::numeric_limits<double>::infinity();
  // for (Vertex v : mesh.vertices()) {
  // if (distribution[v] > biggestVal) {
  // biggestVal = distribution[v];
  // initialGuess = SurfacePoint(v);
  //}
  //}

  // Random initial guess
  initialGuess = SurfacePoint(mesh.vertex(randomIndex(mesh.nVertices())));

  // Compute an approximate mesh diameter to use as parameter in the p=1 case
  double meshDiameter = 1.;
  if (p == 1) {
    double surfaceArea = 0.;
    for (Face f : mesh.faces()) {
      surfaceArea += geom.faceAreas[f];
    }
    meshDiameter = std::sqrt(surfaceArea);
  }


  // = A few parameters
  int maxIters = 100;
  double initialStepSize = 1.0;
  double p1DistanceEps = 1e-9 * meshDiameter; // soften the fraction in p = 1 case to avoid divide by 0
  double convergeThresh = 1 / 3.;             // convergence once step is this fraction of face size
  double faceMapThreshold = 10.;              // switch from using per-vertex maps to sub-vertex
                                              //    accurate maps. units of local face length scale.


  auto evalEnergyAndUpdate = [&](SurfacePoint aboutPoint) -> std::tuple<double, Vector2> {
    // Compute the current log map
    // Solve at the face point
    VertexData<Vector2> logmap = solver.computeLogMap(aboutPoint);

    // Evaluate energy and update step
    double thisEnergy = 0.;
    Vector2 thisUpdate = Vector2::zero();
    double updateWeightSum = 0.;
    for (Vertex v : mesh.vertices()) {
      if (distribution[v] == 0.) continue;

      Vector2 pointCoord = logmap[v];
      double dist2 = pointCoord.norm2();

      if (p == 1) {
        double dist = std::sqrt(dist2);
        double w = distribution[v];

        // Energy
        thisEnergy += dist * w;

        // Update
        thisUpdate += w * pointCoord / (dist + p1DistanceEps);
        updateWeightSum += w / (dist + p1DistanceEps);

      } else if (p == 2) {
        double w = distribution[v];

        // Energy
        thisEnergy += dist2 * w;

        // Update
        thisUpdate += w * pointCoord;
        updateWeightSum += w;
      }
    }
    thisUpdate /= updateWeightSum;

    return std::make_tuple(thisEnergy, thisUpdate);
  };


  // === Perform Weiszfeld iterations
  SurfacePoint currCenter = initialGuess.inSomeFace();

  bool converged = false;
  for (int i = 0; i < maxIters; i++) {

    // Compute energy and update
    double energyBefore;
    Vector2 updateVec;
    std::tie(energyBefore, updateVec) = evalEnergyAndUpdate(currCenter);

    // Line search
    double stepSize = initialStepSize;
    for (int lineSearchIter = 0; lineSearchIter < 8; lineSearchIter++) {

      // Try taking a step
      Vector2 stepVec = updateVec * stepSize;
      TraceOptions options;
      options.includePath = true; // TODO why?
      TraceGeodesicResult traceResult = traceGeodesic(geom, currCenter, stepVec, options);
      SurfacePoint candidatePoint = traceResult.endPoint.inSomeFace();

      // Compute new energy
      double newEnergy = std::get<0>(evalEnergyAndUpdate(candidatePoint));

      // Check for convergence
      double faceScale = std::sqrt(geom.faceAreas[currCenter.inSomeFace().face]);
      if (stepVec.norm() < convergeThresh * faceScale) {
        converged = true;
        break;
      }

      // Accept step if good
      if (newEnergy < energyBefore) {
        currCenter = candidatePoint;
        break;
      }

      // Otherwise decrease step size and repeat
      stepSize *= 0.5;
    }

    if (converged) {
      break;
    }
  }

  geom.unrequireFaceAreas();
  geom.unrequireHalfedgeVectorsInVertex();
  geom.unrequireHalfedgeVectorsInFace();

  return currCenter;
}

} // namespace surface
} // namespace geometrycentral
