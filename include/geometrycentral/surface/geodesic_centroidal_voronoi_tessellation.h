#pragma once

#include "geometrycentral/surface/signpost_intrinsic_triangulation.h"
#include "geometrycentral/surface/trace_geodesic.h"
#include "geometrycentral/surface/vector_heat_method.h"

namespace geometrycentral {
namespace surface {

struct VoronoiResult {
  std::vector<SurfacePoint> siteLocations;           // sites at the centers of the Voronoi cells
  std::vector<VertexData<double>> siteDistributions; // soft indicator functions for each Voronoi cell
  bool hasDistributions = false;                     // is siteDistributions populated?
};

struct VoronoiOptions {
  size_t nSites = 10;                     // number of sites to place
  std::vector<SurfacePoint> initialSites; // desired locations for sites. If blank, locations are chosen randomly
  size_t iterations = 50;                 // number of iterations to run for
  double stepSize = 1;                    // step size for steps towards cell centers
  bool useDelaunay = true;                // solve on an intrinsic Delaunay triangulation of the input
  bool computeDistributions = false;      // return the indicator functions for each cell (`result.siteDistributions`)
  double tCoef = 1;                       // diffusion time for vector heat method
  size_t nSubIterations = 1;              // number of iterations to use when computing Karcher means
};
extern const VoronoiOptions defaultVoronoiOptions;

VoronoiResult computeGeodesicCentroidalVoronoiTessellation(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom,
                                                           VoronoiOptions options = defaultVoronoiOptions);

} // namespace surface
} // namespace geometrycentral
