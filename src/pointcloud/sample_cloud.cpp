#include "geometrycentral/pointcloud/sample_cloud.h"

#include <random>


namespace geometrycentral {

using surface::Face;
using surface::SurfacePoint;

namespace pointcloud {

std::tuple<std::unique_ptr<PointCloud>, PointData<Vector3>, PointData<SurfacePoint>>
uniformlySamplePointsOnSurface(surface::SurfaceMesh& mesh, surface::EmbeddedGeometryInterface& geom, size_t nPts) {


  geom.requireVertexPositions();
  geom.requireFaceAreas();

  // Random number generator
  // if we ever want to seed, do it here
  std::random_device rd;
  std::mt19937 gen(rd());

  // Create a sampler to pick points from faces
  // (uniformly in the geometric sense)
  std::vector<double> areas;
  std::vector<Face> faces;
  areas.reserve(mesh.nFaces());
  faces.reserve(mesh.nFaces());
  for (Face f : mesh.faces()) {
    GC_SAFETY_ASSERT(f.isTriangle(), "can only sample point cloud from triangular mesh");
    areas.push_back(geom.faceAreas[f]);
    faces.push_back(f);
  }
  std::discrete_distribution<size_t> faceDist(areas.begin(), areas.end());

  // Create a real-valued sampler, which we will use for barycentric coordinates within faces
  std::uniform_real_distribution<double> realDist(0., 1.);

  // Store results here
  std::unique_ptr<PointCloud> cloud(new PointCloud(nPts));
  PointData<Vector3> pos(*cloud);
  PointData<SurfacePoint> cloudSources(*cloud);

  // Sample
  for (size_t iSample = 0; iSample < nPts; iSample++) {
    Point p = cloud->point(iSample);

    // Pick a face
    size_t iF = faceDist(gen);
    Face f = faces[iF];

    // Pick barycentric coordinates within the face
    double r1 = realDist(gen);
    double r2 = realDist(gen);
    SurfacePoint surfP(f, Vector3{1. - std::sqrt(r1), std::sqrt(r1) * (1. - r2), std::sqrt(r1) * r2});

    // Interpolate to get the position of the sampled point
    Vector3 newPos = surfP.interpolate(geom.vertexPositions);

    pos[p] = newPos;
    cloudSources[p] = surfP;
  }


  geom.unrequireVertexPositions();
  geom.unrequireFaceAreas();

  return std::make_tuple(std::move(cloud), pos, cloudSources);
}


} // namespace pointcloud
} // namespace geometrycentral
