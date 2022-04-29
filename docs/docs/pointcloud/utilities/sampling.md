Sample random point clouds from triangle meshes.

**Example:** Read a triangle mesh from file and uniformly sample points
```cpp
#include "geometrycentral/pointcloud/point_cloud.h"
#include "geometrycentral/pointcloud/point_position_geometry.h"
#include "geometrycentral/pointcloud/point_cloud_io.h"
#include "geometrycentral/pointcloud/point_cloud.h"
#include "geometrycentral/surface/meshio.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;
using namespace geometrycentral::pointcloud;

// Load a mesh
std::unique_ptr<SurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> meshGeom;
std::tie(mesh, meshGeom) = loadMesh(args::get(inputFilename));

// Sample a point cloud from the mesh
std::unique_ptr<PointCloud> cloud;
PointData<Vector3> pointPos;
PointData<SurfacePoint> cloudSources;
size_t nPts = 5000;
std::tie(cloud, pointPos, cloudSources) = uniformlySamplePointsOnSurface(*mesh, *meshGeom, nPts);

// As an example, use the source points to get face normals
PointData<Vector3> normals(*cloud);
meshGeom->requireFaceNormals();
for (Point p : cloud->points()) {
  normals[p] = meshGeom->faceNormals[cloudSources[p].face];
}

// Construct a geometry object from the positions for subsequent calculations
PointPositionGeometry cloudGeom(*cloud, pointPos);
```

??? func "`#!cpp std::tuple<std::unique_ptr<PointCloud>, PointData<Vector3>, PointData<surface::SurfacePoint>> uniformlySamplePointsOnSurface(surface::SurfaceMesh& mesh, surface::EmbeddedGeometryInterface& geom, size_t nPts)`"
   
    Sample `nPts` points from a triangle mesh. Points are sampled uniformly at random from the underlying surface, independent of tesselation.

    The return tuple holds:

    - The new `PointCloud`
    - A 3D position for each point
    - A `SurfacePoint` on the original mesh corresponding to each point


Note that in addition to the cloud and 3D positions, this routine returns a `SurfacePoint` associated with each point sample. The [surface point](/surface/utilities/surface_point/) is a handy class representing a location on the underlying mesh, and makes it easy to do things like interpolate data defined on the source mesh, or grab face normals, etc.

