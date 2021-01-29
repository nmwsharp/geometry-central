# Heat distance, transport, & logarithmic map on point clouds

Compute geodesic distance, transport tangent vectors, and generate a special parameterization called the logarithmic map using fast solvers based on short-time heat flow.

![point cloud heat solve results](/media/point_heat_solvers.jpg)

These routines implement point cloud versions of the algorithms from:

- [The Heat Method for Distance Computation](http://www.cs.cmu.edu/~kmcrane/Projects/HeatMethod/index.html) (distance)
- [The Vector Heat Method](https://nmwsharp.com/research/vector-heat-method) (parallel transport and log map)
- [A Laplacian for Nonmanifold Triangle Meshes](http://www.cs.cmu.edu/~kmcrane/Projects/NonmanifoldLaplace/NonmanifoldLaplace.pdf) (used to build point cloud Laplacian for both)



All computation is encapsulated by the `PointCloudHeatSolver` class, which maintains prefactored linear systems for the various methods. Some setup work is performed both on construction, and after the first query. Subsequent queries, even with different source points, will be fast.

`#include "geometrycentral/pointcloud/point_cloud_heat_solver.h"`


**Example:** Basic usage

```cpp
#include "geometrycentral/pointcloud/point_cloud.h"
#include "geometrycentral/pointcloud/point_position_geometry.h"
#include "geometrycentral/pointcloud/point_cloud_heat_solver.h"
#include "geometrycentral/pointcloud/point_cloud_io.h"

using namespace geometrycentral;
using namespace geometrycentral::pointcloud;

// Read in a point cloud
std::unique_ptr<PointCloud> cloud;
std::unique_ptr<PointPositionGeometry> geom;
std::tie(cloud, geom) = readPointCloud("my_cloud.ply");

// Create the solver
PointCloudHeatSolver solver(*cloud, *geom);

// Pick a source point or two
Point pSource = cloud->point(7);
Point pSource2 = cloud->point(8);

// Compute geodesic distance
PointData<double> distance = solver.computeDistance(pSource);

// Compute scalar extension
PointData<double> extended = solver.extendScalars({{pSource,  3.},
                                                   {pSource2, -5.}});

// Compute parallel transport
Vector2 sourceVec{1, 2};
PointData<Vector2> transport = solver.transportTangentVector(pSource, sourceVec);

// Compute the logarithmic map
PointData<Vector2> logmap = solver.computeLogMap(pSource);
```

### Constructor

??? func "`#!cpp PointCloudHeatSolver::PointCloudHeatSolver(PointCloud& cloud, PointPositionGeometry& geom, double tCoef = 1.0)`"

    Create a new solver for the heat methods. Precomputation is performed at startup and lazily as needed.

    - `cloud` is the point cloud to operate on.
    - `geom` is the geometry (e.g., vertex positions) for the cloud. Note that nearly any point cloud geometry object can be passed here. For instance, use a `PointPositionFrameGeometry` if you want to represent tangent-valued data in tangent bases of your own choosing.

    - `tCoef` is the time to use for short time heat flow, as a factor `m * h^2`, where `h` is the mean between-point spacing. The default value of `1.0` is almost always sufficient.

    Algorithm options (like `tCoef`) cannot be changed after construction; create a new solver object with the new settings.


## Geodesic distance

_Geodesic distance_ is the distance from a given source along the surface represented by the point cloud. Specifying multiple source points yields the distance to the nearest source.

??? func "`#!cpp PointData<double> PointCloudHeatSolver::computeDistance(const Point& sourcePoint)`"

    Compute the geodesic distance from `sourcePoint` to all other points.

??? func "`#!cpp PointData<double> PointCloudHeatSolver::computeDistance(const std::vector<Point>& sourcePoints)`"

    Like above, but for multiple source points.


## Scalar extension

Given scalar values defined at isolated points in the cloud, extend it to a scalar field on all points. Each point will take the value from the nearest source point, in the sense of geodesic distance. Note that the fast diffusion algorithm means the result is a slightly smoothed-out field.

??? func "`#!cpp PointData<double> PointCloudHeatSolver::extendScalars(const std::vector<std::tuple<Point, double>>& sources)`"

    Given a collection of source points and scalars at those points, extends the scalar field to the whole cloud as a nearest-geodesic-neighbor interpolant.


## Vector extension

Given tangent vectors defined at one or more isolated source locations on a surface, extend transport the vectors across the entire domain according to parallel transport. Each point on the domain will take the value of the nearest source point.  Note that the fast diffusion algorithm means the result is a slightly smoothed-out field.

Despite operating on tangent-valued data, these routines work as expected even on a point cloud with inconsistent normal orientations. See [connection Laplacian](/pointcloud/geometry/#connection-laplacian) for more details.

??? func "`#!cpp PointData<Vector2> PointCloudHeatSolver::transportTangentVector(const Point& sourcePoint, const Vector2& sourceVector)`"

    Shorthand for the general version below when there is just one vector.

    Computes parallel transport of the given vector along shortest geodesics to the rest of the cloud.

    Polar directions are defined in each point's tangent space. [See `PointPositionGeometry::tangentBasis`](/pointcloud/geometry#tangent-basis).

??? func "`#!cpp PointData<Vector2> PointCloudHeatSolver::transportTangentVectors(const std::vector<std::tuple<Point, Vector2>>& sources)`"

    Given a collection of source points and tangent vectors at those points, extends the vector field to the whole cloud as a nearest-geodesic-neighbor interpolant via parallel transport along shortest geodesics.
    
    Polar directions are defined in each point's tangent space. [See `PointPositionGeometry::tangentBasis`](/pointcloud/geometry#tangent-basis).

## Logarithmic Map

The _logarithmic map_ is a very special 2D local parameterization of a surface about a point, where for each point on the surface the magnitude of the log map gives the geodesic distance from the source, and the polar coordinate of the log map gives the direction at which a geodesic must leave the source to arrive at the point.

Despite operating on tangent-valued data, this routine works as expected even on a point cloud with inconsistent normal orientations. See [connection Laplacian](/pointcloud/geometry/#connection-laplacian) for more details.
  
??? func "`#!cpp PointData<Vector2> PointCloudHeatSolver::computeLogMap(const Point& sourcePoint)`"

    Computes the logarithmic map from the given source point.

    Polar directions are defined in each point's tangent space. [See `PointPositionGeometry::tangentBasis`](/pointcloud/geometry#tangent-basis).
