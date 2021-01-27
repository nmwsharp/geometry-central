The `PointCloud` class is the main data structure for representing point clouds in geometry-central. It mimics the design of `SurfaceMesh` for the sake of consistency, though it is generally much simpler. The `PointCloud` class is simply a logical object referring to an abstract collection of points.  Geometric data (aka point positions) are stored separately in the `PointPositionGeometry` class (or related subclasses).

```#include "geometrycentral/pointcloud/point_cloud.h"```

!!! warning "Point clouds are in beta"

    The current point cloud API in geometry-central is preliminary, and may change in future versions.


**Example:** A tour of basic point cloud functionality.
```c++
#include "geometrycentral/pointcloud/point_cloud.h"
#include "geometrycentral/pointcloud/point_position_geometry.h"

size_t nPts = 5000;

// create a new cloud
PointCloud cloud(nPts); 

// access properties of the cloud
std::cout << "cloud has " << cloud.nPoints() << " points.\n";
// Prints "cloud has 5000 points"

// iterate through the points in the cloud
for(Point p : cloud.points()) {
  std::cout << "Hi, I'm point " << p << "\n";
  // Prints "Hi, I'm point p_245" (etc)
}

// store data on a point cloud
PointData<int> values(cloud);
for(Point p : cloud.points()) {
  values[p] = 7;
}

// associate 3D positions with the point cloud
// compute geometric data using the geomety object
PointData<Vector3> positions = /* some positions */;
PointPositionGeometry geom(*cloud, positions);

// for example, compute normals
geom.requireNormals();
for(Point p : cloud.points()) {
  std::cout << "normal for point " << p << " is " << geom.normals[p] << "\n";
}
```


??? func "`#!cpp PointCloud::PointCloud(size_t nPts)`"
    Constructs a new point cloud with the desired number of points.

    See the IO routines for higher-level constructors to simultaneously construct a point cloud and its geometry, or read from file, etc.


??? func "`#!cpp /*iterator type*/ PointCloud::points()`"

    Iterate through the points in a point cloud, like

    ```c++
    for(Point p : cloud.points()) {
      // do science
    }
    ```

??? func "`#!cpp Point PointCloud::point(size_t iP)`"

    Return a handle to the i'th point.

??? func "`#!cpp size_t PointCloud::nPoints() const`"

    Returns the number of points. 

??? func "`#!cpp std::unique_ptr<PointCloud> PointCloud::copy() const`"

    Returns the number of points. 

!!! note "Mutating point clouds"

    TODO The `PointCloud` class contains the groundwork for efficiently adding/removing points from the cloud while updating containers (etc) similar to `SurfaceMesh`. However, this functionality is not fleshed out, tested, or documented yet.
