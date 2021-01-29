# Input & Output

Some basic input/output routines for point clouds.

## Input

**Example:** Read in a point cloud.

```cpp
#include "geometrycentral/pointcloud/point_cloud.h"
#include "geometrycentral/pointcloud/point_position_geometry.h"
#include "geometrycentral/pointcloud/point_cloud_io.h"

using namespace geometrycentral;
using namespace geometrycentral::pointcloud;

std::unique_ptr<PointCloud> cloud;
std::unique_ptr<PointPositionGeometry> geom;
std::tie(cloud, geom) = readPointCloud("my_cloud.ply");

// do science!
```

??? func "`#!cpp std::tuple<std::unique_ptr<PointCloud>, std::unique_ptr<PointPositionGeometry>> readPointCloud(std::string filename, std::string type = "")`"

    Read a point cloud from file, constructing both the cloud and geometry objects.

    Currently accepted file types: `ply`, `obj`. Using the default empty type string will attempt to infer from the filename.

??? func "`#!cpp std::tuple<std::unique_ptr<PointCloud>, std::unique_ptr<PointPositionGeometry>> readPointCloud(std::istream& in, std::string type)`"

    Like above, but reads directly from an `istream`. The type must be specified explicitly.


## Output

**Example:** Write out a point cloud.

```cpp
#include "geometrycentral/pointcloud/point_cloud.h"
#include "geometrycentral/pointcloud/point_position_geometry.h"
#include "geometrycentral/pointcloud/point_cloud_io.h"

using namespace geometrycentral;
using namespace geometrycentral::pointcloud;

// your cloud & point positions, populated somehow
std::unique_ptr<PointCloud> cloud;
std::unique_ptr<PointPositionGeometry> geom;

writePointCloud(*cloud, *geom, "cloud_out.ply");
```

??? func "`#!cpp void writePointCloud(PointCloud& cloud, PointPositionGeometry& geometry, std::string filename, std::string type = "")`"

    Write a point cloud to file.

    Currently accepted file types: `ply`, `obj`. Using the default empty type string will attempt to infer from the filename.

??? func "`#!cpp void writePointCloud(PointCloud& cloud, PointPositionGeometry& geometry, std::ostream& out, std::string type)`"
    
    Like above, but writes directly to an `ostream`. The type must be specified explicitly.
