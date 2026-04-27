A simple utility class for 3D nearest-neighbor queries, backed by a [nanoflann](https://github.com/jlblancoc/nanoflann) KD-tree.

`#!cpp #include "geometrycentral/utilities/knn.h"`

```cpp
#include "geometrycentral/utilities/knn.h"
using namespace geometrycentral;

std::vector<Vector3> points = { {0,0,0}, {1,0,0}, {2,0,0}, {3,0,0} };
NearestNeighborFinder finder(points);

// 2 nearest points to a query location
std::vector<size_t> nearest = finder.kNearest({1.1, 0, 0}, 2);

// 2 nearest neighbors of points[1], excluding itself
std::vector<size_t> neighbors = finder.kNearestNeighbors(1, 2);

// all points within radius 1.5 of a query location
std::vector<size_t> inBall = finder.radiusSearch({0, 0, 0}, 1.5);
```

All query results are indices into the input `points` vector.

??? func "`#!cpp NearestNeighborFinder::NearestNeighborFinder(const std::vector<Vector3>& points)`"

    Construct a finder over the given point set and build the spatial index. All subsequent queries return indices into this vector.

### Queries

??? func "`#!cpp std::vector<size_t> NearestNeighborFinder::kNearest(Vector3 query, size_t k)`"

    Returns the indices of the `k` closest points in the input set to `query`, in order of increasing distance. `query` need not be one of the input points.

??? func "`#!cpp std::vector<size_t> NearestNeighborFinder::kNearestNeighbors(size_t sourceInd, size_t k)`"

    Returns the indices of the `k` closest points to `points[sourceInd]`, excluding `sourceInd` itself. Useful for computing point-to-point neighbor relationships within the input set.

??? func "`#!cpp std::vector<size_t> NearestNeighborFinder::radiusSearch(Vector3 query, double rad)`"

    Returns the indices of all points within distance `rad` of `query`.
    
    Pass the actual desired radius, not the squared radius that some other libraries use.
