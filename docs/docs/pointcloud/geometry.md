The geometry classes associates 3D positions with the points in a point cloud, and can compute various derived quantities.  These use the same paradigm of "requiring quantities" as the [surface mesh geometry interface](/surface/geometry/geometry), see additional documentation there.  


```#include "geometrycentral/pointcloud/point_position_geometry.h"```

!!! warning "Point clouds are in beta"

    The current point cloud API in geometry-central is preliminary, and may change in future versions.

The geometry classes manage buffers of geometric quantities which can be computed, as well as a dependency graph relationship between these quantities. For instance, calling `geom.requireNormals()` will ensure the `geom.normals` buffer is populated with normals, additionally computing & caching any preceding quantities which are needed (here: neighborhoods).  Derived classes can be used to extend with more features, specify additional input data, or override the computation of particular quantities.

**Example:** Basic use of quantities
```c++
#include "geometrycentral/pointcloud/point_cloud.h"
#include "geometrycentral/pointcloud/point_position_geometry.h"

// create a new cloud & geometry object
size_t nPts = 5000;
PointCloud cloud(nPts); 
PointData<Vector3> positions = /* some positions */;
PointPositionGeometry geom(*cloud, positions);

// Per-point normals
geom.requireNormals();
for(Point p : cloud.points()) {
  std::cout << "normal for point " << p << " is " << geom.normals[p] << "\n";
}

// Point cloud Laplace matrix 
geom.requireLaplacian();
Eigen::SparseMatrix<double> L = geom.laplacian;
```



## Point Position Geometry

This is the most basic core geometry, corresponding to a 3D position associated with each vertex.

Remember, for each quantity YYY, call `requireYYY()` to ensure it has been computed, then access it at `geom.YYYs`.

### Construction

??? func "`#!cpp PointPositionGeometry::PointPositionGeometry(PointCloud& cloud)`"

    Constructs a new unitialized point position geometry on the cloud.

    You should immediately set positions like

    ```cpp
    PointPositionGeometry goem(*cloud);
    for(Point p : cloud->points()) {
      geom.positions[p] = /* some Vector3 position */;
    }
    ```

??? func "`#!cpp PointPositionGeometry::PointPositionGeometry(PointCloud& cloud, const PointData<Vector3>& positions)`"

    Constructs a new point position geometry with specified positions.

    The given positions will be copied to the `geom.positions` field.

### Quantities

??? func "neighbors"
    
    ##### neighbors 

    A collection of `k`-nearest neighbors for each point. Represented via a `Neighborhoods` object, which encapsulates the functionality.

    Set the number of neighbors by assigning to `PointPositionGeometry::kNeighborSize`. The default is 30. If you're going to change it, do so _before_ requiring any quantities that involve neighbors (most of them).

    Iterate over neighbors like:
    
    ```cpp
    geom.requireNeighbors();
    for (Point p : cloud->points()) {
      std::vector<Point>& neigh = geom.neighbors->neighbors[p];
      size_t M = neigh.size();
      for (size_t iN = 0; iN < M; iN++) {
        Point pN = neigh[iN];
        std::cout << "Point " << p << " has neighbor " << pN << std::endl;
      }
    }
    ```

    - **member:** `std::unique_ptr<Neighborhoods> PointPositionGeometry::neighbors`
    - **require:** `void PointPositionGeometry::requireNeighbors()`

    !!! note "Changes coming"

        A future update will replace this nested `std::vector<>` interface with a flat list and fancy iterators.

??? func "normals"
    
    ##### normals

    A 3D unit normal at each point.

    Normals are comptuted via PCA over neighbors. By default, nothing is done to orient normals, so orientations will be arbitrary.

    - **member:** `PointData<Vector3> PointPositionGeometry::normals`
    - **require:** `void PointPositionGeometry::requireNormals()`

??? func "tangent basis"
    
    ##### tangent basis

    A pair of X-Y tangent basis axes at each point. Guaranteed to form an orthonormal right-handed coordinate frame with the normal vector.

    - **member:** `PointData<std::array<Vector3,2>> PointPositionGeometry::tangentBasis`
    - **require:** `void PointPositionGeometry::requireTangentBasis()`


??? func "neighborhood tangent coordinates"
    
    ##### tangent coordinates

    Local 2D tangent coordinates associated with each neighboring point, corresponding to projection in to the axes in `geom.tangentBasis`.

    - **member:** `PointData<std::vector<Vector2>> PointPositionGeometry::tangentCoordinates`
    - **require:** `void PointPositionGeometry::requireTangentCoordinates()`

??? func "neighborhood tangent transport"
    
    ##### tangent transport

    Parallel transport coefficients to rotate tangent vectors between neighboring frames.`tangentTransport[i][j]` holds the rotation which maps a vector in the tangent space of i to that of j.

    - **member:** `PointData<std::vector<Vector2>> PointPositionGeometry::tangentTransport`
    - **require:** `void PointPositionGeometry::requireTangentTransport()`

??? func "tufted triangulation"
    
    ##### tufted triangulation

    A tufted intrinsic triangulation associated with the point cloud. Intuitively this is a special triangulation atop the points in the pointcloud, represented only by its connectivity and edge lengths. It is a very effective numerical data structure to compute downstream quantities, like a highly-quality Laplace matrix. To be clear, this is _not_ a "nice" triangulation like you might get from 3D reconstruction; instead it is a crazy nonmanifold triangulation which happens to have a very useful structure for subsequent numerical computations.  See the publication ["A Laplacian for Nonmanifold Triangle Meshes"](http://www.cs.cmu.edu/~kmcrane/Projects/NonmanifoldLaplace/NonmanifoldLaplace.pdf) for formal details.

    Intuitively, 

    - **member:** `PointData<std::unique_ptr<surface::SurfaceMesh>> PointPositionGeometry::tuftedMesh`
    - **member:** `PointData<std::unique_ptr<surface::EdgeLengthGeometry>> PointPositionGeometry::tuftedGeom`
    - **require:** `void PointPositionGeometry::requireTuftedTriangulation()`

??? func "Laplacian"
    
    ##### Laplacian

    A Laplace matrix for the point cloud. Computed internally using the tufted triangulation as described in [A Laplacian for Nonmanifold Triangle Meshes"](http://www.cs.cmu.edu/~kmcrane/Projects/NonmanifoldLaplace/NonmanifoldLaplace.pdf).

    - **member:** `Eigen::SparseMatrix<double> PointPositionGeometry::laplacian`
    - **require:** `void PointPositionGeometry::requireLaplacian()`

??? func "connection Laplacian"
    
    ##### connection Laplacian

    A connection Laplace matrix for the point cloud, which is similar to normal scalar Laplacian, but operates on tangent vectors at each point rather than on scalar values. See ["The Vector Heat Method"](https://nmwsharp.com/research/vector-heat-method) for a detailed introduction to connection Laplacians.

    **Note:** This connection Laplacian includes an orientation flip to handle inconsistent normals. Rather than a complex matrix, it builds a `2N x 2N` expanded real matrix, with sign flips on imaginary elements to conjugate (a necessary trick, because conjugation is not a complex-linear operation). This allows us to apply the connection even when normal orientations are inconsistent.

    - **member:** `Eigen::SparseMatrix<double> PointPositionGeometry::connectionLaplacian`
    - **require:** `void PointPositionGeometry::requireConnectionLaplacian()`

??? func "gradient"
    
    ##### gradient

    A gradient operator `G` for the point cloud. Given a vector `v` of scalar values at each point in the point cloud, multiplying `G*v` yields a tanget vector at each point (as a complex number), representing the spatial gradient of `v` at that point.

    Constructed via a least-squares approximation in the neighborhood tangent coordinate space at each point.

    - **member:** `Eigen::SparseMatrix<std::complex<double>> PointPositionGeometry::gradient`
    - **require:** `void PointPositionGeometry::requireGradient()`


## Point Position & Normal Geometry

The class `PointPositionNormalGeometry` extends the basic `PointPositionGeometry` class to manually specify a set of known normals at points. These normals will be used for all computations, rather than computing normals from scratch.

The `PointPositionNormalGeometry` **is-a** `PointPositionGeometry` via inheritance, so you can directly pass it anywhere a `PointPositionGeometry` is expected.

```#include "geometrycentral/pointcloud/point_position_normal_geometry.h"```

**Example:**
```cpp
std::unique_ptr<PointCloud> cloud; // your cloud
PointData<Vector3> positions(*cloud); // populate with your positions
PointData<Vector3> knownNormals(*cloud); // populate with your normals

PointPositionNormalGeometry geom(*cloud, positions, knownNormals);

// we can now use this just like a `PointPositionGeometry` as above
```

??? func "`#!cpp PointPositionNormalGeometry::PointPositionNormalGeometry(PointCloud& cloud, const PointData<Vector3>& positions, const PointData<Vector3>& normals)`"

    Constructs a new point position geometry with specified normals.

    The given normals will be copied to the `geom.normals` field.


## Point Position & Tangent Basis Geometry

The class `PointPositionFrameGeometry` extends the basic `PointPositionGeometry` class to manually specify a set of known normals and tangent frames at points. These frames / normals will be used for all computations, rather than computing frames / normals from scratch.

The input frames are specified via three vectors at each point `std::array<Vector3, 3>`. These vectors should form a right-handed, orthonormal basis at each vertex, with the two vectors forming the X/Y tangent basis, and the third as the normal direction.

The `PointPositionFrameGeometry` **is-a** `PointPositionGeometry` via inheritance, so you can directly pass it anywhere a `PointPositionGeometry` is expected.

```#include "geometrycentral/pointcloud/point_position_frame_geometry.h"```

**Example:**
```cpp
std::unique_ptr<PointCloud> cloud; // your cloud
PointData<Vector3> positions(*cloud); // populate with your positions
// your {tangent-X, tangent-Y, normal} at each vertex
PointData<std::array<Vector3, 3>> knownFrames(*cloud); 

PointPositionFrameGeometrygeom(*cloud, positions, knownNormals);

// we can now use this just like a `PointPositionGeometry` as above
```

??? func "`#!cpp PointPositionFrameGeometry::PointPositionFrameGeometry(PointCloud& cloud, const PointData<Vector3>& positions, const PointData<std::array<Vector3, 3>>& frames)`"

    Constructs a new point position geometry with specified normals.

    The given tangent bases and normals will be copied to the `geom.normals` and `geom.tangentBasis` fields.
