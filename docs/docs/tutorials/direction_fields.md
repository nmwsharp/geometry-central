This tutorial shows a higher-level example in geometry-central, where we use a built-in routine to visualize a direction field.

[View full, runnable source code in the tutorial repository.](https://github.com/nmwsharp/geometry-central-tutorials)


### Basic setup

To begin we include headers, bring in namespaces, See the first tutorial for more info on these steps.

```cpp
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/direction_fields.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;
```

Again, load an input mesh from file, and visualize it using Polyscope. We'll require that the mesh be manifold so we have obvious vertex tangent spaces.
  
```cpp 
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;
std::tie(mesh, geometry) = readManifoldSurfaceMesh(args::get(inputFilename));

polyscope::init();
auto* psMesh =
    polyscope::registerSurfaceMesh("input mesh",
        geometry->vertexPositions, 
        mesh->getFaceVertexList());
```

Now, computing the direction field itself amoutns to a single function call
```cpp
VertexData<Vector2> field = computeSmoothestVertexDirectionField(*geometry);
```


Before we visualize this direction field, we need to do a little more setup. The tangent vectors are defined in an abritrary, intrinsic tangent space at each vertex. We need to tell Polyscope how these tangent spaces sit in 3D space.

```cpp
geometry->requireVertexTangentBasis();
VertexData<Vector3> vBasisX(*mesh);
for (Vertex v : mesh->vertices()) {
  vBasisX[v] = geometry->vertexTangentBasis[v][0];
}
psMesh->setVertexTangentBasisX(vBasisX);
```

Now we can add the vector field to Polyscope and inspect it!
```cpp
psMesh->addVertexIntrinsicVectorQuantity("vectors", field);
polyscope::show();
```


Result:
![direction field](/media/tutorials/dir_field.jpeg)
