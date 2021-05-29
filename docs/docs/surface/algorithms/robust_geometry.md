Many data structures and algorithms in geometry-central help to "robustify" algorithms, and help them automatically work better on meshes with poor geometry and/or connectivity. Whenever possible, this done in a way which is transparent to the algorithm designer and end user.


### Intrinsic Mollifiction

`#include "geometrycentral/surface/intrinsic_mollification.h"`

Poor quality meshes might have triangles whose shape is so close to being degenerate that even basic floating point arithmetic breaks down. Such features are difficult to resolve with by perturbing vertex positions, because a motion that makes one triangle better might make another triangle worse. Fortunately, when working with intrinsic geometry, there is a simple strategy which is always guaranteed to work.

To ensure that all triangles are nondegenerate, we want to have

$$
l_{ij} + l_{jk} > l_{ki} + \delta
$$

for all triples of triangle edge lengths, according to some user-specified numerical tolerance $\delta$, which specified how far from degenerate the triangles should be.  This is easily achieved by simply adding some small value $\epsilon$ to all edge lengths, chosen to be the smallest $\epsilon$ which will make the above inequality hold. 

This strategy ensures all triangles are non-degenerate, yet yields a negligible change to the mesh's geometry. For a mesh which already has no degenerate triangles, it will have no effect.

Example:

```cpp
#include "geometrycentral/surface/intrinsic_mollification.h"
using namespace geometrycentral;
using namespace surface;

// your mesh and intrinsic geometry
std::unique_ptr<SurfaceMesh> mesh;
std::unique_ptr<EdgeLengthGeometry> geometry;

// mollify the edge lengths
mollifyIntrinsic(*mesh, geometry->edgeLengths);

// ensure that any existing quantities are updated for the new edge lengths
geometry->refreshQuantities();

// continue running algorithms, etc
geometry->requireVertexGaussianCurvatures();
for(Vertex v : mesh->vertices()) {
  std::cout << "Gaussian curvature of " << v << " is " 
    << geometry->vertexGaussianCurvatures[v] << std::endl;
}
```



??? func "`#!cpp void mollifyIntrinsic(SurfaceMesh& mesh, EdgeData<double>& edgeLengths, double relativeFactor = 1e-6)`"
 
    Mollify the edge lengths in `edgeLengths`, to ensure that all triangles are nondegenerate by a factor $\delta$, where $\delta$ is computed as `relativeFactor` times the mean edge length in the input.


??? func "`#!cpp void mollifyIntrinsicAbsolute(SurfaceMesh& mesh, EdgeData<double>& edgeLengths, double absoluteFactor)`"

    Similar to above, but mollifies with an absolute factor given by `absoluteFactor`, rather than a relative factor.


This strategy is described in [A Laplacian for Nonmanifold Triangle Meshes](http://www.cs.cmu.edu/~kmcrane/Projects/NonmanifoldLaplace/NonmanifoldLaplace.pdf). The appropriate citation is:

```bib
@article{Sharp:2020:LNT,
  author={Nicholas Sharp and Keenan Crane},
  title={{A Laplacian for Nonmanifold Triangle Meshes}},
  journal={Computer Graphics Forum (SGP)},
  volume={39},
  number={5},
  year={2020}
}
```
