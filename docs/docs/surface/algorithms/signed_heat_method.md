This section describes the _Signed Heat Method_ in geometry-central, which computes signed and unsigned distance to possibly broken geometry using heat flow.

Note that these quantities all depend on the _intrinsic_ geometry of a surface (via the `IntrinsicGeometryInterface`). Therefore, these routines can be run on abstract geometric domains as well as traditional surfaces in 3D.

These algorithms are described in [A Heat Method for Generalized Signed Distance](https://nzfeng.github.io/research/SignedHeatMethod/SignedDistance.pdf). 

`#include "geometrycentral/surface/signed_heat_method.h"`


## Signed Heat Solver

The stateful class `SignedHeatSolver` shares precomputation for all of the routines below. What this means is that later solves may become significantly faster, as precomputation occurs on earlier solves.

??? func "`#!cpp SignedHeatSolver::SignedHeatSolver(IntrinsicGeometryInterface& geom, double tCoef=1.0)`"

    Create a new solver for the Signed Heat Method. Precomputation is performed lazily as needed.

    - `geom` is the geometry (and hence mesh) on which to compute. Note that nearly any geometry object (`VertexPositionGeometry`, etc) can be passed here.

    - `tCoef` is the time to use for short time heat flow `tCoef * h^2`, where `h` is the mean distance between nodes of the mesh. The default value of `1.0` is almost always sufficient.


## Signed & Unsigned Geodesic Distance

![generalized signed distance on a panda and armchair](/media/signed_heat_method.png)

Example:
```cpp
#include "geometrycentral/surface/signed_heat_method.h"

// your mesh and geometry
VertexPositionGeometry geometry;
SurfaceMesh mesh;

// construct a solver
SignedHeatSolver signedHeatSolver(geometry);

// specify some source geometry
std::vector<Curve> curves;
curves.emplace_back();
curves.back().nodes.emplace_back(mesh.vertex(0));
curves.back().nodes.emplace_back(mesh.edge(5), 0.3);

// solve!
VertexData<double> distance = signedHeatSolver->computeDistance(curves);
```

??? func "`#!cpp VertexData<double> SignedHeatSolver::computeDistance(const std::vector<Curve>& curves, const std::vector<SurfacePoint>& points, const SignedHeatOptions& options = SignedHeatOptions())`"

    Compute distance to a collection of curves and points. Curves may be a source of either signed or unsigned distance (see `Curve` object below.) Points are given as a list of [surface points](../../utilities/surface_point/) to which we compute unsigned distance.

??? func "`#!cpp VertexData<double> SignedHeatSolver::computeDistance(const std::vector<Curve>& curves, const SignedHeatOptions& options = SignedHeatOptions())`"

    Compute distance to a collection of curves. Curves may be a source of either signed or unsigned distance (see `Curve` object below.) 

??? func "`#!cpp VertexData<double> SignedHeatSolver::computeDistance(const std::vector<SurfacePoint>& points, const SignedHeatOptions& options = SignedHeatOptions())`"

    Compute unsigned distance to a collection of isolated point sources. Point sources are given as a list of [surface points](../../utilities/surface_point/).

The diffusion time can also be changed using the following function.

??? func "`#!cpp void SignedHeatSolver::setDiffusionTimeCoefficient(double tCoef)`"

    Re-computes the time used for short time heat flow `tCoef * h^2`, where `h` is the mean distance between nodes of the mesh.

## Helper Types

### Curves

| Field | Default value |Meaning|
|---|---|---|
| `#!cpp std::vector<SurfacePoint> nodes`| `std::vector<SurfacePoint>()` | The nodes of the curve, given as an ordered sequence of [surface points](../../utilities/surface_point/). |
| `#!cpp bool isSigned`| `true` | Whether the curve is oriented or not. If `isSigned` is `true`, then (generalized) signed distance will be computed to the curve; unsigned distance otherwise. |

### Options
Options are passed in to `computeDistance` via a `SignedHeatOptions` struct, which has the following fields.

| Field | Default value |Meaning|
|---|---|---|
| `#!cpp bool preserveSourceNormals`| `false` | If `true`, preserve the initial curve normals at the source curve during vector diffusion. |
| `#!cpp LevelSetConstraint levelSetConstraint`| `LevelSetConstraint::ZeroSet` | Specifies how/if level sets should be preserved. Can be set to `LevelSetConstraint::ZeroSet`, `LevelSetConstraint::Multiple`, or `LevelSetConstraint::None`, corresponding to preserving the zero set, mulitple level sets (one for each curve component), or no level sets, respectively. |
| `#!cpp double softLevelSetWeight`| `-1` | If greater than 0, gives the weight with which the given level set constraint is "softly" enforced. |

## Citation

If these algorithms contribute to academic work, please cite the following paper:

```bib
@article{Feng:2024:SHM,
  author = {Feng, Nicole and Crane, Keenan},
  title = {A Heat Method for Generalized Signed Distance},
  year = {2024},
  issue_date = {August 2024},
  publisher = {Association for Computing Machinery},
  address = {New York, NY, USA},
  volume = {43},
  number = {4},
  issn = {0730-0301},
  url = {https://doi.org/10.1145/3658220},
  doi = {10.1145/3658220},
  journal = {ACM Trans. Graph.},
  month = {jul},
  articleno = {92},
  numpages = {19}
}
```
