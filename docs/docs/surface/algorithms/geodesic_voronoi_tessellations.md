A Voronoi tessellation partitions a domain $M$ based on a set of sites $s_i \in M$. To each site $s_i$, we associate a region (or "Voronoi cell") consisting of the points closer to $s_i$ than to any other site. Such a tessellation is called _centroidal_ if each site is at the center of its corresponding cell. This routine computes centroidal Voronoi tessellations on surface meshes.

!!! warning "Output is not unique!"

    By default, this procedure picks some random points and optimizes them to compute a centroidal Voronoi tessellation. Since the initialization is random, the results will generally be different each time the procedure is run.

`#include "geometrycentral/surface/geodesic_centroidal_voronoi_tessellation.h"`

![karcher means](/media/geodesic-centroidal-voronoi.png){: style="max-width: 25em; display: block; margin-left: auto; margin-right: auto;"}


??? func "`#!cpp VoronoiResult computeGeodesicCentroidalVoronoiTessellation(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom, VoronoiOptions options = defaultVoronoiOptions)`"

    Compute a geodesic centroidal Voronoi tessellation on the input mesh. Options are passed as a [VoronoiOptions](#options) object, and the output is returned as a [VoronoiResult](#result) object.

## Helper Types
### Options
Options are passed in to `computeGeodesicCentroidalVoronoiTessellation` via a `VoronoiOptions` object.

| Field | Default value |Meaning|
|---|---|---|
| `#!cpp size_t nSites;`| `10` | the number of sites to place |
| `#!cpp std::vector<SurfacePoint> initialSites;`| `{}` | desired locations for sites. If blank, initial locations are chosen randomly |
| `#!cpp size_t iterations;`| `50` | number of iterations to run for |
| `#!cpp bool useDelaunay;`| `true` | solve on an [intrinsic Delaunay triangulation](/surface/intrinsic_triangulations/basics) of the input |
| `#!cpp double tCoef;`| `1` | diffusion time for the [vector heat method](/surface/algorithms/vector_heat_method) |
| `#!cpp size_t nSubIterations;`| `1` | number of iterations to use when computing [surface centers](/surface/algorithms/surface_centers) during optimization |


### Result
The result is returned as a `VoronoiResult`, which has 3 fields:

| Field | Meaning|
|---|---|
| `#!cpp std::vector<SurfacePoint> siteLocations;`| the sites at the centers of the Voronoi cells |
| `#!cpp std::vector<VertexData<double>> siteDistributions;`| indicator functions which each Voronoi cell |
| `#!cpp bool hasDistributions;`| is `siteDistributions` populated? |

## Citation

This algorithms is described in [The Vector Heat Method](http://www.cs.cmu.edu/~kmcrane/Projects/VectorHeatMethod/paper.pdf), the appropriate citation is:

```bib
@article{sharp2019vector,
  title={The Vector Heat Method},
  author={Sharp, Nicholas and Soliman, Yousuf and Crane, Keenan},
  journal={ACM Transactions on Graphics (TOG)},
  volume={38},
  number={3},
  pages={24},
  year={2019},
  publisher={ACM}
}
```
