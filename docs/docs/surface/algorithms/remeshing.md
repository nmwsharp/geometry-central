
| Original | Remeshed |
| - | - |
|![finely_triangulated_spot](/media/spot-triangulation-original.png){: style="max-width: 20em; display: block; margin-left: auto; margin-right: auto;"}|![coarsely_triangulated_spot](/media/spot-triangulation-remeshed.png){: style="max-width: 20em; display: block; margin-left: auto; margin-right: auto;"}

These routines try to improve mesh quality using repeated rounds of [vertex position smoothing](#tangential-vertex-smoothing), [edge flipping](#extrinsic-delaunay-flipping), and [edge splits/collapses](#edge-length-adjustment).

??? func "`#!cpp void remesh(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom, RemeshOptions options = defaultRemeshOptions);`"

    Remesh a mesh using vertex smoothing along with edge splits, collapses, and flips. Options are passed as a [RemeshOptions](#options) object.
    
??? func "`#!cpp void remesh(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom, MutationManager& mm, RemeshOptions options = defaultRemeshOptions);`"

    Remesh a mesh using vertex smoothing along with edge splits, collapses, and flips. Options are passed as a [RemeshOptions](#options) object.
    All mesh mutations are performed through a `MutationManager` object, which can e.g. keep special vertices fixed or run callback functions when certain mesh mutations occur.

## Tangential Vertex Smoothing
Vertex smoothing moves each vertex towards the average of its neighborhood. This average can be computed in two ways: in *circumentric* smoothing, every vertex is moved towards the area-weighted average of the circumcenters of its neighboring faces (as in [[Chen & Holst 2011]](https://www.sciencedirect.com/science/article/abs/pii/S0045782510003117)), whereas in *Laplacian* smoothing, every vertex is moved towards the average of the positions of its neighboring vertices.

In either case we perform *tangential smoothing*, meaning that we only move vertices tangentially along the surface. 

??? func "`#!cpp double smoothByCircumcenter(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, double stepSize = 1, RemeshBoundaryCondition bc = RemeshBoundaryCondition::Tangential);`"
                         
    Move each vertex tangentially towards the area-weighted average of its neighboring face circumcenters.
    Returns the average amount that each vertex was moved by.
    
??? func "`#!cpp double smoothByCircumcenter(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, MutationManager& mm, double stepSize = 1, RemeshBoundaryCondition bc = RemeshBoundaryCondition::Tangential);`"
    
    Move each vertex tangentially towards the area-weighted average of its neighboring face circumcenters
    , using the provided `MutationManager` to move the vertices
    Returns the average amount that each vertex was moved by.

??? func "`#!cpp double smoothByLaplacian(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, double stepSize = 1, RemeshBoundaryCondition bc = RemeshBoundaryCondition::Tangential);`"
                         
    Move each vertex tangentially towards the average of its neighbors
    Returns the average amount that each vertex was moved by.
??? func "`#!cpp double smoothByLaplacian(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, MutationManager& mm, double stepSize = 1, RemeshBoundaryCondition bc = RemeshBoundaryCondition::Tangential);`"
    
    Move each vertex tangentially towards the average of its neighbors, using the provided `MutationManager` to move the vertices
    Returns the average amount that each vertex was moved by.

### Boundary Conditions
The boundary conditions for smoothing can be set to `RemeshBoundaryCondition::Tangential`, `RemeshBoundaryCondition::Fixed`, or `RemeshBoundaryCondition::Free`. Tangential smoothing allows boundary vertices to move along the tangent direction to the boundary curve, fixed smoothing fixes the boundary vertices, and free smoothing allows boundary vertices to move along any direction in the surface's tangent plane. Leaving boundary vertices free allows the mesh to degenerate, so it is not recommended unless you impose additional constraints.

## Extrinsic Delaunay Flipping
Delaunay flipping performs edge flips to improve triangle quality.

!!! warning "No guarantees"

    Unlike the [intrinsic Delaunay flipping routines](/surface/intrinsic_triangulations/basics/#high-level-mutators), extrinsic flipping algorithms are not guaranteed to produce a Delaunay mesh. Nonetheless, these edge flips generally improve triangle quality in practice.
    
??? func "`#!cpp size_t fixDelaunay(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom);`"

    Try to make all triangles Delaunay using extrinsic edge flips.
    Returns the number of flips performed.
    
??? func "`#!cpp size_t fixDelaunay(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, MutationManager& mm);`"
    Try to make all triangles Delaunay using extrinsic edge flips, using the provided `MutationManager` to flip edges.
    Returns the number of flips performed.
    
## Edge Length Adjustment
These routines perform edge splits and collapses to drive mesh edge lengths towards a target value. If curvature adaptation is enabled, this target length is made shorter in high-curvature areas, leading to more resolution there.

??? func "`#!cpp bool adjustEdgeLengths(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, RemeshOptions options = defaultRemeshOptions);`"

    Apply splits and collapses to adjust edge lengths.
   
    Reads the following options from `options`, with the following default values and meanings:
    
    * `#!cpp double options.targetEdgeLength = -1`: the target edge length in flat regions. If `targetEdgeLength` is negative, the target length is set relative to the input mesh's mean length.
    * `#!cpp double options.curvatureAdaptation = 0`: how much target edge length should vary due to curvature. Set `curvatureAdaptation` to zero if you want edge lengths to be approximately `targetEdgeLength` everywhere.
    * `#!cpp double options.minRelativeLength = 0.05`: the minimum possible edge length in the output mesh. Defined relative to `targetEdgeLength`.

??? func "`#!cpp bool adjustEdgeLengths(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, MutationManager& mm, RemeshOptions options = defaultRemeshOptions);`"

    Apply splits and collapses to adjust edge lengths.
    All splits and collapses are performed using the provided `MutationManager`.
   
    Reads the following options from `options`, with the following default values and meanings:
    
    * `#!cpp double options.targetEdgeLength = -1`: the target edge length in flat regions. If `targetEdgeLength` is negative, the target length is set relative to the input mesh's mean length.
    * `#!cpp double options.curvatureAdaptation = 0`: how much target edge length should vary due to curvature. Set `curvatureAdaptation` to zero if you want edge lengths to be approximately `targetEdgeLength` everywhere.
    * `#!cpp double options.minRelativeLength = 0.05`: the minimum possible edge length in the output mesh. Defined relative to `targetEdgeLength`.

## Helper Types
### Options
Options are passed in to `remesh` via a `RemeshOptions` object.

| Field                                              | Default value                         | Meaning                                                                                                                                                    |
|----------------------------------------------------|---------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `#!cpp double targetEdgeLength`                   | `-1`                                  | the target edge length in flat regions. If `targetEdgeLength` is negative, the target edge length is set to relative the input mesh's mean edge length     |
| `#!cpp size_t maxIterations`                      | `10`                                  | the maximum number of iterations to run for                                                                                                                |
| `#!cpp double curvatureAdaptation`                | `0`                                   | how much target length should vary due to curvature. Set curvatureAdaptation to 0 if you want edge lengths to be approximately targetEdgeLength everywhere |
| `#!cpp double minRelativeLength`                  | `0.05`                                | the minimum possible edge length allowed in the output mesh. Defined relative to targetEdgeLength                                                          |
| `#!cpp RemeshSmoothStyle smoothStyle`             | `RemeshSmoothStyle::Circumcentric`    | the type of vertex smoothing to use (either `RemeshSmoothStyle::Circumcentric` or `RemeshSmoothStyle::Laplacian`)                                          |
| `#!cpp RemeshBoundaryCondition boundaryCondition` | `RemeshBoundaryCondition::Tangential` | the type of motions allowed for boundary vertices (either `RemeshBoundaryCondition::Fixed`, `RemeshBoundaryCondition::Tangential` or `RemeshBoundaryCondition::Free`)                            |

!!! warning "'Fixed' boundary may still move slightly"

    Setting `boundaryCondition` to `RemeshBoundaryCondition::Fixed` only fixes boundary vertices during vertex smoothing. Boundary edges may still be split or collapsed during edge length adjustment, which can also cause the boundary to move slightly. To stop all motion of the boundary, you can pass in a `MutationManager` which marks all boundary edges as not splittable or collapsible.
