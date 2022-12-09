
| Original | Remeshed |
| - | - |
|![finely_triangulated_spot](/media/spot-triangulation-original.png){: style="max-width: 20em; display: block; margin-left: auto; margin-right: auto;"}|![coarsely_triangulated_spot](/media/spot-triangulation-remeshed.png){: style="max-width: 20em; display: block; margin-left: auto; margin-right: auto;"}


??? func "`#!cpp void remesh(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom, RemeshOptions options = defaultRemeshOptions)`"

    Remesh a mesh using vertex smoothing along with edge splits, collapses, and flips. Options are passed as a [RemeshOptions](#options) object.
    
??? func "`#!cpp void remesh(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom, MutationManager& mm, RemeshOptions options = defaultRemeshOptions)`"

    Remesh a mesh using vertex smoothing along with edge splits, collapses, and flips. Options are passed as a [RemeshOptions](#options) object.
    All mesh mutations are performed through a `MutationManager` object, which can e.g. keep special vertices fixed or run callback functions when certain mesh mutations occur.

## Vertex Smoothing
Vertex 

## Boundary Conditions

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


