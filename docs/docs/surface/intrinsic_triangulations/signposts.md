# Signpost Intrinsic Triangulation

![signpost example](/media/signposts.svg){: style="max-width: 25em; display: block; margin-left: auto; margin-right: auto;"}

The signpost data structure encodes an intrinsic triangulation by storing "signposts" at mesh vertices. Explicitly, for each intrinsic edge it stores the edge's length and direction at the two incident vertices. This information fully specifies how the intrinsic triangulation sits above the input mesh. For more details, see [Navigating Intrinsic Triangulations](http://www.cs.cmu.edu/~kmcrane/Projects/NavigatingIntrinsicTriangulations/paper.pdf).

The `SignpostIntrinsicTriangulation` is one concrete implementation (subclass) of [intrinsic triangulations](../basics). **Most of its functionality is documented there.** Here, we document additional routines which are specific to signposts.

### Tradeoffs

The signpost intrinsic triangulation is one of the two main intrinsic triangulation representations currently available in geometry-central, the other being [integer coordinates](../integer_coordinates). For many tasks, either reprentation will be highly effective, but there are some tradeoffs.

**Pros:**

- **Performance.** Signposts are moderately faster in terms of runtime than integer coordinates (although both are often on the order of milliseconds).

- **Tangent vector data.** Signposts naturally offer tangent space coordinate systems which are consistent with the input mesh, making it easy to work with tangent-valued data at vertices of the intrinsic triangulation.

**Cons:**

- **Robustness.** The representation relies heavily on tracing intrinsic edge paths across the input surface, which can be error-prone in floating-point arithmetic. In particular, reconstructing the common subdivision may fail if the input contains degenerate triangles.


### Constructors

??? func "`#!cpp SignpostIntrinsicTriangulation::SignpostIntrinsicTriangulation(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& inputGeom)`"

    Initialize an intrinsic triangulation sitting on top of `mesh`. Recall that `IntrinsicGeometryInterface` can be almost any geometry object, including a `VertexPositionGeometry`.

    Initially, the intrinsic triangulation will be identical to the input mesh; it can be modified with the routines below.

### Methods

??? func "`#!cpp std::vector<SurfacePoint> SignpostIntrinsicTriangulation::traceIntrinsicHalfedgeAlongInput(Halfedge intrinsicHe, bool trimEnd)`"

    This function is generally the same as `traceIntrinsicHalfedgeAlongInput()`.

    When edges paths from signposts, the path often does not _exactly_ hit the destination vertex, but rather ends somewhere very close in the adjacent 1-ring. If `trimEnd=true`, a simple heuristic is used to clean up the path so it exactly hits the target vertex; with `trimEnd=false` the result of tracing is directly reported.

    By default `traceIntrinsicHalfedgeAlongInput(he)` is equivalent to `traceIntrinsicHalfedgeAlongInput(he, true)`.


## Citation

```bib
@article{sharp2019navigating,
  title={Navigating intrinsic triangulations},
  author={Sharp, Nicholas and Soliman, Yousuf and Crane, Keenan},
  journal={ACM Transactions on Graphics (TOG)},
  volume={38},
  number={4},
  pages={55},
  year={2019},
  publisher={ACM}
}
```
