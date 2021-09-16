# Integer Coordinates Intrinsic Triangulation

![integer coordinates example](/media/integer_coordinates.png){: style="max-width: 20em; display: block; margin-left: auto; margin-right: auto;"}

The integer coordinates data structure encodes an intrinsic triangulation using _normal coordinates_ and _roundabouts_, which are integer values stored on edges. Since it avoids floating point data, this data structure is more robust than signposts, but operations such as vertex insertion can take longer.
For more details, see [Integer Coordinates for Intrinsic Geometry Processing](https://arxiv.org/pdf/2106.00220.pdf).

The `IntegerCoordinatesIntrinsicTriangulation` is one concrete implementation (subclass) of [intrinsic triangulations](../basics). **Most of its functionality is documented there.**

### Tradeoffs

The signpost intrinsic triangulation is one of the two main intrinsic triangulation representations currently available in geometry-central, the other being [signposts](../signposts). For many tasks, either reprentation will be highly effective, but there are some tradeoffs.

**Pros:**

- **Robustness.** Integer coordinates are guaranteed to always maintain a topologically valid representation of the common subdivision, through any sequence of operations.

**Cons:**

- **Performance.** Integer coordinates may be moderately more expensive in terms of runtime than signposts (although both are often on the order of milliseconds). Some operations, such as tracing out a single edge of the intrinsic triangulation, are only implemented by first extracting the entire common subdivision, which has some overhead.

!!! warning "Not yet implemented"

    The methods `equivalentPointOnIntrinsic()` and `splitEdge()` from `IntrinsicTriangulation` are not yet implemented for the integer coordinates representation. For `splitEdge()`, an alternate version is provided which returns a vertex instead.


### Constructors

??? func "`#!cpp IntegerCoordinatesIntrinsicTriangulation::IntegerCoordinatesIntrinsicTriangulation(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& inputGeom, double mollifyEPS=1e-5)`"

    Initialize an intrinsic triangulation sitting on top of `mesh`. Recall that `IntrinsicGeometryInterface` can be almost any geometry object, including a `VertexPositionGeometry`.

    Initially, the intrinsic triangulation will be identical to the input mesh; it can be modified with the routines below.
   
    The `mollifyEPS` parameter performs initial [mollification](/surface/algorithms/robust_geometry/#intrinsic-mollification) on the intrinsic triangulation, which greatly improves floating-point robustness while generally, while having a negligible impact on accuracy. Set `mollifyEPS=0` to disable.

## Citation

```bib
@article{gillespie2021integer,
  title={Integer Coordinates for Intrinsic Geometry Processing},
  author={Gillespie, Mark and Sharp, Nicholas and Crane, Keenan},
  journal={arXiv preprint arXiv:2106.00220},
  year={2021}
}
```
