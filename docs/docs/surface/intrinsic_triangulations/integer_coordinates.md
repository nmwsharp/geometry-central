# Integer Coordinates Intrinsic Triangulation

![integer coordinates example](/media/integer_coordinates.png){: style="max-width: 20em; display: block; margin-left: auto; margin-right: auto;"}

The integer coordinates data structure encodes an intrinsic triangulation using _normal coordinates_ and _roundabouts_, which are integer values stored on edges. Since it avoids floating point data, this data structure is more robust than signposts, but operations such as vertex insertion can take longer.
For more details, see [Integer Coordinates for Intrinsic Geometry Processing](https://arxiv.org/pdf/2106.00220.pdf).

The `IntegerCoordinatesIntrinsicTriangulation` is one concrete implementation (subclass) of [intrinsic triangulations](../basics). **Most of its functionality is documented there.**

### Constructors

??? func "`#!cpp IntegerCoordinatesIntrinsicTriangulation::IntegerCoordinatesIntrinsicTriangulation(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& inputGeom)`"

    Initialize an intrinsic triangulation sitting on top of `mesh`. Recall that `IntrinsicGeometryInterface` can be almost any geometry object, including a `VertexPositionGeometry`.

    Initially, the intrinsic triangulation will be identical to the input mesh; it can be modified with the routines below.
    
    Performs [mollification](/surface/algorithms/robust_geometry/#intrinsic-mollification) on the intrinsic trianglulation, using a relative factor of `2e-5`.


## Citation

```bib
@article{gillespie2021integer,
  title={Integer Coordinates for Intrinsic Geometry Processing},
  author={Gillespie, Mark and Sharp, Nicholas and Crane, Keenan},
  journal={arXiv preprint arXiv:2106.00220},
  year={2021}
}
```
