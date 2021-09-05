# Signpost Intrinsic Triangulation

![signpost example](/media/signposts.svg){: style="max-width: 25em; display: block; margin-left: auto; margin-right: auto;"}

The signpost data structure encodes an intrinsic triangulation by storing "signposts" at mesh vertices. Explicitly, for each intrinsic edge it stores the edge's length and direction at the two incident vertices. This information fully specifies how the intrinsic triangulation sits above the input mesh. For more details, see [Navigating Intrinsic Triangulations](http://www.cs.cmu.edu/~kmcrane/Projects/NavigatingIntrinsicTriangulations/paper.pdf).

The `SignpostIntrinsicTriangulation` is one concrete implementation (subclass) of [intrinsic triangulations](../basics). **Most of its functionality is documented there.**


### Constructors

??? func "`#!cpp SignpostIntrinsicTriangulation::SignpostIntrinsicTriangulation(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& inputGeom)`"

    Initialize an intrinsic triangulation sitting on top of `mesh`. Recall that `IntrinsicGeometryInterface` can be almost any geometry object, including a `VertexPositionGeometry`.

    Initially, the intrinsic triangulation will be identical to the input mesh; it can be modified with the routines below.

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
