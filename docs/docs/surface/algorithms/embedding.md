![A circle and a triangle glued together to form a convex polyhedron](/media/embed_convex.jpg)

Isometric embedding algorithms attempt to find 3D vertex positions for a mesh (i.e., an extrinsic embedding) that match a given set of edge lengths (i.e., an intrinsic metric).  Loosely speaking, these algorithms enable the 3D shape of an object to be recovered from its 2D UV coordinates—as long as the mesh was flattened without any stretching or distortion.  In general, an exact isometric embedding may not exist, or may require that the triangulation be modified.  For instance, the example above shows the shape obtained by gluing together a circular disk and a triangle of equal perimeter along their boundary.  In this case, we can successfully find an embedding, but must flip many edges in order for the disk to bend as needed.

At present, Geometry Central implements only one embedding algorithm, for _convex_ metrics.  In this case a solution always exists, and can always be found via an efficient algorithm.  A convex metric means that every vertex has non-negative [Gaussian curvature](/surface/geometry/quantities/#vertex-gaussian-curvature). [Alexandrov's uniqueness theorem](https://en.wikipedia.org/wiki/Alexandrov%27s_uniqueness_theorem) guarantees that any convex metric has an isometric embedding in three-dimensional space, which is unique up to rigid motions.  Geometry Central uses a method based on [Bobenko and Izmestiev's algorithm](https://arxiv.org/abs/math/0609447) to find such embeddings.

The input metric is purely _intrinsic_, and can hence be provided using any `IntrinsicGeometryInterface`.  It can also be provided by specifying edge lengths or UV coordinates, as described below.  In all cases, (i) the input mesh must be triangulated, (ii) edge lengths must be positive, (iii) edge lengths must satisfy the [triangle inequality](https://en.wikipedia.org/wiki/Triangle_inequality) in each face, and (iv) [vertex Gaussian curvature](/surface/geometry/quantities/#vertex-gaussian-curvature) must be positive.

The output is an ordinary triangle mesh with vertex positions in three-dimensional space, encoded by a `ManifoldSurfaceMesh` and `VertexPositionGeometry`.  Optionally (e.g., for debugging and visualization), one can also recover vertex positions at corners for intermediate stages of the embedding, and the correspondence between the input and embedded mesh, encoded by an [intrinsic triangulation](/intrinsic_triangulations/basics/).

!!! warning "uv coordinates"

    When a metric is read from UV coordinates (e.g., in an OBJ file), lengths should be equal on opposite sides of each seam edge.  These lengths will automatically be averaged if they are not exactly equal.  Also be careful that coordinates written to disk are not truncated to a small number of digits, which can make flat vertices appear to have slight negative curvature.

The basic idea of Bobenko and Izmestiev's algorithm is to create an intrinsic tetrahedral mesh, where every triangle of the given surface mesh is connected to a central "apex" by radial edges that are initially very long, creating long skinny tetrahedra.  The radii of these edges are then optimized until the dihedral angle sum around all of them is equal to 2π, allowing the tetrahedra to be laid out in three-dimensional space without overlap.  During optimization, edge flips are applied to ensure that edges of the boundary (i.e., the surface mesh) remain convex rather than concave.

Note that the Geometry Central implementation of this algorithm is based in part on the [Java reference implementation](https://gitlab.discretization.de/sechel/f1software/-/tree/master/src/alexandrov) by Stefan Sechelmann.  However, our implementation makes some significant modifications not advocated by the original authors.  Users interested in making comparisons (e.g., for academic research papers) are advised to consult the reference implementation.

### Basic Usage

The easiest way to compute an embedding is through "one shot" wrapper functions that take some description of the intrinsic metric as input, and produce an embedded mesh as output.

Example
```cpp
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/embed_convex.h"

// Load a mesh with UVs
// Note: vertex positions will be ignored (e.g., can all be set to zero),
// since our goal is to solve for vertex positions that exhibit the UV lengths.
std::unique_ptr<SurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;
std::unique_ptr<CornerData<Vector2>> uvs;
tie(mesh, geometry, uvs) = readParameterizedManifoldSurfaceMesh(filename);

// Compute the embedding
EmbedConvexResult result = embedConvex( *mesh, *uvs );

if( result.success ) {
    // do something with result.mesh and result.geometry
}
```

Additional options can be passed to the embedder by passing an optional `EmbedConvexOptions` struct to the embedding function:

```cpp
EmbedConvexOptions options;
options.tolerance = 1e-2; // use a looser tolerance than default
EmbedConvexResult result = embedConvex( *mesh, *uvs, options );
```

??? func "`#!cpp EmbedConvexResult embedConvex(ManifoldSurfaceMesh& mesh, EdgeData<double>& edgeLengths, EmbedConvexOptions options = defaultEmbedConvexOptions);`"

    Embed a metric described by a set of edge lengths.  Options are passed as a [RemeshOptions](#options) object.  If embedding fails, the flag `EmbedConvexResult::success` will be set to `false`.

??? func "`#!cpp EmbedConvexResult embedConvex(ManifoldSurfaceMesh& mesh, CornerData<Vector2>& textureCoordinates, EmbedConvexOptions options = defaultEmbedConvexOptions);`"

    Embed a metric described by a set of UV coordinates.  Options are passed as a [RemeshOptions](#options) object.  If embedding fails, the flag `EmbedConvexResult::success` will be set to `false`.  Note that since UVs are in general discontinuous, the two lengths across any seam edge must agree (and will be averaged if not).  Note also that shape recovery from UVs will work only if the mesh was unfolded isometrically, i.e., without any stretching or distortion of edges.
    
??? func "`#!cpp EmbedConvexResult embedConvex(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& intrinsicGeometry, EmbedConvexOptions options = defaultEmbedConvexOptions);`"

    Embed a metric described by an `IntrinsicGeometryInterface`.  Note that this could be a purely intrinsic description such as `EdgeLengthGeometry`, but also an extrinsic description such as `VertexPositionGeometry` (e.g., for debugging/sanity checking).  Options are passed as a [RemeshOptions](#options) object.  If embedding fails, the flag `EmbedConvexResult::success` will be set to `false`.

### Advanced Usage

![Animation of embedding progress](/media/embed_convex.gif)

More options are available when performing an embedding via the `ConvexEmbedder` class.  The basic pattern is:

```cpp
ConvexEmbedder embedder( *mesh, *uvs, options );
if( embedder.init() ) { // try initializing the embedding
    if( embedder.embed() ) { // try finding the embedding
        embedder.refreshVertexCoordinates(); // compute final vertex coordinates
        // now do something with members of embedder
    }
}
```

The separate `init()` and `embed()` method help provide diagnostics---if initialization fails, the mesh failed to meet the preconditions of the algorithm (e.g., non-triangular faces); if embedding fails, the embedder likely ran into numerical difficulty (and perhaps looser tolerances or more iterations are required).  Rather than embedding the whole surface all at once, you can also take one step at a time:

```cpp
for( int i = 0; i < n; i++ ) {
    embedder.stepEmbedding();
    embedder.refreshVertexCoordinates();
    // e.g., animate the partial results at each step
    // using embedder->localLayout
}
```

!!! warning "final vertex coordinates"

    Note that the embedding is primarily represented by the radii stored in `ConvexEmbedder::r`.  When using the class interface, neither vertex positions nor corner positions are computed automatically (since they are not needed by intermediate stages of the algorithm).  One must therefore call `ConvexEmbedder::refreshVertexCoordinates()` before attempting to access `ConvexEmbedder::embedding` or `ConvexEmbedder::localLayout`.

#### Constructors

??? func "`#!cpp ConvexEmbedder::ConvexEmbedder(ManifoldSurfaceMesh& mesh, EdgeData<double>& edgeLengths, EmbedConvexOptions options = defaultEmbedConvexOptions);`"

    Create an embedder from a metric described by a set of edge lengths.  Options are passed as a [RemeshOptions](#options) object.

??? func "`#!cpp ConvexEmbedder::ConvexEmbedder(ManifoldSurfaceMesh& mesh, CornerData<Vector2>& textureCoordinates, EmbedConvexOptions options = defaultEmbedConvexOptions);`"

    Create an embedder from a metric described by a set of UV coordinates.  Options are passed as a [RemeshOptions](#options) object.  Note that since UVs are in general discontinuous, the two lengths across any seam edge must agree (and will be averaged if not).  Note also that shape recovery from UVs will work only if the mesh was unfolded isometrically, i.e., without any stretching or distortion of edges.
    
??? func "`#!cpp ConvexEmbedder::ConvexEmbedder(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& intrinsicGeometry, EmbedConvexOptions options = defaultEmbedConvexOptions);`"

    Create an embedder from a metric described by an `IntrinsicGeometryInterface`.  Note that this could be a purely intrinsic description such as `EdgeLengthGeometry`, but also an extrinsic description such as `VertexPositionGeometry` (e.g., for debugging/sanity checking).  Options are passed as a [RemeshOptions](#options) object.

#### Embedding functions

??? func "`#!cpp bool ConvexEmbedder::init();`"

    Initialize the embedding, returning false if the given triangulation and metric do not meet the criteria for convex embedding (no boundary, triangulated, genus zero, and positive angle defect at all vertices).

??? func "`#!cpp bool ConvexEmbedder::embed();`"

    Embed the given metric, returning true if embedding was successful.  Note that this method "fails gracefully" in the sense that if the embedder can't find an embedding that meets the specified tolerance, it will still construct vertex positions that are as close as possible, by averaging the most recent corner coordinates.

??? func "`#!cpp void ConvexEmbedder::stepEmbedding();`"

    Take a single embedding step.

??? func "`#!cpp void ConvexEmbedder::refreshVertexCoordinates();`"

    Compute extrinsic vertex positions at corners from the current intrinsic embedding, stored in `localLayout`.  Also compute average coordinates at vertices, stored in `embeding`.

#### Embedding data

    TODO

## Helper Types
### Options
Options are passed in to `remesh` via a `EmbedConvexOptions` object.

TODO UPDATE THIS TABLE FOR CONVEX EMBEDDER

| Field                                              | Default value                         | Meaning                                                                                                                                                    |
|----------------------------------------------------|---------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `#!cpp double targetEdgeLength`                   | `-1`                                  | the target edge length in flat regions. If `targetEdgeLength` is negative, the target edge length is set to relative the input mesh's mean edge length     |
| `#!cpp size_t maxIterations`                      | `10`                                  | the maximum number of iterations to run for                                                                                                                |
| `#!cpp double curvatureAdaptation`                | `0`                                   | how much target length should vary due to curvature. Set curvatureAdaptation to 0 if you want edge lengths to be approximately targetEdgeLength everywhere |
| `#!cpp double minRelativeLength`                  | `0.05`                                | the minimum possible edge length allowed in the output mesh. Defined relative to targetEdgeLength                                                          |
| `#!cpp RemeshSmoothStyle smoothStyle`             | `RemeshSmoothStyle::Circumcentric`    | the type of vertex smoothing to use (either `RemeshSmoothStyle::Circumcentric` or `RemeshSmoothStyle::Laplacian`)                                          |
| `#!cpp RemeshBoundaryCondition boundaryCondition` | `RemeshBoundaryCondition::Tangential` | the type of motions allowed for boundary vertices (either `RemeshBoundaryCondition::Fixed`, `RemeshBoundaryCondition::Tangential` or `RemeshBoundaryCondition::Free`)                            |

