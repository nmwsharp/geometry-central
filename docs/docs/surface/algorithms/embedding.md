![A circle and a triangle glued together to form a convex polyhedron](/media/embed_convex.jpg)

Isometric embedding algorithms attempt to find 3D vertex positions for a mesh (i.e., an extrinsic embedding) that match a given set of edge lengths (i.e., an intrinsic metric).  Loosely speaking, these algorithms enable the 3D shape of an object to be recovered from its 2D UV coordinates—as long as the mesh was flattened without any stretching or distortion.  In general, an exact isometric embedding may not exist, or may require that the triangulation be modified.  For instance, the example above shows the shape obtained by gluing together a circular disk and a triangle of equal perimeter along their boundary.  In this case, we can successfully find an embedding, but must flip many edges in order for the disk to bend as needed.

At present, Geometry Central implements only one embedding algorithm, for _convex_ metrics.  In this case a solution always exists, and can always be found via an efficient algorithm.  A convex metric means that every vertex has non-negative [Gaussian curvature](/surface/geometry/quantities/#vertex-gaussian-curvature). [Alexandrov's uniqueness theorem](https://en.wikipedia.org/wiki/Alexandrov%27s_uniqueness_theorem) guarantees that any convex metric has an isometric embedding in three-dimensional space, which is unique up to rigid motions.  Geometry Central uses a method based on [Bobenko and Izmestiev's algorithm](https://arxiv.org/abs/math/0609447) to find such embeddings.

The input metric is purely _intrinsic_, and can hence be provided using any `IntrinsicGeometryInterface`.  It can also be provided by specifying edge lengths or UV coordinates, as described below.  In all cases, (i) the input mesh must be triangulated, (ii) edge lengths must be positive, (iii) edge lengths must satisfy the [triangle inequality](https://en.wikipedia.org/wiki/Triangle_inequality) in each face, and (iv) [vertex Gaussian curvature](/surface/geometry/quantities/#vertex-gaussian-curvature) must be positive.

The output is an ordinary triangle mesh with vertex positions in three-dimensional space, encoded by a `ManifoldSurfaceMesh` and `VertexPositionGeometry`.  Optionally (e.g., for debugging and visualization), one can also recover vertex positions at corners for intermediate stages of the embedding, and the correspondence between the input and embedded mesh, encoded by an [intrinsic triangulation](/intrinsic_triangulations/basics/).

!!! warning "uv coordinates"

    When a metric is read from UV coordinates (e.g., in an OBJ file), lengths should be equal on opposite sides of each seam edge.  These lengths will automatically be averaged if they are not exactly equal.  Also be careful that coordinates written to disk are not truncated to a small number of digits, which can make flat vertices appear to have slight negative curvature.

The basic idea of Bobenko and Izmestiev's algorithm is to create an intrinsic tetrahedral mesh, where every triangle of the given surface mesh is connected to a central "apex" by radial edges that are initially very long, creating long skinny tetrahedra.  The radii of these edges are then optimized until the dihedral angle sum around all of them is equal to 2π, allowing the tetrahedra to be laid out in three-dimensional space without overlap.  During optimization, edge flips are applied to ensure that edges of the boundary (i.e., the surface mesh) remain convex rather than concave.

Note that the Geometry Central implementation of this algorithm is based in part on the [Java reference implementation](https://gitlab.discretization.de/sechel/f1software/-/tree/master/src/alexandrov) by Stefan Sechelmann.  However, our implementation makes some significant modifications not advocated by the original authors.  Users interested in making comparisons (e.g., for academic research papers) are advised to consult the reference implementation.

!!! warning "convergence"

    Although the convex embedding algorithm is guaranteed to work in exact arithmetic, it can still fail to converge due to floating point error (e.g., when triangulations have near-degenerate elements, or are slightly non-convex due to truncation error in the input data).  In this case, it can sometimes help to loosen the tolerances and/or increase the number of iterations in the [embedder options](#options).  Additional diagnostic information is provided by setting `EmbedConvexOptions::verbose = true`.

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

    Embed a metric described by a set of edge lengths.  Options are passed as a [RemeshOptions](#options) object.  If embedding fails, the flag `EmbedConvexResult::success` will be set to `false`.  Input mesh `mesh` and edge lengths `edgeLengths` are unchanged.

??? func "`#!cpp EmbedConvexResult embedConvex(ManifoldSurfaceMesh& mesh, CornerData<Vector2>& textureCoordinates, EmbedConvexOptions options = defaultEmbedConvexOptions);`"

    Embed a metric described by a set of UV coordinates.  Options are passed as a [RemeshOptions](#options) object.  If embedding fails, the flag `EmbedConvexResult::success` will be set to `false`.  Note that since UVs are in general discontinuous, the two lengths across any seam edge must agree (and will be averaged if not).  Note also that shape recovery from UVs will work only if the mesh was unfolded isometrically, i.e., without any stretching or distortion of edges.  Input mesh `mesh` and UVs `textureCoordinates` are unchanged.
    
??? func "`#!cpp EmbedConvexResult embedConvex(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& intrinsicGeometry, EmbedConvexOptions options = defaultEmbedConvexOptions);`"

    Embed a metric described by an `IntrinsicGeometryInterface`.  Note that this could be a purely intrinsic description such as `EdgeLengthGeometry`, but also an extrinsic description such as `VertexPositionGeometry` (e.g., for debugging/sanity checking).  Options are passed as a [RemeshOptions](#options) object.  If embedding fails, the flag `EmbedConvexResult::success` will be set to `false`.  Input mesh `mesh` and metric `IntrinsicGeometry` are unchanged.

### Advanced Usage

![Animation of embedding progress](/media/convex_embed.gif)

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

    Create an embedder from a metric described by a set of edge lengths.  Options are passed as a [RemeshOptions](#options) object.  Input mesh `mesh` and edge lengths `edgeLengths` are unchanged.

??? func "`#!cpp ConvexEmbedder::ConvexEmbedder(ManifoldSurfaceMesh& mesh, CornerData<Vector2>& textureCoordinates, EmbedConvexOptions options = defaultEmbedConvexOptions);`"

    Create an embedder from a metric described by a set of UV coordinates.  Options are passed as a [RemeshOptions](#options) object.  Note that since UVs are in general discontinuous, the two lengths across any seam edge must agree (and will be averaged if not).  Note also that shape recovery from UVs will work only if the mesh was unfolded isometrically, i.e., without any stretching or distortion of edges.  Input mesh `mesh` and UVs `textureCoordinates` are unchanged.
    
??? func "`#!cpp ConvexEmbedder::ConvexEmbedder(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& intrinsicGeometry, EmbedConvexOptions options = defaultEmbedConvexOptions);`"

    Create an embedder from a metric described by an `IntrinsicGeometryInterface`.  Note that this could be a purely intrinsic description such as `EdgeLengthGeometry`, but also an extrinsic description such as `VertexPositionGeometry` (e.g., for debugging/sanity checking).  Options are passed as a [RemeshOptions](#options) object.  Input mesh `mesh` and metric `IntrinsicGeometry` are unchanged.

#### Embedding functions

??? func "`#!cpp bool ConvexEmbedder::init();`"

    Initialize the embedding, returning false if the given triangulation and metric do not meet the criteria for convex embedding (no boundary, triangulated, genus zero, and positive angle defect at all vertices).

??? func "`#!cpp bool ConvexEmbedder::embed();`"

    Embed the given metric, returning true if embedding was successful.  Note that this method "fails gracefully" in the sense that if the embedder can't find an embedding that meets the specified tolerance, it will still construct vertex positions that are as close as possible, by averaging the most recent corner coordinates.

??? func "`#!cpp void ConvexEmbedder::stepEmbedding();`"

    Take a single embedding step.

??? func "`#!cpp void ConvexEmbedder::refreshVertexCoordinates();`"

    Compute extrinsic vertex positions at corners from the current intrinsic embedding, stored in `localLayout`.  Also compute average coordinates at vertices, stored in `embeding`.

#### Embedder data

The current state of the embedding is encoded by edge lengths on a tetrahedral mesh connecting all surface triangles to a central "apex" vertex (`ConvexEmbedder::apex`).  The edge lengths for this tetrahedral mesh are split up into lengths for each edge of a surface mesh (`ConvexEmbedder::intrinsicTriangulation`), as well as radii associated with each vertex of the surface mesh (`ConvexEmbedder::r`).  At any moment, one can update corner and vertex coordinates by calling `ConvexEmbedder::refreshVertexCoordinates()`.  The input metric is also stored in order to, e.g., visualize correspondence between the input metric and embedded triangulation.

??? func "Input metric"

    * `ManifoldSurfaceMesh ConvexEmbedder::originalMesh` — the triangulation used to define the input metric
    * `IntrinsicGeometryInterface* ConvexEmbedder::originalMetric` - the input metric, which determines edge lengths for `originalMesh`

??? func "Current metric"
    
    * `IntegerCoordinatesIntrinsicTriangulation* intrinsicTriangulation` — the triangulation and edge lengths defining the current metric

??? func "Embedding"
    
    * `VertexData<double> r` — length of radial edges connecting each vertex to the central apex
    * `CornerData<Vector3> localLayout` — embedding of current metric (will be discontinuous for non-flat metrics). Update by calling `ConvexEmbedder::refreshVertexCoordinates()`.
    * `VertexData<Vector3> embedding` — embedding of current metric (average of discontinuous values).  Update by calling `ConvexEmbedder::refreshVertexCoordinates()`.
    * `Vector3 apex` — location of central vertex, which makes tetrahedra with each triangle.  Update by calling `ConvexEmbedder::refreshVertexCoordinates()`.

## Helper Types

### Options

Options are passed in to `embedConvex()` and constructors for `ConvexEmbedder` via a `EmbedConvexOptions` object.

| Field                                   | Default value   | Meaning                                                     |
|-----------------------------------------|-----------------|-------------------------------------------------------------|
| `#!cpp double metricConvexityTolerance` | `1e-5`          | how negative do we allow input Gaussian curvature to be?    |
| `#!cpp double edgeConvexityTolerance`   | `1e-3`          | how negative do we allow mean curvature to be?              |
| `#!cpp double initialStepSize`          | `1`             | initial step size for an embedding step                     |
| `#!cpp int maxSteps`                    | `20`            | maximum number of embedding steps                           |
| `#!cpp int maxLineSearchSteps`          | `100`           | maximum number of line search steps                         |
| `#!cpp double embeddingTolerance`       | `1e-3`          | maximum angle defect for any radial edge                    |
| `#!cpp int maxNewtonIterations`         | `20`            | maximum number of steps for inner Newton solver             |
| `#!cpp int maxNewtonLineSearchSteps`    | `32`            | maximum number of line search steps for inner Newton solver |
| `#!cpp double newtonTolerance`          | `1e-4`          | l2 norm of residual for convergence of Newton's method      |
| `#!cpp bool verbose`                    | `false`         | whether to display diagnostic output

