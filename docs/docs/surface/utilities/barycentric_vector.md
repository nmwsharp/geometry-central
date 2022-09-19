A `BarycentricVector` is a vector that lies along a surface. This vector can lie within a face, along an edge, or at a single vertex (in which case it is simply a zero vector.) Barycentric vectors represent the difference between two barycentric points on the surface, hence their coordinates always sum to 0. (Barycentric points are implemented as [SurfacePoint](../surface_point) in geometry-central.) Currently, barycentric points and vectors are only supported on triangle meshes.

![example of a barycentric vector within a face of an intrinsic triangulation](/media/barycentric_vector.svg)

Although a barycentric vector can be constructed as the difference of two barycentric points, a barycentric vector technically does not define a single unique vector along the surface but rather a constant vector field within a face.

Using barycentric vectors, one can easily do vector arithmetic on a surface. Barycentric vectors are especially useful when working with an intrinsic representation of a surface; they can be used to do computations on vectors that depend only on intrinsic geometry, such as inner products.

`#include "geometrycentral/surface/barycentric_vector.h"`

The field `BarycentricVector::type` is an enum:
```cpp
enum class BarycentricVectorType { Face = 0, Edge, Vertex };
```

which indicates what kind of vector it is.

- If the barycentric vector is **inside a face**, the field `BarycentricVector::face` indicates which face. Otherwise it is the null default face. The field `BarycentricVector::faceCoords` stores the coordinates of the vector, expressed in barycentric coordinates of the face. If one thinks of a barycentric vector as being the difference of two barycentric points, the order of the vector coordinates can be thought of as corresponding to the usual iteration order of vertices about the face. More technically, the vector coordinates represent vectors in the tangent plane of the _standard triangle_ (see Section 2.3 of [these notes](https://markjgillespie.com/Research/int-tri-course/int_tri_course.pdf).)

- If the barycentric vector is **along an edge**, the field `BarycentricVector::edge` indicates which edge. Otherwise it is the null default edge. The field `BarycentricVector::edgeCoords` stores the coordinates of the vector, expressed in barycentric coordinates of the edge (ordered according to the usual order of vertices on the edge, i.e. the first component corresponds to `edge.halfedge().vertex() == edge.firstVertex()`.)

- If the barycentric vector is on a **vertex**, the field `BarycentricVector::vertex` indicates which vertex. Otherwise it is the null default vertex. There are no corresponding coordinates, since a vector on a vertex is always considered to have zero length.

### Construction

??? func "`#!cpp BarycentricVector::BarycentricVector(Face f, Vector3 faceCoords)`"

    Construct a barycentric vector that lies along the given face `f`, with the given barycentric coordinates `faceCoords`.

??? func "`#!cpp BarycentricVector::BarycentricVector(Edge e, Vector2 edgeCoords)`"

    Construct a barycentric vector that lies along the given edge `e`, with the given barycentric coordinates `edgeCoords`.

??? func "`#!cpp BarycentricVector::BarycentricVector(Vertex v)`"

    Construct a (zero-length) barycentric vector that lies on the given vertex `v`.

??? func "`#!cpp BarycentricVector::BarycentricVector(SurfacePoint pA, SurfacePoint pB)`"

    Construct a barycentric vector from the given [surface points](../surface_point), i.e. barycentric points. The vector direction goes from `pA` to `pB`.

    Note that there do not exist unique barycentric endpoints for a given barycentric vector, except for the degenerate case of a vertex-type vector where both endpoints are the vertex itself. Perhaps a more straightforward way to say this is that barycentric vectors really are **vectors** (displacements), **not** rays.

### Arithmetic
Barycentric vectors support addition, subtraction, scalar multiplication, and scalar division.

### Member operations

??? func "`#!cpp double BarycentricVector::norm(IntrinsicGeometryInterface& geom) const`"

    Returns the norm (length) of the vector. Because the norm depends on the geometry, you must pass in the geometry (either intrinsic or extrinsic) on which the vector is defined.

??? func "`#!cpp double BarycentricVector::norm2(IntrinsicGeometryInterface& geom) const`"

    Returns the squared norm of the vector. Because the norm depends on the geometry, you must pass in the geometry (either intrinsic or extrinsic) on which the vector is defined.

??? func "`#!cpp BarycentricVector BarycentricVector::rotated90(IntrinsicGeometryInterface& geom) const`"

    Rotate the vector 90 degrees counter-clockwise within the face it belongs to. This requires the geometry, since the meaning of "90 degrees" depends on the geometry of the triangle.

Barycentric vectors also have a few utility functions:

??? func "`#!cpp BarycentricVector BarycentricVector::BarycentricVector inSomeFace() const`"

    All barycentric vectors (vertex, edge, face) have an equivalent vector in one or many adjacent faces. This function reeturns one of the equivalent barycentric vectors in a face (chosen arbitrarily). If the vector is already a face vector, the output is a copy of the vector.

??? func "`#!cpp BarycentricVector BarycentricVector::BarycentricVector inFace(Face f) const`"

    Returns the barycentric vector as a face vector in face `f`. If the vector is not on or adjacent to the requested face, throws an error.

??? func "`#!cpp BarycentricVector BarycentricVector::BarycentricVector inEdge(Edge e) const`"

    Returns the barycentric vector as an edge-type vector in edge `e`. If the vector is already an edge vector, the output is a copy of the vector. If the the barycentric vector is not on `e` or one of its adjacent vertices, throws an error.

### Function operations

??? func "`#!cpp double norm(IntrinsicGeometryInterface& geom, const BarycentricVector& v)`"

    Returns the norm of the input vector `v`.

    Also available as `vec.norm()`.

??? func "`#!cpp double norm2(IntrinsicGeometryInterface& geom, const BarycentricVector& v)`"

    Returns the squared norm of the input vector `v`.

    Also available as `vec.norm2()`.

??? func "`#!cpp double dot(IntrinsicGeometryInterface& geom, const BarycentricVector& u, const BarycentricVector& v)`"

    Returns the inner product between `u` and `v`.

??? func "`#!cpp Face sharedFace(const BarycentricVector& u, const BarycentricVector& v)`"

    Returns an arbitrary face shared by the two vectors `u` and `v`, if one exists; returns the null default face if none.

??? func "`#!cpp Edge sharedEdge(const BarycentricVector& u, const BarycentricVector& v)`"

    Returns the edge shared by the two vectors `u` and `v`, if one exists; return the null default edge if none.