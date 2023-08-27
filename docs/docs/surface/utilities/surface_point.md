A `SurfacePoint` is a generic location on a surface, which might be at a vertex, along an edge, or inside a face. Surface points are used throughout geometry-central for methods that input or output arbitrary locations on surfaces.

`#include "geometrycentral/surface/surface_point.h"`

The field `SurfacePoint::type` is an enum:
```cpp
enum class SurfacePointType { Vertex = 0, Edge, Face };
```

which indicates what kind of point this is.

- if the surface point is a **vertex**, the field `SurfacePoint::vertex` indicates which vertex. Otherwise it is the null default vertex.

- if the surface point is **along an edge**, the field `SurfacePoint::edge` indicates which edge. Otherwise it is the null default edge. The field `SurfacePoint::tEdge` indicates the location along that edge, in the range `[0,1]`, with `0` at `edge.halfedge().vertex()`.

- if the surface point is **inside a face**, the field `SurfacePoint::face` indicates which face. Otherwise it is the null default face. The field `SurfacePoint::faceCoords` indicates the location inside that face, as barycentric coordinates (numbered according to the iteration order of vertices about the face, as usual).

Surface points have a few useful utility methods:

??? func "`#!cpp T SurfacePoint::interpolate(const VertexData<T>& data)`"

    Given data of tempalte type `T` defined at vertices, linearly interpolates to a value at this location.


??? func "`#!cpp SurfacePoint SurfacePoint::inSomeFace()`"

    All surface points (vertex, edge, face) have an equivalent point in one or many adjacent faces. For instance, a vertex could be equivalently a point in any of the incident faces, with a single `1` barycentric coordinate, or a point on an edge could be a point in either of the two adjacent faces.

    This function returns one of the equivalent surface points in a face (chosen arbitrarily). If this point is a face point, the output is a copy of this point.

??? func "`#!cpp Vertex SurfacePoint::nearestVertex()`"

    Returns the nearest vertex which is adjacent to this point.

    For surface points which are vertices, it will return the same vertex.  For surface points which are along edges, it will return one of the two incident vertices.  For surface points which are inside faces, it will return one of the three incident vertices.

In addition, surface points can be used to construct [barycentric vectors](../barycentric_vector).