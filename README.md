# Welcome to Geometry Central

Geometry Central is a modern C++ codebase, providing the low-level tools to implement algorithms in geometry processing, scientific computing, and computer graphics/vision, with a particular focus on the geometry of surfaces.

Features include:
- A polished *halfedge mesh* class, and a system of containers for associating data with mesh elements.
- Implementations of canonical *geometric quantities* on surfaces, such as normals and curvatures.
- Tools for generating and manipulating *surface parameterizations*.
- The basic building blocks of *discrete exterior calculus*.
- A coherent set of sparse *linear algebra tools*, currently based on Eigen, and augmented with solvers that automatically utilize SuiteSparse if available on your system.


## What is it not?
- **A user interface**. Geometry Central does not include any facilities for user interaction; it is an algorithms and data structures library on which you might build user-facing tools. This philosphy keeps the library lightweight, and avoids dependencies on rendering and windowing systems. For a UI built on top of Geometry Central, see [Polyscope](https://github.com/nmwsharp/polyscope).
- **A research code dump**. Geometry Central was built by researchers, and is used in many prototype research projects. However, we strive to ensure that this library contains only polished, broadly useful algorithms, rather than a mashup of one-off research ideas. See the list below for projects built on Geometry Central.

## Examples

TODO

## Authors:
- [Nick Sharp](http://nmwsharp.com)
- [Keenan Crane](http://keenan.is/here)
- [Rohan Sawheny](http://rohansawhney.io/)
- And the rest of the [Geometry Collective](http://geometry.cs.cmu.edu) at Carnegie Mellon University

Development of this software was funded in part by NSF Award 1717320, an NSF graduate research fellowship, and gifts from Adobe Research and Autodesk, Inc.

## Projects that use Geometry Central

TODO
