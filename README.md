# Documention is hosted at [geometry-central.net](http://geometry-central.net)
---

# Welcome to Geometry Central

[![actions status linux](https://github.com/nmwsharp/geometry-central/workflows/linux/badge.svg)](https://github.com/nmwsharp/geometry-central/actions)
[![actions status macOS](https://github.com/nmwsharp/geometry-central/workflows/macOS/badge.svg)](https://github.com/nmwsharp/geometry-central/actions)
[![actions status windows](https://github.com/nmwsharp/geometry-central/workflows/windows/badge.svg)](https://github.com/nmwsharp/geometry-central/actions)

Geometry-central is a modern C++ library of data structures and algorithms for geometry processing, with a particular focus on surface meshes.

Features include:

- A polished **surface mesh** class, with efficient support for mesh modification, and a system of containers for associating data with mesh elements.
- Implementations of canonical **geometric quantities** on surfaces, ranging from normals and curvatures to tangent vector bases to operators from discrete differential geometry.
- A suite of **powerful algorithms**, including computing distances on surface, generating direction fields, and manipulating intrinsic Delaunay triangulations.
- A coherent set of sparse **linear algebra tools**, based on Eigen and augmented to automatically utilize better solvers if available on your system.


**Sample:**

```cpp
// Load a mesh
std::unique_ptr<SurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;
std::tie(mesh, geometry) = readSurfaceMesh("spot.obj"); 

// Compute vertex areas
VertexData<double> vertexAreas(*mesh);

geometry->requireFaceAreas();
for(Vertex v : mesh->vertices()) {
  double A = 0.;
  for(Face f : v.adjacentFaces()) {
    A += geometry->faceAreas[f] / v.degree();
  }
  vertexAreas[v] = A;
}
```

Check out the docs, tutorials, and build instructions at [geometry-central.net](http://geometry-central.net).  Use the [sample project](https://github.com/nmwsharp/gc-polyscope-project-template/) to get started with a build system and a gui.


**Related alternatives:** 
[CGAL](https://www.cgal.org/),
[libIGL](https://github.com/libigl/libigl),
[OpenMesh](http://www.openmesh.org/),
[Polygon Mesh Processing Library](https://www.pmp-library.org/),
[CinoLib](https://github.com/mlivesu/cinolib)

---

**Credits**

Geometry-central is developed by [Nicholas Sharp](http://nmwsharp.com), with many contributions from 
[Keenan Crane](http://keenan.is/here), 
[Yousuf Soliman](http://www.its.caltech.edu/~ysoliman/),
[Mark Gillespie](http://markjgillespie.com/),
[Rohan Sawhney](http://rohansawhney.io/), and many others.



If geometry-central contributes to an academic publication, cite it as:
```bib
@misc{geometrycentral,
  title = {geometry-central},
  author = {Nicholas Sharp and Keenan Crane and others},
  note = {www.geometry-central.net},
  year = {2019}
}
```

Development of this software was funded in part by NSF Award 1717320, an NSF graduate research fellowship, and gifts from Adobe Research and Autodesk, Inc.
