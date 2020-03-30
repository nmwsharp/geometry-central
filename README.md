# Documention is hosted at [geometry-central.net](http://geometry-central.net)
---

# Welcome to Geometry Central

Geometry-central is a modern C++ library of data structures and algorithms for geometry processing, with a particular focus on surface meshes.

Features include:

- A polished **halfedge mesh** class, with efficient support for mesh modification, and a system of containers for associating data with mesh elements.
- Implementations of canonical **geometric quantities** on surfaces, ranging from normals and curvatures to tangent vector bases to operators from discrete differential geometry.
- A suite of **powerful algorithms**, including computing distances on surface, generating direction fields, and manipulating intrinsic Delaunay triangulations.
- A coherent set of sparse **linear algebra tools**, based on Eigen and augmented to automatically utilize better solvers if available on your system.


**Sample:**

```cpp
// Load a mesh
std::unique_ptr<HalfedgeMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;
std::tie(mesh, geometry) = loadMesh("spot.obj");

// Compute vertex areas
VertexData<double> vertArea(*mesh, 0.);
geometry->requireFaceAreas();
for(Vertex v : mesh->vertices()) {
  for(Face f : v.adjacentVertices()) {
    vertArea[v] += geometry->faceAreas[f] / f.degree();
  }
}
```

For more, see the [tutorials](../tutorials/load_mesh). To get started with the code, see [building](../build/build).

**What is geometry-central not?**

- **A user interface**. Geometry-central does not include any facilities for user interaction; it is an algorithms and data structures library on which you might build user-facing tools. This philosophy keeps the library lightweight, and avoids dependencies on rendering and windowing systems. For a UI that interoperates well with geometry-central, see [Polyscope](https://polyscope.run).
- **A research code dump**. Geometry-central was built by researchers, and is used to prototype research projects. However, we strive to ensure that this library contains only polished & tested, broadly useful algorithms.


**Related alternatives:**
[CGAL](https://www.cgal.org/),
[libIGL](https://github.com/libigl/libigl),
[OpenMesh](http://www.openmesh.org/),
[Polygon Mesh Processing Library](https://www.pmp-library.org/),
[CinoLib](https://github.com/mlivesu/cinolib)

---

**Credits**

Geometry-central is developed primarily by [Nicholas Sharp](http://nmwsharp.com), with contributions from
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
