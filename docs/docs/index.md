# Welcome to Geometry Central

Geometry-central is a modern C++ library of data structures and algorithms for geometry processing, with a particular focus on surface meshes.

[![actions status linux](https://github.com/nmwsharp/geometry-central/workflows/linux/badge.svg)](https://github.com/nmwsharp/geometry-central/actions)
[![actions status macOS](https://github.com/nmwsharp/geometry-central/workflows/macOS/badge.svg)](https://github.com/nmwsharp/geometry-central/actions)
[![actions status windows](https://github.com/nmwsharp/geometry-central/workflows/windows/badge.svg)](https://github.com/nmwsharp/geometry-central/actions)

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

For more, see the [tutorials](/tutorials/load_mesh). To get started with the code, see [building](/build/building). Use the [sample project](https://github.com/nmwsharp/gc-polyscope-project-template/) to get started with a build system and a gui.

A introductory talk on geometry-central was given at SGP 2020, check it out to get started: [www.youtube.com/watch?v=mw5Xz9CFZ7A](https://www.youtube.com/watch?v=mw5Xz9CFZ7A)


**Bindings & Plugins:**

- **Python:** [Potpourri3d](https://github.com/nmwsharp/potpourri3d)
- **Grasshopper/Rhino:** [Lionfish](https://www.food4rhino.com/app/lion-fish) by Math Whittaker

If you're interested in creating additional bindings/plugins, feel free to reach out!


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
[Rohan Sawhney](http://rohansawhney.io/), 
[Chris Yu](https://www.cs.cmu.edu/~christoy/),
and many others.

If geometry-central contributes to an academic publication, cite it as:
```bib
@article{geometrycentral,
  title={GeometryCentral: A modern C++ library of data structures and algorithms for geometry processing},
  author={Nicholas Sharp and Keenan Crane and others},
  howpublished="\url{https://geometry-central.net/}",
  year={2019}
}
```

Development of this software was funded in part by NSF Award 1717320, an NSF graduate research fellowship, and gifts from Adobe Research and Autodesk, Inc.
