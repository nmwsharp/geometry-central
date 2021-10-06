This tutorial walks through some basics with geometry-central, showing how to load a mesh from file and iterate through its elements.

[View full, runnable source code in the tutorial repository.](https://github.com/nmwsharp/geometry-central-tutorials)

### Basic setup

To begin, we include the relevant headers, including some for visualization using [Polyscope](https://polyscope.run).

```cpp
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/surface_mesh.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
```

All functionality of geometry-central is contained within the `geometrycentral` namespace; surface meshes live in `geometrycentral::surface`. We will bring both in to scope so we can just type `SurfaceMesh` instead of `geometrycentral::surface::SurfaceMesh`, etc.

```cpp
using namespace geometrycentral;
using namespace geometrycentral::surface;
```

Now, we use the [mesh loaders](/surface/utilities/io) to construct a surface mesh's connectivity and geometry from file. Many common file formats like `.obj`, `.ply`, and `.stl` are supported. By default, this mesh class represents very general polygonal meshes, including nonmanifold meshes.

```cpp
  std::unique_ptr<SurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geometry;
  std::tie(mesh, geometry) = readSurfaceMesh("spot.obj");
``` 

Alternately, you could construct a mesh which is required to be manifold.
```cpp
  std::tie(mesh, geometry) = readManifoldSurfaceMesh("spot.obj");
``` 

If you'r not already familiar with `std::unique_ptr<>`, the following gives a bit more context (click to expand).

??? note "Why use `std::unique_ptr<>`?"

    The mesh loader, like many functions in geometry-central, returns constructed objects via a `unique_ptr`. Unique pointers are an important tool for memory management in modern C++; if you haven't used them before, we suggest you give them a try!

    In most ways, a `unique_ptr` acts just like a normal C++ pointer. You can dereference it with `*uPtr`, and access its members and function like `uPtr->function()`. However, the `unique_ptr` helps prevent common memory-management mistakes, and communicates the programmer's intent about object lifetime. This is accomplished with two properties:

      - You don't need to call `delete` on a `unique_ptr`, it happens automatically when the pointer is destructed, e.g. when it goes out of scope at the end of a function, or when the object it is a member of gets deleted. This helps prevent memory leaks where you forget to deallocate the object.

      - You cannot copy the `unique_ptr`; hence it is "unique"! You can still pass around references, or `std::move()` the pointer, which are sufficient for most reasonable uses. This helps prevent you from creating a copy, and then accidentally deleting the pointer twice.

    ---

    The general paradigm in geometry-central (and a recommended style in modern C++) is to return long-lived, allocated objects with a `unique_ptr`, and pass these objects in to functions and dependent classes by reference.

    For instance, we might write a function which takes a mesh as an argument like

    ```cpp
    void processMesh(SurfaceMesh& inputMesh) { /* do stuff */}
    ```
   
    and call it by dereferencing the unique pointer to pass a reference

    ```cpp
    std::unique_ptr<SurfaceMesh> mesh;
    std::unique_ptr<VertexPositionGeometry> geometry;
    std::tie(mesh, geometry) = readSurfaceMesh("spot.obj"); 
    
    processMesh(*mesh);
    ```

    For more details about unique pointers, see the [language documentation](https://en.cppreference.com/w/cpp/memory/unique_ptr), or many tutorials around the web.

    --- 

    If you really don't want to use unique pointers, you can simply release the unique pointer to an ordinary pointer:

    ```cpp
    std::unique_ptr<SurfaceMesh> mesh /* populated as above */;
    SurfaceMesh* meshPtr = mesh.release();
    ```
    
    The `meshPtr` now points the mesh object, and you are responsible for eventually deleting this pointer. After calling `release()`, the unique pointer points to nothing and will no longer deallocate the object.


### Traverse the mesh
  
As a simple demonstration of the mesh data structure, lets iterate through the vertices of the mesh, and for each vertex print the adjacent faces. 
```cpp 
  for (Vertex v : mesh->vertices()) {
    std::cout << "Vertex " << v << " has degree " << v.degree() << "\n";
    for (Face fn : v.adjacentFaces()) {
      std::cout << "  incident on face " << fn << "\n";
    }
  }
```

This prints something like:

```
...
Vertex v_2907 has degree 6
  incident on face f_5815
  incident on face f_2885
  incident on face f_2887
  incident on face f_5812
  incident on face f_5813
  incident on face f_5814
Vertex v_2908 has degree 6
  incident on face f_5822
  incident on face f_2888
...
```



### Visualize the mesh with Polyscope

Finally, we can easily visualize the mesh we loaded via Polyscope.

```cpp
  polyscope::init(); // initialize the gui

  // add the mesh to the gui
  polyscope::registerSurfaceMesh("my mesh", 
      geometry->vertexPositions, mesh->getFaceVertexList());

  polyscope::show(); // pass control to the gui until the user exits
```

