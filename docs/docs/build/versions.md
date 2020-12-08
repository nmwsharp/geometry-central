

### Changes in main branch since last version

- Change the semantics of `MeshData<>::size()` to match what `size()` usually means

### July 3, 2020

This version:

- Generalize the main halfedge mesh type to support nonmanifold meshes in routines where they make sense. The old `HalfedgeMesh` is now `ManifoldSurfaceMesh`, which is a subclass of the new more general `SurfaceMesh`, offering many of the same operations. The header `halfedge_mesh.h` typedef's `HalfedgeMesh` as `ManifoldSurfaceMesh` so existing code will mostly still work.
- Renamed `PolygonSoupMesh` to `SimplePolygonMesh`, and simplified some methods of this class. For now, the old type `PolygonSoupMesh` is typedef'd to `SimplePolygonMesh`, and the header `polygon_soup_mesh.h` includes `simple_polygon_mesh.h` so existing code should work. Please use `SimplePolygonMesh` in any new code.
- Renamed `PlyHalfedgeMeshData` to `RichSurfaceMeshData`, and changed its workings to apply to more general meshes.
- Changed underlying storage of `MeshData<>` containers from `std::vector<>` to `Eigen::VectorX_`.
- Moved `halfedge_containers.h` to `utilities/mesh_data.h`, along with reorganizing various mesh element headers (you shouldn't need to include any of these headers in user code anyway, just including `surface_mesh.h` is sufficient)
