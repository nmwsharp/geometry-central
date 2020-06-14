


???+ note "Version 2.0"

    - Renamed `PolygonSoupMesh` to `SimplePolygonMesh`, and simplified some methods of this class. For now, the old type `PolygonSoupMesh` is typedef'd to `SimplePolygonMesh`, and the header `polygon_soup_mesh.h` includes `simple_polygon_mesh.h` so existing code should work. Please use `SimplePolygonMesh` in any new code.
    - Changed underlying storage of `MeshData<>` containers from `std::vector<>` to `Eigen::VectorX_`.
    - Renamed `halfedge_containers.h` to `utilities/mesh_data.h`, along with reorganizing various mesh element headers (you shouldn't need to include any of these headers in user code anyway, just including `halfedge_mesh.h` is sufficient)
