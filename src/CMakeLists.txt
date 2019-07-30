cmake_minimum_required(VERSION 3.8.0)

# Maybe stop from CMAKEing in the wrong place
if (CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR)
    message(FATAL_ERROR "Source and build directories cannot be the same. Go use the /build directory.")
endif()

# List of .cpp files to be compiled
SET(SRCS

  surface/halfedge_mesh.cpp
  surface/halfedge_factories.cpp
  surface/meshio.cpp
  surface/polygon_soup_mesh.cpp
  surface/ply_halfedge_mesh_data.cpp
  
  surface/base_geometry_interface.cpp
  surface/intrinsic_geometry_interface.cpp
  surface/extrinsic_geometry_interface.cpp
  surface/embedded_geometry_interface.cpp
  surface/edge_length_geometry.cpp
  surface/vertex_position_geometry.cpp
  surface/direction_fields.cpp
  surface/heat_method_distance.cpp
  surface/vector_heat_method.cpp
  surface/trace_geodesic.cpp
  surface/surface_centers.cpp
  surface/signpost_intrinsic_triangulation.cpp
  #surface/mesh_graph_algorithms.cpp
  #surface/detect_symmetry.cpp
  #surface/mesh_ray_tracer.cpp
  #surface/exact_polyhedral_geodesics.cpp
  #surface/fast_marching_method.cpp

  numerical/linear_algebra_utilities.cpp
  numerical/suitesparse_utilities.cpp
  numerical/linear_solvers.cpp
  numerical/eigenproblem_solvers.cpp
  numerical/qr_solvers.cpp
  numerical/square_solvers.cpp
  numerical/positive_definite_solvers.cpp

  utilities/utilities.cpp
  utilities/quaternion.cpp
  utilities/disjoint_sets.cpp
)


# Create a single library for the project
add_library(
        geometry-central 
        ${SRCS}
        )

target_include_directories(geometry-central PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../include")
target_include_directories(geometry-central PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../deps/eigen-git-mirror")
target_include_directories(geometry-central PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../deps/nanort/include")
target_include_directories(geometry-central PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../deps/nanoflann/include")
target_include_directories(geometry-central PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../deps/happly")
if(GC_HAVE_SUITESPARSE)
  target_include_directories(geometry-central PUBLIC "${SUITESPARSE_INCLUDE_DIRS}")
endif()

# If we (optionally) resolved suitesparse support, we want to link against that as well
if(GC_HAVE_SUITESPARSE)
    target_link_libraries(geometry-central ${SUITESPARSE_LIBRARIES})
endif(GC_HAVE_SUITESPARSE)


# Set compiler properties for the library
set_property(TARGET geometry-central PROPERTY CXX_STANDARD 11)
