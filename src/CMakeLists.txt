cmake_minimum_required(VERSION 3.8.0)

# Maybe stop from CMAKEing in the wrong place
if (CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR)
    message(FATAL_ERROR "Source and build directories cannot be the same. Go use the /build directory.")
endif()

# List of .cpp files to be compiled
SET(SRCS

  surface/surface_mesh.cpp
  surface/manifold_surface_mesh.cpp
  surface/halfedge_factories.cpp  
  surface/surface_mesh_factories.cpp
  surface/meshio.cpp
  surface/simple_polygon_mesh.cpp
  surface/rich_surface_mesh_data.cpp

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
  surface/simple_idt.cpp
  surface/mesh_graph_algorithms.cpp
  surface/surface_point.cpp
  surface/fast_marching_method.cpp
  surface/uniformize.cpp
  surface/parameterize.cpp
  surface/surgery.cpp
  surface/simple_idt.cpp
  surface/exact_polyhedral_geodesics.cpp
  surface/intrinsic_mollification.cpp
  surface/tufted_laplacian.cpp
  surface/flip_geodesics.cpp
  #surface/detect_symmetry.cpp
  #surface/mesh_ray_tracer.cpp

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
  utilities/knn.cpp
)

SET(INCLUDE_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../include/geometrycentral/")
SET(HEADERS
  ${INCLUDE_ROOT}/numerical/linear_algebra_utilities.h
  ${INCLUDE_ROOT}/numerical/linear_algebra_utilities.ipp
  ${INCLUDE_ROOT}/numerical/linear_solvers.h
  ${INCLUDE_ROOT}/numerical/suitesparse_utilities.h

  ${INCLUDE_ROOT}/surface/barycentric_coordinate_helpers.h
  ${INCLUDE_ROOT}/surface/barycentric_coordinate_helpers.ipp
  ${INCLUDE_ROOT}/surface/base_geometry_interface.h
  ${INCLUDE_ROOT}/surface/detect_symmetry.h
  ${INCLUDE_ROOT}/surface/direction_fields.h
  ${INCLUDE_ROOT}/surface/edge_length_geometry.h
  ${INCLUDE_ROOT}/surface/edge_length_geometry.ipp
  ${INCLUDE_ROOT}/surface/embedded_geometry_interface.h
  ${INCLUDE_ROOT}/surface/exact_polyhedral_geodesics.h
  ${INCLUDE_ROOT}/surface/extrinsic_geometry_interface.h
  ${INCLUDE_ROOT}/surface/fast_marching_method.h
  ${INCLUDE_ROOT}/surface/halfedge_element_types.h
  ${INCLUDE_ROOT}/surface/halfedge_element_types.ipp
  ${INCLUDE_ROOT}/surface/halfedge_factories.h
  ${INCLUDE_ROOT}/surface/halfedge_logic_templates.ipp
  ${INCLUDE_ROOT}/surface/halfedge_mesh.h
  ${INCLUDE_ROOT}/surface/heat_method_distance.h
  ${INCLUDE_ROOT}/surface/intrinsic_geometry_interface.h
  ${INCLUDE_ROOT}/surface/intrinsic_mollification.h
  ${INCLUDE_ROOT}/surface/manifold_surface_mesh.h
  ${INCLUDE_ROOT}/surface/meshio.h
  ${INCLUDE_ROOT}/surface/mesh_graph_algorithms.h
  ${INCLUDE_ROOT}/surface/mesh_ray_tracer.h
  ${INCLUDE_ROOT}/surface/parameterize.h
  ${INCLUDE_ROOT}/surface/rich_surface_mesh_data.h
  ${INCLUDE_ROOT}/surface/rich_surface_mesh_data.ipp
  ${INCLUDE_ROOT}/surface/polygon_soup_mesh.h
  ${INCLUDE_ROOT}/surface/signpost_intrinsic_triangulation.h
  ${INCLUDE_ROOT}/surface/signpost_intrinsic_triangulation.ipp
  ${INCLUDE_ROOT}/surface/simple_idt.h
  ${INCLUDE_ROOT}/surface/simple_polygon_mesh.h
  ${INCLUDE_ROOT}/surface/surface_centers.h
  ${INCLUDE_ROOT}/surface/surface_mesh.h
  ${INCLUDE_ROOT}/surface/surface_mesh.ipp
  ${INCLUDE_ROOT}/surface/surface_point.h
  ${INCLUDE_ROOT}/surface/surface_point.ipp
  ${INCLUDE_ROOT}/surface/surgery.h
  ${INCLUDE_ROOT}/surface/trace_geodesic.h
  ${INCLUDE_ROOT}/surface/uniformize.h
  ${INCLUDE_ROOT}/surface/vector_heat_method.h
  ${INCLUDE_ROOT}/surface/vertex_position_geometry.h
  ${INCLUDE_ROOT}/surface/vertex_position_geometry.ipp

  ${INCLUDE_ROOT}/utilities/combining_hash_functions.h
  ${INCLUDE_ROOT}/utilities/curve.h
  ${INCLUDE_ROOT}/utilities/curve.ipp
  ${INCLUDE_ROOT}/utilities/dependent_quantity.h
  ${INCLUDE_ROOT}/utilities/dependent_quantity.ipp
  ${INCLUDE_ROOT}/utilities/disjoint_sets.h
  ${INCLUDE_ROOT}/utilities/eigen_interop_helpers.h
  ${INCLUDE_ROOT}/utilities/knn.h
  ${INCLUDE_ROOT}/utilities/mesh_data.h
  ${INCLUDE_ROOT}/utilities/mesh_data.ipp
  ${INCLUDE_ROOT}/utilities/quaternion.h
  ${INCLUDE_ROOT}/utilities/timing.h
  ${INCLUDE_ROOT}/utilities/utilities.h
  ${INCLUDE_ROOT}/utilities/vector2.h
  ${INCLUDE_ROOT}/utilities/vector2.ipp
  ${INCLUDE_ROOT}/utilities/vector3.h
  ${INCLUDE_ROOT}/utilities/vector3.ipp
)

# Create a single library for the project
add_library(geometry-central ${SRCS} ${HEADERS})

# Includes from this project
target_include_directories(geometry-central PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../include")

# Add all includes and link libraries from dependencies, which were populated in deps/CMakeLists.txt
target_link_libraries(geometry-central PUBLIC ${GC_DEP_LIBS})

# Set compiler properties for the library
set_property(TARGET geometry-central PROPERTY CXX_STANDARD 11)
set_property(TARGET geometry-central PROPERTY CXX_STANDARD_REQUIRED TRUE)
target_compile_definitions(geometry-central PUBLIC NOMINMAX _USE_MATH_DEFINES)

# Define CMAKE flag used in these sources (but should be kept OUT of headers)
if(GC_HAVE_SUITESPARSE)
  target_compile_definitions(geometry-central PUBLIC GC_HAVE_SUITESPARSE)
  target_include_directories(geometry-central PRIVATE ${SUITESPARSE_INCLUDE_DIRS})
  target_link_libraries(geometry-central PRIVATE ${SUITESPARSE_LIBRARIES})
endif()
