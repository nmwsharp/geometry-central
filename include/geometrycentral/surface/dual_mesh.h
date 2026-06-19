#pragma once


#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"


namespace geometrycentral
{
namespace surface
{

/**
 * @brief       Builds the dual manifold mesh of the input manifold mesh.
 * 
 * @details     This function computes the dual manifold mesh of the given primal 
 *              manifold mesh.\n 
 *              Each face of the primal mesh corresponds to a vertex in the dual mesh.
 *              If two faces are adjacent (i.e., share an edge) in the primal mesh,
 *              the corresponding dual vertices are connected by an edge in the dual mesh.\n 
 *              As a consequence, each vertex in the primal mesh corresponds to a face
 *              in the dual mesh (i.e., the dual vertices corresponding to fan of primal
 *              faces around the primal vertex).\n 
 *              If the primal manifold mesh has boundaries, it is impossible to preserve
 *              the vertex-face bijection in both directions.\n 
 *              If the parameter keepBoundaries is set to false, the bijection between
 *              primal faces and dual vertices is preserved, but the primal vertices at the
 *              boundary will not have corresponding faces, and the 3D embeddings of the
 *              primal and dual mesh would have different boundaries.\n 
 *              If the parameter keepBoundaries is set to true, additional dual vertices will
 *              be created at the primal boundary vertices and at the midpoints of the
 *              primal boundary edges. This preserves the bijection between primal vertices
 *              and dual faces, but the dual vertices at boundary will not have corresponding
 *              primal faces.
 * 
 * @param mesh The mesh whose the dual must be computed.
 * @param keepBoundaries Whether or not to create additional dual vertices to preserve the boundary embedding.
 * @return std::unique_ptr<geometrycentral::surface::ManifoldSurfaceMesh> The resulting dual mesh.
 */
std::unique_ptr<geometrycentral::surface::ManifoldSurfaceMesh>
dual_mesh(geometrycentral::surface::ManifoldSurfaceMesh& mesh,
          bool keepBoundaries);
    
/**
 * @brief       Builds the dual manifold mesh of the input manifold mesh, with the corresponding
 *              vertex geometry.
 * 
 * @details     This function computes the dual manifold mesh of the given primal 
 *              manifold mesh, with the corresponding vertex geometry.\n 
 *              Each face of the primal mesh corresponds to a vertex in the dual mesh, whose
 *              position in 3D space is the barycenter of the primal face.
 *              If two faces are adjacent (i.e., share an edge) in the primal mesh,
 *              the corresponding dual vertices are connected by an edge in the dual mesh.\n 
 *              As a consequence, each vertex in the primal mesh corresponds to a face
 *              in the dual mesh (i.e., the dual vertices corresponding to fan of primal
 *              faces around the primal vertex).\n 
 *              If the primal manifold mesh has boundaries, it is impossible to preserve
 *              the vertex-face bijection in both directions.\n 
 *              If the parameter keepBoundaries is set to false, the bijection between
 *              primal faces and dual vertices is preserved, but the primal vertices at the
 *              boundary will not have corresponding faces, and the 3D embeddings of the
 *              primal and dual mesh would have different boundaries.\n 
 *              If the parameter keepBoundaries is set to true, additional dual vertices will
 *              be created at the primal boundary vertices and at the midpoints of the
 *              primal boundary edges. This preserves the bijection between primal vertices
 *              and dual faces, but the dual vertices at boundary will not have corresponding
 *              primal faces.
 * 
 * @param mesh The mesh whose the dual must be computed.
 * @param geometry The geometry of the primal mesh.
 * @param keepBoundaries Whether or not to create additional dual vertices to preserve the boundary embedding.
 * @return std::tuple<std::unique_ptr<geometrycentral::surface::ManifoldSurfaceMesh>,
 * std::unique_ptr<geometrycentral::surface::VertexPositionGeometry>> The tuple containing the dual mesh and its embedded geometry.
 */
std::tuple<std::unique_ptr<geometrycentral::surface::ManifoldSurfaceMesh>, 
           std::unique_ptr<geometrycentral::surface::VertexPositionGeometry>>
dual_mesh(geometrycentral::surface::ManifoldSurfaceMesh& mesh,
          geometrycentral::surface::VertexPositionGeometry& geometry,
          bool keepBoundaries);

} // namespace surface
} // namespace geometrycentral