#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/mesh_graph_algorithms.h"
#include "geometrycentral/utilities/disjoint_sets.h"

namespace geometrycentral {
namespace surface {


// Returns a brand new mesh cut along the edges, and a map from new halfedges --> old halfedges which can be used to
// evaluate correspondence. All halfedges on the resulting mesh will have come from some halfedge on the original mesh,
// with the exception of new exterior halfedges along the cut, which will be set to Halfedge().
std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, HalfedgeData<Halfedge>> cutAlongEdges(ManifoldSurfaceMesh& mesh,
                                                                                       const EdgeData<char>& cut);

} // namespace surface
} // namespace geometrycentral
