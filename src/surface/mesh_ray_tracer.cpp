#include "geometrycentral/surface/mesh_ray_tracer.h"

#include <vector>

using std::cout;
using std::endl;

namespace geometrycentral {
namespace surface {

MeshRayTracer::MeshRayTracer(Geometry<Euclidean>* geometry_) {
  mesh = geometry_->getMesh();
  geometry = geometry_;

  buildBVH();
}

void MeshRayTracer::buildBVH() {
  cout << "Building BVH for mesh..." << endl;

  nanort::BVHBuildOptions<double> options; // Use default options

  if (!mesh->isTriangular()) {
    throw std::runtime_error("Can only trace rays on triangle meshes.");
  }

  // Build face and vertex arrays
  rawPositions.resize(mesh->nVertices() * 3);
  rawFaces.resize(mesh->nFaces() * 3);
  VertexData<size_t> vInd = mesh->getVertexIndices();
  for (Vertex v : mesh->vertices()) {
    unsigned int i = 3 * vInd[v];
    Vector3 p = geometry->position(v);
    for (unsigned int j = 0; j < 3; j++) rawPositions[i + j] = p[j];
  }
  FaceData<size_t> fInd = mesh->getFaceIndices();
  for (Face f : mesh->faces()) {
    unsigned int i = 3 * fInd[f];
    unsigned int j = 0;
    for (Vertex v : f.adjacentVertices()) {
      rawFaces[i + j] = vInd[v];
      j++;
    }
  }

  // Construct nanort mesh objects
  nanort::TriangleMesh<double> triangle_mesh(rawPositions.data(), rawFaces.data(), sizeof(double) * 3);
  nanort::TriangleSAHPred<double> triangle_pred(rawPositions.data(), rawFaces.data(),
                                                sizeof(double) * 3); // still have no idea what this does
  bool ret = accel.Build(mesh->nFaces(), options, triangle_mesh, triangle_pred);
  if (!ret) {
    throw std::runtime_error("BVH construction failed");
  }

  nanort::BVHBuildStatistics stats = accel.GetStatistics();

  cout << "BVH statistics:" << endl;
  cout << "    # of leaf   nodes: " << stats.num_leaf_nodes << endl;
  cout << "    # of branch nodes: " << stats.num_branch_nodes << endl;
  cout << "    Max tree depth   : " << stats.max_tree_depth << endl;

  double lengthScale = geometry->lengthScale();
  tFar = lengthScale * 1e3;
}

RayHitResult MeshRayTracer::trace(Vector3 start, Vector3 dir) {
  // Create the ray
  nanort::Ray<double> ray;
  ray.min_t = 0.0;
  ray.max_t = tFar;
  for (int i = 0; i < 3; i++) ray.org[i] = start[i];
  dir = unit(dir);
  for (int i = 0; i < 3; i++) ray.dir[i] = dir[i];

  // Compute the intersection
  nanort::BVHTraceOptions trace_options;
  nanort::TriangleIntersector<double> triangle_intersector(rawPositions.data(), rawFaces.data(), sizeof(double) * 3);
  bool hit = accel.Traverse(ray, trace_options, triangle_intersector);

  // Return the result
  if (hit) {
    RayHitResult result;
    result.hit = true;
    result.tHit = triangle_intersector.intersection.t;
    result.face = mesh->face(triangle_intersector.intersection.prim_id);

    // Convert barycentric formats
    double U = triangle_intersector.intersection.u;
    double V = triangle_intersector.intersection.v;
    result.baryCoords = Vector3{1.0 - U - V, U, V};

    return result;
  } else {
    return RayHitResult{false, std::numeric_limits<double>::quiet_NaN(), Face(), Vector3{-1.0, -1.0, -1.0}};
  }
}

} // namespace surface
}; // namespace geometrycentral
