#include "geometrycentral/surface/mesh_ray_tracer.h"

#include "nanort.h"

namespace geometrycentral {
namespace surface {

class MeshRayTracer::Impl {
public:
  Impl(EmbeddedGeometryInterface& geom) {
    geom.requireVertexPositions();

    std::vector<Vector3> positions;
    positions.reserve(geom.mesh.nVertices());
    for (Vertex v : geom.mesh.vertices()) {
      positions.push_back(geom.vertexPositions[v]);
    }

    std::vector<std::vector<size_t>> faces;
    faces.reserve(geom.mesh.nFaces());
    for (Face f : geom.mesh.faces()) {
      std::vector<size_t> tri;
      for (Vertex v : f.adjacentVertices()) tri.push_back(v.getIndex());
      if (tri.size() != 3) throw std::runtime_error("MeshRayTracer: mesh must be triangulated");
      faces.push_back(std::move(tri));
    }

    geom.unrequireVertexPositions();
    build(positions, faces);
  }

  Impl(const SimplePolygonMesh& mesh) { build(mesh.vertexCoordinates, mesh.polygons); }

  void build(const std::vector<Vector3>& positions, const std::vector<std::vector<size_t>>& faces) {
    double INF = std::numeric_limits<double>::infinity();
    Vector3 bboxMin{INF, INF, INF};
    Vector3 bboxMax{-INF, -INF, -INF};

    rawPositions.reserve(positions.size() * 3);
    for (const Vector3& p : positions) {
      rawPositions.push_back(p.x);
      rawPositions.push_back(p.y);
      rawPositions.push_back(p.z);
      bboxMin = componentwiseMin(bboxMin, p);
      bboxMax = componentwiseMax(bboxMax, p);
    }
    lengthScale = norm(bboxMax - bboxMin);

    rawFaces.reserve(faces.size() * 3);
    for (const std::vector<size_t>& f : faces) {
      if (f.size() != 3) throw std::runtime_error("MeshRayTracer: mesh must be triangulated");
      for (size_t idx : f) rawFaces.push_back(static_cast<unsigned int>(idx));
    }

    if (rawFaces.empty()) return; // empty mesh — all queries will miss

    nanort::TriangleMesh<double> nanortMesh(rawPositions.data(), rawFaces.data(), sizeof(double) * 3);
    nanort::TriangleSAHPred<double> nanortPred(rawPositions.data(), rawFaces.data(), sizeof(double) * 3);
    bool ok = accel.Build(static_cast<unsigned int>(rawFaces.size() / 3), nanortMesh, nanortPred);
    if (!ok) throw std::runtime_error("MeshRayTracer: BVH construction failed");
  }

  RayHitResult trace(Vector3 origin, Vector3 dir) const {
    if (rawFaces.empty()) return {};

    nanort::Ray<double> ray;
    ray.min_t = 1e-6 * lengthScale;
    ray.max_t = 1e1 * lengthScale;
    for (int i = 0; i < 3; i++) ray.org[i] = origin[i];
    for (int i = 0; i < 3; i++) ray.dir[i] = dir[i];

    nanort::TriangleIntersection<double> isect;
    nanort::TriangleIntersector<double> intersector(rawPositions.data(), rawFaces.data(), sizeof(double) * 3);
    bool hit = accel.Traverse(ray, intersector, &isect);

    if (!hit) return {};

    RayHitResult result;
    result.hit = true;
    result.t = isect.t;
    result.faceIndex = static_cast<size_t>(isect.prim_id);
    result.faceCoords = {1.0 - isect.u - isect.v, isect.u, isect.v};
    return result;
  }

  std::vector<double> rawPositions;
  std::vector<unsigned int> rawFaces;
  nanort::BVHAccel<double> accel;
  double lengthScale = 1.0;
};


MeshRayTracer::MeshRayTracer(EmbeddedGeometryInterface& geom) : impl(new Impl(geom)) {}
MeshRayTracer::MeshRayTracer(const SimplePolygonMesh& mesh) : impl(new Impl(mesh)) {}
MeshRayTracer::~MeshRayTracer() = default;

RayHitResult MeshRayTracer::trace(Vector3 origin, Vector3 dir) const { return impl->trace(origin, dir); }

SurfacePoint toSurfacePoint(const RayHitResult& hit, SurfaceMesh& mesh) {
  if (!hit.hit) throw std::runtime_error("toSurfacePoint: RayHitResult is not a hit");
  return SurfacePoint(mesh.face(hit.faceIndex), hit.faceCoords);
}

} // namespace surface
} // namespace geometrycentral
