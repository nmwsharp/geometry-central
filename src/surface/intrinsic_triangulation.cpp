#include "geometrycentral/surface/intrinsic_triangulation.h"


namespace geometrycentral {
namespace surface {

IntrinsicTriangulation::IntrinsicTriangulation(ManifoldSurfaceMesh& mesh_, IntrinsicGeometryInterface& inputGeom_)
    : EdgeLengthGeometry(*mesh_.copy().release()), inputMesh(mesh_), inputGeom(inputGeom_),
      intrinsicMesh(dynamic_cast<ManifoldSurfaceMesh*>(&mesh)) {


  // do this here, rather than in the constructer, since we need to call require() first
  inputGeom.requireEdgeLengths();
  edgeLengths = inputGeom.edgeLengths;

  // Make sure the input mesh is triangular
  if (!mesh.isTriangular()) {
    throw std::runtime_error("intrinsic triangulation requires triangle mesh as input");
  }

  // Initialize vertex locations
  vertexLocations = VertexData<SurfacePoint>(mesh);
  for (size_t iV = 0; iV < mesh.nVertices(); iV++) {
    vertexLocations[iV] = SurfacePoint(inputMesh.vertex(iV));
  }

  // == Register the default callback which maintains marked edges
  auto updateMarkedEdges = [&](Edge oldE, Halfedge newHe1, Halfedge newHe2) {
    if (markedEdges.size() > 0 && markedEdges[oldE]) {
      markedEdges[newHe1.edge()] = true;
      markedEdges[newHe2.edge()] = true;
    }
  };
  edgeSplitCallbackList.push_back(updateMarkedEdges);
}

// ======================================================
// ======== Queries & Accessors
// ======================================================


EdgeData<std::vector<SurfacePoint>> IntrinsicTriangulation::traceEdges() {

  EdgeData<std::vector<SurfacePoint>> tracedEdges(mesh);

  for (Edge e : mesh.edges()) {
    Halfedge he = e.halfedge();
    tracedEdges[e] = traceHalfedge(he, false);
  }

  return tracedEdges;
}

bool IntrinsicTriangulation::isDelaunay(double delaunayEPS) {
  for (Edge e : mesh.edges()) {
    if (!isDelaunay(e)) {
      return false;
    }
  }
  return true;
}

bool IntrinsicTriangulation::isDelaunay(Edge e, double delaunayEPS) {
  if (!isFixed(e) && edgeCotanWeight(e) < -delaunayEPS) {
    return false;
  }
  return true;
}

double IntrinsicTriangulation::minAngleDegrees() {
  double minAngle = std::numeric_limits<double>::infinity();
  for (Corner c : mesh.corners()) {
    minAngle = std::min(minAngle, cornerAngle(c));
  }
  return minAngle * 180. / M_PI;
}

// ======================================================
// ======== Mutators
// ======================================================

bool SignpostIntrinsicTriangulation::flipEdgeIfNotDelaunay(Edge e, double delaunayEPS) {

  // Can't flip
  if (isFixed(e)) return false;

  // Don't want to flip
  double cWeight = edgeCotanWeight(e);
  if (cWeight > -delaunayEPS) return false;

  // Get geometric data
  Halfedge he = e.halfedge();
  std::array<Vector2, 4> layoutPositions = layoutDiamond(he);

  // Combinatorial flip
  bool flipped = intrinsicMesh->flip(e, false);

  // Should always be possible, something unusual is going on if we end up here
  if (!flipped) {
    return false;
  }

  // Compute the new edge length
  double newLength = (layoutPositions[1] - layoutPositions[3]).norm();

  // If we're going to create a non-finite edge length, abort the flip
  // (only happens if you're in a bad numerical place)
  if (!std::isfinite(newLength)) {
    intrinsicMesh->flip(e, false);
    return false;
  }

  // Assign the new edge lengths
  edgeLengths[e] = newLength;

  // do any extra work required by subclass
  flipEdgeInternal(e, newLength);

  invokeEdgeFlipCallbacks(e);
  return true;
}

bool IntrinsicTriangulation::flipEdgeIfPossible(Edge e, double possibleEPS) {

  // Can't flip
  if (isFixed(e)) return false;

  // Get geometric data
  Halfedge he = e.halfedge();
  std::array<Vector2, 4> layoutPositions = layoutDiamond(he);

  // Test if geometryically flippable flippable (both signed areas of new triangles are positive)
  double A1 = cross(layoutPositions[1] - layoutPositions[0], layoutPositions[3] - layoutPositions[0]);
  double A2 = cross(layoutPositions[3] - layoutPositions[2], layoutPositions[1] - layoutPositions[2]);
  double areaEPS = possibleEPS * (A1 + A2);
  if (A1 < areaEPS || A2 < areaEPS) {
    return false;
  }


  // Combinatorial flip
  bool flipped = intrinsicMesh->flip(e, false);

  // Might not have been flippable for connectivity reasons
  if (!flipped) {
    return false;
  }

  // Compute the new edge length
  double newLength = (layoutPositions[1] - layoutPositions[3]).norm();

  // If we're going to create a non-finite edge length, abort the flip
  // (only happens if you're in a bad numerical place)
  if (!std::isfinite(newLength)) {
    intrinsicMesh->flip(e, false);
    return false;
  }

  // Assign the new edge lengths
  // TODO project to satisfy triangle inequality?
  edgeLengths[e] = newLength;

  // do any extra work required by subclass
  flipEdgeInternal(e, newLength);

  invokeEdgeFlipCallbacks(e);
  return true;
}


} // namespace surface
} // namespace geometrycentral
