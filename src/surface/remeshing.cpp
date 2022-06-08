#include "geometrycentral/surface/remeshing.h"

namespace geometrycentral {
namespace surface {

Vector3 vertexNormal(VertexPositionGeometry& geometry, Vertex v, MutationManager& mm) {
  Vector3 totalNormal = Vector3::zero();
  Vector3 fixedFaceNormals = Vector3::zero();
  bool foundFixedFace = false;
  for (Corner c : v.adjacentCorners()) {
    size_t nFixedVertices = 0;
    for (Vertex w : c.face().adjacentVertices()) {
      if (!mm.mayRepositionVertex(w)) nFixedVertices++;
    }

    Vector3 cornerNormal = geometry.cornerAngle(c) * geometry.faceNormal(c.face());
    if (nFixedVertices >= 2) {
      fixedFaceNormals += cornerNormal;
      foundFixedFace = true;
    }
    totalNormal += cornerNormal;
  }
  if (foundFixedFace) {
    return normalize(fixedFaceNormals);
  } else {
    return normalize(totalNormal);
  }
}

inline Vector3 projectToPlane(Vector3 v, Vector3 norm) { return v - norm * dot(norm, v); }

inline Vector3 edgeMidpoint(SurfaceMesh& mesh, VertexPositionGeometry& geometry, Edge e) {
  Vector3 endPos1 = geometry.vertexPositions[e.halfedge().tailVertex()];
  Vector3 endPos2 = geometry.vertexPositions[e.halfedge().tipVertex()];
  return (endPos1 + endPos2) / 2;
}

Vector3 findCircumcenter(Vector3 p1, Vector3 p2, Vector3 p3) {
  // barycentric coordinates of circumcenter
  double a = (p3 - p2).norm();
  double b = (p3 - p1).norm();
  double c = (p2 - p1).norm();
  double a2 = a * a;
  double b2 = b * b;
  double c2 = c * c;
  Vector3 circumcenterLoc{a2 * (b2 + c2 - a2), b2 * (c2 + a2 - b2), c2 * (a2 + b2 - c2)};
  // normalize to sum of 1
  circumcenterLoc = normalizeBarycentric(circumcenterLoc);

  // change back to space
  return circumcenterLoc[0] * p1 + circumcenterLoc[1] * p2 + circumcenterLoc[2] * p3;
}

Vector3 findCircumcenter(VertexPositionGeometry& geometry, Face f) {
  // retrieve the face's vertices
  int index = 0;
  Vector3 p[3];
  for (Vertex v0 : f.adjacentVertices()) {
    p[index] = geometry.vertexPositions[v0];
    index++;
  }
  return findCircumcenter(p[0], p[1], p[2]);
}

// Returns the barycenter for faces incident on a nonflippable edge (e.g. a boundary edge), and the circumcenter for all
// other faces
Vector3 findODTCenter(VertexPositionGeometry& geom, Face f, MutationManager& mm) {
  Vector3 p0 = geom.vertexPositions[f.halfedge().tailVertex()];
  Vector3 p1 = geom.vertexPositions[f.halfedge().tipVertex()];
  Vector3 p2 = geom.vertexPositions[f.halfedge().next().tipVertex()];

  for (Edge e : f.adjacentEdges()) {
    if (e.isBoundary() || !mm.mayFlipEdge(e)) {
      // e is not flippable. return barycenter
      return (p0 + p1 + p2) / 3.;
    }
  }
  return findCircumcenter(p0, p1, p2);
}

bool isDelaunay(VertexPositionGeometry& geometry, Edge e) {
  float angle1 = geometry.cornerAngle(e.halfedge().next().next().corner());
  float angle2 = geometry.cornerAngle(e.halfedge().twin().next().next().corner());
  return angle1 + angle2 <= PI;
}

inline double diamondAngle(Vector3 a, Vector3 b, Vector3 c, Vector3 d) // dihedral angle at edge a-b
{
  Vector3 n1 = cross(b - a, c - a);
  Vector3 n2 = cross(b - d, a - d);
  return PI - angle(n1, n2);
}

inline bool checkFoldover(Vector3 a, Vector3 b, Vector3 c, Vector3 x, double angle) {
  return diamondAngle(a, b, c, x) < angle;
}

bool shouldCollapse(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geometry, Edge e) {
  std::vector<Halfedge> edgesToCheck;
  Vertex v1 = e.halfedge().vertex();
  Vertex v2 = e.halfedge().twin().vertex();

  // find (halfedge) link around the edge, starting with those surrounding v1
  for (Halfedge he : v1.outgoingHalfedges()) {
    if (he.next().tailVertex() != v2 && he.next().tipVertex() != v2) {
      edgesToCheck.push_back(he.next());
    }
  }

  // link around v2
  for (Halfedge he : v2.outgoingHalfedges()) {
    if (he.next().tailVertex() != v1 && he.next().tipVertex() != v1) {
      edgesToCheck.push_back(he.next());
    }
  }

  // see if the point that would form after a collapse would cause a major foldover with surrounding edges
  Vector3 midpoint = edgeMidpoint(mesh, geometry, e);
  for (Halfedge he0 : edgesToCheck) {
    Halfedge heT = he0.twin();
    Vertex v1 = heT.tailVertex();
    Vertex v2 = heT.tipVertex();
    Vertex v3 = heT.next().tipVertex();
    Vector3 a = geometry.vertexPositions[v1];
    Vector3 b = geometry.vertexPositions[v2];
    Vector3 c = geometry.vertexPositions[v3];
    if (checkFoldover(a, b, c, midpoint, 2)) {
      return false;
    }
  }
  return true;
}

// Warning: requires that geometry.vertexDualAreas and geometry.vertexMeanCurvatures are filled in with accurate data
double getSmoothMeanCurvature(VertexPositionGeometry& geometry, Vertex v) {
  double A = geometry.vertexDualAreas[v];
  double S = geometry.vertexMeanCurvatures[v];
  double K = S / A;
  return K;
}

// flatLength: specifies how long the target edge length should be in flat regions
// curvatureAdaptation: controls how much variation in target length occurs due to curvature
double findMeanTargetL(SurfaceMesh& mesh, VertexPositionGeometry& geometry, Edge e, double flatLength,
                       double curvatureAdaptation) {
  double averageH = 0;
  for (Vertex v : e.adjacentVertices()) {
    averageH += getSmoothMeanCurvature(geometry, v);
  }
  averageH /= 2;
  double L = flatLength / (fabs(averageH) * curvatureAdaptation + 1);
  return L;
}


size_t fixDelaunay(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geometry) {
  MutationManager mm(mesh, geometry);
  return fixDelaunay(mesh, geometry, mm);
}

size_t fixDelaunay(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geometry, MutationManager& mm) {
  // Logic duplicated from surface/intrinsic_triangulation.cpp

  std::deque<Edge> edgesToCheck;      // queue of edges to check if Delaunay
  EdgeData<bool> inQueue(mesh, true); // true if edge is currently in edgesToCheck

  // start with all edges
  for (Edge e : mesh.edges()) {
    edgesToCheck.push_back(e);
  }

  // counter and limit for number of flips
  size_t flipMax = 100 * mesh.nVertices();
  size_t nFlips = 0;
  while (!edgesToCheck.empty() && nFlips < flipMax) {
    Edge e = edgesToCheck.front();
    edgesToCheck.pop_front();
    inQueue[e] = false;

    if (e.isBoundary() || isDelaunay(geometry, e)) continue;

    // if not Delaunay, try to flip edge
    bool wasFlipped = mm.flipEdge(e);

    if (!wasFlipped) continue;

    nFlips++;

    // Add neighbors to queue, as they may need flipping now
    Halfedge he = e.halfedge();
    std::array<Edge, 4> neighboringEdges{he.next().edge(), he.next().next().edge(), he.twin().next().edge(),
                                         he.twin().next().next().edge()};
    for (Edge nE : neighboringEdges) {
      if (!inQueue[nE]) {
        edgesToCheck.push_back(nE);
        inQueue[nE] = true;
      }
    }
  }
  return nFlips;
}

double smoothByLaplacian(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geometry, double stepSize) {
  MutationManager mm(mesh, geometry);
  return smoothByLaplacian(mesh, geometry, mm, stepSize);
}

double smoothByLaplacian(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geometry, MutationManager& mm,
                         double stepSize) {
  VertexData<Vector3> vertexOffsets(mesh);
  for (Vertex v : mesh.vertices()) {
    if (v.isBoundary()) {
      vertexOffsets[v] = Vector3::zero();
    } else {
      // calculate average of surrounding vertices
      Vector3 avgNeighbor = Vector3::zero();
      for (Vertex j : v.adjacentVertices()) {
        avgNeighbor += geometry.vertexPositions[j];
      }
      avgNeighbor /= v.degree();

      // and project the average to the tangent plane
      Vector3 updateDirection = avgNeighbor - geometry.vertexPositions[v];
      Vector3 normal = vertexNormal(geometry, v, mm);
      vertexOffsets[v] = stepSize * projectToPlane(updateDirection, normal);
    }
  }

  // update final vertices
  double totalMovement = 0;
  for (Vertex v : mesh.vertices()) {
    bool didMove = mm.repositionVertex(v, vertexOffsets[v]);
    if (didMove) {
      totalMovement += vertexOffsets[v].norm();
    }
  }
  return totalMovement / mesh.nVertices();
}

double smoothByCircumcenter(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geometry, double stepSize) {
  MutationManager mm(mesh, geometry);
  return smoothByCircumcenter(mesh, geometry, mm, stepSize);
}

double smoothByCircumcenter(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geometry, MutationManager& mm,
                            double stepSize) {
  geometry.requireFaceAreas();
  VertexData<Vector3> vertexOffsets(mesh);
  for (Vertex v : mesh.vertices()) {
    if (v.isBoundary()) {
      vertexOffsets[v] = Vector3::zero();
    } else {
      vertexOffsets[v] = Vector3::zero();
      for (Face f : v.adjacentFaces()) {
        // add the circumcenter weighted by face area to the update direction
        Vector3 circum = findODTCenter(geometry, f, mm);
        vertexOffsets[v] += geometry.faceArea(f) * (circum - geometry.vertexPositions[v]);
      }
      vertexOffsets[v] /= (3 * geometry.vertexDualArea(v));
      // project update direction to tangent plane
      vertexOffsets[v] = stepSize * projectToPlane(vertexOffsets[v], vertexNormal(geometry, v, mm)) / 2.;
    }
  }

  // update final vertices
  double totalMovement = 0;
  for (Vertex v : mesh.vertices()) {
    bool didMove = mm.repositionVertex(v, vertexOffsets[v]);
    if (didMove) {
      totalMovement += vertexOffsets[v].norm();
    }
  }
  return totalMovement / mesh.nVertices();
}


bool adjustEdgeLengths(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geometry, double flatLength,
                       double curvatureAdaptation, double minRelativeLength) {
  MutationManager mm(mesh, geometry);
  return adjustEdgeLengths(mesh, geometry, mm, flatLength, curvatureAdaptation, minRelativeLength);
}

bool adjustEdgeLengths(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geometry, MutationManager& mm,
                       double flatLength, double curvatureAdaptation, double minRelativeLength) {
  geometry.requireVertexDualAreas();
  geometry.requireVertexMeanCurvatures();

  bool useCurvatureAdaptation = curvatureAdaptation > 0;
  double minLength = minRelativeLength * flatLength;

  bool didSplitOrCollapse = false;
  // queues of edges to CHECK to change
  std::vector<Edge> toSplit;
  std::vector<Edge> toCollapse;

  for (Edge e : mesh.edges()) {
    toSplit.push_back(e);
  }

  // actually splitting
  while (!toSplit.empty()) {
    Edge e = toSplit.back();
    toSplit.pop_back();
    double length_e = geometry.edgeLength(e);
    double threshold =
        (useCurvatureAdaptation) ? findMeanTargetL(mesh, geometry, e, flatLength, curvatureAdaptation) : flatLength;
    if (length_e > minLength && length_e > threshold * 1.5) {
      Vector3 newPos = edgeMidpoint(mesh, geometry, e);
      Halfedge he = mm.splitEdge(e, newPos);
      if (he != Halfedge()) {
        didSplitOrCollapse = true;
      }
    } else {
      toCollapse.push_back(e);
    }
  }

  // actually collapsing
  while (!toCollapse.empty()) {
    Edge e = toCollapse.back();
    toCollapse.pop_back();
    if (e == Edge() || e.isDead()) continue; // make sure it exists

    double threshold =
        (useCurvatureAdaptation) ? findMeanTargetL(mesh, geometry, e, flatLength, curvatureAdaptation) : flatLength;
    if (geometry.edgeLength(e) < threshold * 0.5) {
      Vector3 newPos = edgeMidpoint(mesh, geometry, e);
      if (shouldCollapse(mesh, geometry, e)) {
        Vertex v = mm.collapseEdge(e, newPos);
        if (v != Vertex()) {
          didSplitOrCollapse = true;
        }
      }
    }
  }

  geometry.unrequireVertexDualAreas();
  geometry.unrequireVertexMeanCurvatures();

  mesh.compress();
  return didSplitOrCollapse;
}

void remesh(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geometry, double targetEdgeLength, size_t maxIterations,
            double curvatureAdaptation, double minRelativeLength) {
  MutationManager mm(mesh, geometry);
  return remesh(mesh, geometry, mm, targetEdgeLength, maxIterations, curvatureAdaptation, minRelativeLength);
}

void remesh(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geometry, MutationManager& mm, double targetEdgeLength,
            size_t maxIterations, double curvatureAdaptation, double minRelativeLength) {

  bool doConnectivityChanges = true;

  for (size_t iIt = 0; iIt < maxIterations; iIt++) {
    if (doConnectivityChanges) {
      doConnectivityChanges =
          adjustEdgeLengths(mesh, geometry, mm, targetEdgeLength, curvatureAdaptation, minRelativeLength);
    }

    size_t nFlips = fixDelaunay(mesh, geometry, mm);
    double flowDist = smoothByCircumcenter(mesh, geometry, mm);

    // std::cout << iIt << " : " << changedConnectivity << " " << nFlips << " " << flowDist << std::endl;
    if ((nFlips == 0) && (flowDist < 0.01)) break;
  }
  geometry.refreshQuantities();
}

} // namespace surface
} // namespace geometrycentral
