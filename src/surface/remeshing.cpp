#include "geometrycentral/surface/remeshing.h"

namespace geometrycentral {
namespace surface {

// The default trace options
const RemeshOptions defaultRemeshOptions;

void remesh(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, RemeshOptions options) {
  MutationManager mm(mesh, geom);
  return remesh(mesh, geom, mm, options);
}

void remesh(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, MutationManager& mm, RemeshOptions options) {
  if (options.targetEdgeLength < 0) {
    double meanLength = 0;
    geom.requireEdgeLengths();
    for (Edge e : mesh.edges()) {
      meanLength += geom.edgeLengths[e];
    }
    geom.unrequireEdgeLengths();
    meanLength /= mesh.nEdges();

    options.targetEdgeLength = abs(options.targetEdgeLength * meanLength);
  }

  bool doConnectivityChanges = true;

  for (size_t iIt = 0; iIt < options.maxIterations; iIt++) {
    if (doConnectivityChanges) {
      doConnectivityChanges = adjustEdgeLengths(mesh, geom, mm, options);
    }

    size_t nFlips = fixDelaunay(mesh, geom, mm);
    double flowDist;
    switch (options.smoothStyle) {
    case RemeshSmoothStyle::Circumcentric:
      flowDist = smoothByCircumcenter(mesh, geom, mm, 1, options.boundaryCondition);
      break;
    case RemeshSmoothStyle::Laplacian:
      flowDist = smoothByLaplacian(mesh, geom, mm, 1, options.boundaryCondition);
      break;
    }

    // std::cout << iIt << " : " << changedConnectivity << " " << nFlips << " " << flowDist << std::endl;
    if ((nFlips == 0) && (flowDist < 0.01)) break;
  }
  geom.refreshQuantities();
}

Vector3 vertexNormal(VertexPositionGeometry& geom, Vertex v, MutationManager& mm) {
  Vector3 totalNormal = Vector3::zero();
  for (Corner c : v.adjacentCorners()) {
    Vector3 cornerNormal = geom.cornerAngle(c) * geom.faceNormal(c.face());
    totalNormal += cornerNormal;
  }
  return normalize(totalNormal);
}

Vector3 boundaryVertexTangent(VertexPositionGeometry& geom, Vertex v, MutationManager& mm) {
  if (v.isBoundary()) {
    auto edgeVec = [&](Edge e) -> Vector3 {
      return (geom.vertexPositions[e.halfedge().tipVertex()] - geom.vertexPositions[e.halfedge().tailVertex()])
          .normalize();
    };

    Vector3 totalTangent = Vector3::zero();
    for (Edge e : v.adjacentEdges()) {
      if (e.isBoundary()) {
        totalTangent += edgeVec(e);
      }
    }
    return totalTangent.normalize();
  } else {
    return Vector3::zero();
  }
}

inline Vector3 projectToPlane(Vector3 v, Vector3 norm) { return v - norm * dot(norm, v); }
inline Vector3 projectToLine(Vector3 v, Vector3 tangent) { return tangent * dot(tangent, v); }

inline Vector3 edgeMidpoint(SurfaceMesh& mesh, VertexPositionGeometry& geom, Edge e) {
  Vector3 endPos1 = geom.vertexPositions[e.halfedge().tailVertex()];
  Vector3 endPos2 = geom.vertexPositions[e.halfedge().tipVertex()];
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

Vector3 findCircumcenter(VertexPositionGeometry& geom, Face f) {
  // retrieve the face's vertices
  int index = 0;
  Vector3 p[3];
  for (Vertex v0 : f.adjacentVertices()) {
    p[index] = geom.vertexPositions[v0];
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

bool isDelaunay(VertexPositionGeometry& geom, Edge e) {
  float angle1 = geom.cornerAngle(e.halfedge().next().next().corner());
  float angle2 = geom.cornerAngle(e.halfedge().twin().next().next().corner());
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

bool shouldCollapse(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, Edge e) {
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
  Vector3 midpoint = edgeMidpoint(mesh, geom, e);
  for (Halfedge he0 : edgesToCheck) {
    Halfedge heT = he0.twin();
    Vertex v1 = heT.tailVertex();
    Vertex v2 = heT.tipVertex();
    Vertex v3 = heT.next().tipVertex();
    Vector3 a = geom.vertexPositions[v1];
    Vector3 b = geom.vertexPositions[v2];
    Vector3 c = geom.vertexPositions[v3];
    if (checkFoldover(a, b, c, midpoint, 2)) {
      return false;
    }
  }
  return true;
}

// Warning: requires that geom.vertexDualAreas and geom.vertexMeanCurvatures are filled in with accurate data
double getSmoothMeanCurvature(VertexPositionGeometry& geom, Vertex v) {
  double A = geom.vertexDualAreas[v];
  double S = geom.vertexMeanCurvatures[v];
  double K = S / A;
  return K;
}

// flatLength: specifies how long the target edge length should be in flat regions
// curvatureAdaptation: controls how much variation in target length occurs due to curvature
double findMeanTargetL(SurfaceMesh& mesh, VertexPositionGeometry& geom, Edge e, double flatLength,
                       double curvatureAdaptation) {
  double averageH = 0;
  for (Vertex v : e.adjacentVertices()) {
    averageH += getSmoothMeanCurvature(geom, v);
  }
  averageH /= 2;
  double L = flatLength / (fabs(averageH) * curvatureAdaptation + 1);
  return L;
}


size_t fixDelaunay(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom) {
  MutationManager mm(mesh, geom);
  return fixDelaunay(mesh, geom, mm);
}

size_t fixDelaunay(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, MutationManager& mm) {
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

    if (e.isBoundary() || isDelaunay(geom, e)) continue;

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

double smoothByLaplacian(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, double stepSize,
                         RemeshBoundaryCondition bc) {
  MutationManager mm(mesh, geom);
  return smoothByLaplacian(mesh, geom, mm, stepSize, bc);
}

double smoothByLaplacian(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, MutationManager& mm, double stepSize,
                         RemeshBoundaryCondition bc) {
  VertexData<Vector3> vertexOffsets(mesh);
  for (Vertex v : mesh.vertices()) {
    // calculate average of surrounding vertices
    Vector3 avgNeighbor = Vector3::zero();
    for (Vertex j : v.adjacentVertices()) {
      avgNeighbor += geom.vertexPositions[j];
    }
    avgNeighbor /= v.degree();

    Vector3 updateDirection = avgNeighbor - geom.vertexPositions[v];

    // project updateDirection onto space of allowed movements
    Vector3 stepDir;
    if (v.isBoundary()) {
      switch (bc) {
      case RemeshBoundaryCondition::Fixed:
        stepDir = Vector3::zero();
        break;
      case RemeshBoundaryCondition::Tangential:
        // for free boundary vertices, project the average to the boundary tangent line
        stepDir = projectToLine(updateDirection, boundaryVertexTangent(geom, v, mm));
        break;
      case RemeshBoundaryCondition::Free:
        // for free boundary vertices, project the average to the surface tangent plane
        stepDir = projectToPlane(updateDirection, vertexNormal(geom, v, mm));
        break;
      }
    } else {
      // for interior vertices, project the average to the tangent plane
      stepDir = projectToPlane(updateDirection, vertexNormal(geom, v, mm));
    }
    vertexOffsets[v] = stepSize * stepDir;
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

double smoothByCircumcenter(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, double stepSize,
                            RemeshBoundaryCondition bc) {
  MutationManager mm(mesh, geom);
  return smoothByCircumcenter(mesh, geom, mm, stepSize, bc);
}

double smoothByCircumcenter(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, MutationManager& mm,
                            double stepSize, RemeshBoundaryCondition bc) {
  geom.requireFaceAreas();
  VertexData<Vector3> vertexOffsets(mesh);
  for (Vertex v : mesh.vertices()) {
    Vector3 updateDirection = Vector3::zero();
    for (Face f : v.adjacentFaces()) {
      // add the circumcenter weighted by face area to the update direction
      Vector3 circum = findODTCenter(geom, f, mm);
      updateDirection += geom.faceArea(f) * (circum - geom.vertexPositions[v]);
    }
    updateDirection /= (6 * geom.vertexDualArea(v));

    // project updateDirection onto space of allowed movements
    Vector3 stepDir;
    if (v.isBoundary()) {
      switch (bc) {
      case RemeshBoundaryCondition::Fixed:
        stepDir = Vector3::zero();
        break;
      case RemeshBoundaryCondition::Tangential:
        // for free boundary vertices, project the average to the boundary tangent line
        stepDir = projectToLine(updateDirection, boundaryVertexTangent(geom, v, mm));
        break;
      case RemeshBoundaryCondition::Free:
        // for free boundary vertices, project the average to the surface tangent plane
        stepDir = projectToPlane(updateDirection, vertexNormal(geom, v, mm));
        break;
      }
    } else {
      // for interior vertices, project the average to the tangent plane
      stepDir = projectToPlane(updateDirection, vertexNormal(geom, v, mm));
    }
    vertexOffsets[v] = stepSize * stepDir;
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
} // namespace surface


bool adjustEdgeLengths(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, RemeshOptions options) {
  MutationManager mm(mesh, geom);
  return adjustEdgeLengths(mesh, geom, mm, options);
}

bool adjustEdgeLengths(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, MutationManager& mm,
                       RemeshOptions options) {
  geom.requireVertexDualAreas();
  geom.requireVertexMeanCurvatures();

  double curvatureAdaptation = options.curvatureAdaptation;
  bool useCurvatureAdaptation = curvatureAdaptation > 0;
  double flatLength = options.targetEdgeLength;
  double minLength = options.minRelativeLength * options.targetEdgeLength;

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
    double length_e = geom.edgeLength(e);
    double threshold =
        (useCurvatureAdaptation) ? findMeanTargetL(mesh, geom, e, flatLength, curvatureAdaptation) : flatLength;
    if (length_e > minLength && length_e > threshold * 1.5) {
      Vector3 newPos = edgeMidpoint(mesh, geom, e);
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
        (useCurvatureAdaptation) ? findMeanTargetL(mesh, geom, e, flatLength, curvatureAdaptation) : flatLength;
    if (geom.edgeLength(e) < threshold * 0.5) {
      Vector3 newPos = edgeMidpoint(mesh, geom, e);
      if (shouldCollapse(mesh, geom, e)) {
        Vertex v = mm.collapseEdge(e, newPos);
        if (v != Vertex()) {
          didSplitOrCollapse = true;
        }
      }
    }
  }

  geom.unrequireVertexDualAreas();
  geom.unrequireVertexMeanCurvatures();

  mesh.compress();
  return didSplitOrCollapse;
}

} // namespace surface
} // namespace geometrycentral
