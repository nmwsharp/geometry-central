#pragma once

#include "geometrycentral/geometry.h"

#include <array>
#include <unordered_set>

namespace geometrycentral {

// Forward declarations
struct MutableHalfedge;

struct MutableVertex {

  MutableVertex(Vector3 p_);

  Vector3 position;
  MutableHalfedge* halfedge = nullptr;
  size_t degree = 0;
  bool isBoundary = false;

  // Not always populated
  size_t index;
  double desiredVertexDensity; // used for adaptive remeshing
  bool touched = false;

  Vector3 normal();
};

struct MutableFace {

  MutableHalfedge* halfedge = nullptr; // NOTE: Not guaranteed to be real FORNOW, unlike HalfedgeMesh
                                       // (can be done, but not implemented)

  Vector3 normal();
  Vector3 areaNormal();
  Vector3 barycenter();
  Vector3 circumcenter();
  double desiredDensity();
};

struct MutableHalfedge {

  MutableHalfedge* twin = nullptr;
  MutableHalfedge* next = nullptr;
  MutableVertex* vertex = nullptr;
  MutableFace* face = nullptr;

  bool isOnBoundary();
  bool isReal();

  Vector3 vector();
  double oppositeAngle();
};


class MutableMesh {

public:
  MutableMesh(Geometry<Euclidean>* geometry);
  ~MutableMesh();

  // Parameters
  double normalThreshold = PI / 2.0; // how much a modification is allowed to change a face normal
  double skinnyThresh = 0.01;        // relative threshold for deciding if a face is skinny

  // Modifications
  bool tryFlipEdge(MutableHalfedge* he1);
  bool tryCollapseEdge(MutableHalfedge* he1);
  bool trySplitEdge(MutableHalfedge* he1);
  bool tryMoveVertex(MutableVertex* v, Vector3 newPos);

  Vector3 splineMidpoint(Vector3 p1, Vector3 n1, Vector3 p2, Vector3 n2);

  Geometry<Euclidean>* toHalfedgeMesh();
  void indexVertices();
  void computeSizingField(double absoluteEpsilon);
  void setVerticesUntouched();

  void checkDegree();

  // Members
  std::unordered_set<MutableVertex*> vertices;
  std::unordered_set<MutableFace*> faces;
  std::unordered_set<MutableHalfedge*> halfedges;
};

} // namespace geometrycentral
