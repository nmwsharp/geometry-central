#pragma once

#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/utilities/vector3.h"

#include <iostream>
#include <limits>
#include <sstream>

namespace geometrycentral {
namespace surface {

// Represents a generic point on the surface of a triangle mesh, which could be a
// vertex, a point along an edge, or a point inside a face.
// NOTE: implicitly assumes triangluar without checks

enum class SurfacePointType { Vertex = 0, Edge, Face };

struct SurfacePoint {

  // === Constructors
  SurfacePoint();                           // default: yields invalid SurfacePoint with Type::Vertex and null vertex
  SurfacePoint(Vertex v);                   // at vertex
  SurfacePoint(Edge e, double tEdge);       // in edge
  SurfacePoint(Face f, Vector3 faceCoords); // in face


  // === Identifying data
  SurfacePointType type = SurfacePointType::Vertex;

  // if vertex
  Vertex vertex = Vertex();

  // if edge
  Edge edge = Edge();
  double tEdge = std::numeric_limits<double>::quiet_NaN();

  // if face
  Face face = Face();
  Vector3 faceCoords = Vector3::undefined();


  // === Methods

  // All surface points (vertex, edge, face) have an equivalent point in one or many adjacent faces. This function
  // returns one of the equivalent surface points in a face (chosen arbitrarily). If this point is a face point, the
  // output is a copy of this point.
  inline SurfacePoint inSomeFace() const;
  
  // Return the nearest vertex to this surface point
  inline Vertex nearestVertex() const;


  // Linearly interpolate data at vertices to this point.
  // T must support addition and multiplication by a double.
  template <typename T>
  inline T interpolate(const VertexData<T>& data) const;


  // Throws an exception if the surface point is invalid in any way
  inline void validate() const;


  // === Operators
  bool operator==(const SurfacePoint& other) const;
  bool operator!=(const SurfacePoint& other) const;
};

// Printing
::std::ostream& operator<<(std::ostream& output, const SurfacePoint& p);

} // namespace surface
} // namespace geometrycentral

namespace std {
std::string to_string(geometrycentral::surface::SurfacePoint p);
}

#include "geometrycentral/surface/surface_point.ipp"

