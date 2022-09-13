#pragma once

#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/surface_point.h"

namespace geometrycentral {
namespace surface {

// A BarycentricVector represents a barycentric vector on the surface of a triangle mesh, which is the difference of two
// barycentric points (SurfacePoints.)

enum class BarycentricVectorType { Face = 0, Edge, Vertex };

struct BarycentricVector {

  // === Constructors
  BarycentricVector();       // default: yields invalid vector with null Face
  BarycentricVector(Face f); // yields face-type BarycentricVector with zero vector
  BarycentricVector(SurfacePoint pA, SurfacePoint pB);
  BarycentricVector(Vertex v);
  BarycentricVector(Edge e, Vector2 edgeCoords);
  BarycentricVector(Face f, Vector3 faceCoords);

  BarycentricVectorType type = BarycentricVectorType::Face;

  // if face-type
  Face face = Face();
  Vector3 faceCoords = Vector3::undefined(); // must sum to zero

  // if edge-type
  Edge edge = Edge();
  Vector2 edgeCoords = Vector2::undefined(); // must sum to zero

  // if vertex-type
  Vertex vertex = Vertex();

  // === Methods

  // All BarycentricVectors (vertex, edge, face) have an equivalent vector in one or many adjacent faces. This function
  // returns one of the equivalent BarycentricVectors in a face (chosen arbitrarily). If this vector is a face vector,
  // the output is a copy of this vector.
  BarycentricVector inSomeFace() const;

  // Returns the BarycentricVector as a face vector in face f. If the the BarycentricVector is not on or adjacent to the
  // requested face, throws an error.
  BarycentricVector inFace(Face f) const;

  // Returns the BarycentricVector as an edge vector in edge e. If this vector is already an edge vector, the output is
  // a copy of the vector. If the the BarycentricVector is not on the edge or one of its adjacent vertices, throws an
  // error.
  BarycentricVector inEdge(Edge e) const;

  void validate() const;

  // Shift values so the components sum to zero.
  void normalizeDisplacement();                     // modify in-place
  BarycentricVector normalizedDisplacement() const; // return a new vector

  // Compute norms; must define the geometry w.r.t. which norms are evaluated.
  double norm(IntrinsicGeometryInterface& geom) const;
  double norm2(IntrinsicGeometryInterface& geom) const;

  // Overloaded operators
  BarycentricVector operator+(const BarycentricVector& v) const;
  BarycentricVector operator-(const BarycentricVector& v) const;
  BarycentricVector operator*(double s) const;
  BarycentricVector operator/(double s) const;
  BarycentricVector& operator+=(const BarycentricVector& other);
  BarycentricVector& operator-=(const BarycentricVector& other);
  BarycentricVector& operator*=(const double& s);
  BarycentricVector& operator/=(const double& s);
  bool operator==(const BarycentricVector& other) const;
  bool operator!=(const BarycentricVector& other) const;
  const BarycentricVector operator-() const;

  // Rotate the vector 90 degrees CCW within the face it belongs to.
  // This requires the geometry, since the meaning of "90 degrees" depends on the geometry of the triangle.
  BarycentricVector rotated90(IntrinsicGeometryInterface& geom) const; // return a new vector
};

// Returns an arbitrary face shared by two vectors, if one exists; returns Face() if none.
Face sharedFace(const BarycentricVector& u, const BarycentricVector& v);
// Returns the edge shared by two vectors, if one exists; return Edge() if none.
Edge sharedEdge(const BarycentricVector& u, const BarycentricVector& v);

// Norms & inner products (require geometry)
double norm(IntrinsicGeometryInterface& geom, const BarycentricVector& v);
double norm2(IntrinsicGeometryInterface& geom, const BarycentricVector& v);
double dot(IntrinsicGeometryInterface& geom, const BarycentricVector& u, const BarycentricVector& v);

// Scalar multiplication
template <typename T>
BarycentricVector operator*(const T s, const BarycentricVector& v);

// double angle(const BarycentricVector& u, const BarycentricVector& v);

// Printing
::std::ostream& operator<<(::std::ostream& output, const BarycentricVector& v);

} // namespace surface
} // namespace geometrycentral

namespace std {
std::string to_string(geometrycentral::surface::BarycentricVector v);
}

#include "geometrycentral/surface/barycentric_vector.ipp"