#pragma once

// The Geometry class specifies the geometry of a given mesh, i.e.,
// the location of its vertices in space, where "in space" could mean
// a variety of things (in 3D, in the plane, on the sphere).  It is
// also used to answer queries about that geometry, e.g., what is the
// area of a given triangle or the length of a given edge.
//
// There are two principal ways to evaluate a given geometric quantity:
//
//    1. via methods that evaluate attributes directly, or
//    2. by caching attributes in a Data vector (Geometry::get*())
//
// The former is useful in scenarios where the geometry is constantly
// changing (e.g., when smoothing or editing a surface), hence cached
// values quickly become stale and there is little gained by storing
// them.  The latter is useful in scenarios where the geometry remains
// fixed throughout the entire algorithm (e.g., flattening or
// remeshing), hence there is no reason to repeatedly recompute values.
// Caching may also simplify syntax and reduce function call overhead.
// Finally, caching only the attributes that are needed (rather than
// all possible attributes) reduces memory usage.
//
// Example usage:
//
//    // compute the total surface area directly
//    double sum = 0.;
//    for(FacePtr f : mesh->faces)
//    {
//       sum += geometry->area(f);
//    }
//
//    // compute the total surface area using cached values
//    FaceData<double> area;
//    geometry->getFaceAreas(area); // fill the cache
//    double sum = 0.;
//    for(FacePtr f : mesh->faces())
//    {
//       sum += area[f];
//    }
//
// Note that the Geometry class also includes convenience methods for
// computing standard quantities like total surface area, etc.
//

// TODO: Right now many parts of Geometry implicitly assume that all real faces
// are triangular.
//       If you want to store non-triangular faces, be careful.

#include <iostream>

#include "geometrycentral/halfedge_mesh.h"
#include "geometrycentral/sparse_matrix.h"
#include "geometrycentral/unit_vector3.h"
#include "geometrycentral/vector2.h"
#include "geometrycentral/vector3.h"

namespace geometrycentral {

// Possible geometry types
typedef Vector2 Planar;  // TODO change to Complex
typedef Vector3 Euclidean;
typedef UnitVector3 Spherical;
// TODO Hyperbolic

enum class DualType { Barycentric, Circumcentric, Incentric, PowerDiagram };

// TODO In the future, could be extended to other types of mesh data structures
// (e.g., via an additional template argument)
template <class T>
class Geometry : public VertexData<T> {
  // TODO add flags for dual type (circumcentric, barycentric, incentric...)
  // TODO add subclass WeightedGeometry that has additional value (weight) per
  // vertex

 public:
  Geometry(HalfedgeMesh& mesh_)
      : VertexData<T>(&mesh_), mesh(mesh_), dualMesh(mesh_), p(*this) {}

  HalfedgeMesh* getMesh(void);  // Returns a pointer to the domain

  DualType dualType =
      DualType::Barycentric;  // determines geometry of dual cells

  // Vertex attributes
  // --- Primal ---
  T& position(VertexPtr p);  // TODO get rid of this method; all write access
                             // should be done through operator[], to
                             // distinguish it from all other named accessors
                             // (which are read-only)
  T position(VertexPtr p) const;
  double volume(VertexPtr v);       // always equal to 1
  double angleDefect(VertexPtr v);  // 2π minus sum of incident angles
  Vector3 normal(
      VertexPtr v);  // area-weighted average of incident face normals
  Vector3 boundaryNormal(VertexPtr v);  // length-weighted normal vector to the
                                        // two neighboring edges
  Vector3 projectToTangentSpace(VertexPtr v, const Vector3& inVec);
  Complex tangentVectorToComplexAngle(VertexPtr v, const Vector3& inVec);
  Vector3 complexAngleToTangentVector(VertexPtr v, Complex inAngle);
  Complex principalDirection(
      VertexPtr v);  // the 2-symmetric complex vector aligned with k1
  // --- Dual ---
  T position(DualVertexPtr p);
  double volume(DualVertexPtr v);

  // Edge attributes
  // --- Primal ---
  T midpoint(EdgePtr e);
  double length(EdgePtr e);
  double cotanWeight(EdgePtr e,
                     bool faceAreaWeighted = false);  // **triangles only**
  double dihedralAngle(EdgePtr e);
  // --- Dual ---
  double length(DualEdgePtr e);  // **triangles only**

  // Face attributes
  // --- Primal ---
  double area(FacePtr f);
  Vector3 normal(FacePtr f);
  Vector3 areaVector(FacePtr f);
  T barycenter(FacePtr f);
  T circumcenter(FacePtr f);
  // --- Dual ---
  double area(DualFacePtr f);

  // Halfedge attributes
  T vector(HalfedgePtr h);
  double angle(HalfedgePtr h);              // **triangles only**
  double angle(CornerPtr c);                // **triangles only**
  double angularCoordinate(HalfedgePtr h);  // **triangles only** Measured CCW
                                            // against the tail vertex's
                                            // arbitrary halfedge
  double cotan(HalfedgePtr h);              // **triangles only**

  // Global attributes
  double totalArea(void);  // Total surface area (assuming all triangles)
  T center(void);          // Center of mass (assuming constant density)
  void boundingBox(T& bboxMin,
                   T& bboxMax);  // Corners of axis-aligned bounding box
  T extent(void);  // Width, height, and depth of the axis-aligned bounding box
  double lengthScale(void);  // A length scale for the geometry

  // Methods for caching current attributes
  void getVertexPositions(VertexData<T>& vertexPosition);
  void getVertexNormals(VertexData<Vector3>& vertexNormal);
  void getVertexAngleDefects(VertexData<double>& vertexAngleDefect);
  void getPrincipalDirections(VertexData<Complex>& principalDirections);
  void getPrincipalDirections(VertexData<Complex>& principalDirections,
                              HalfedgeData<double>& angularCoordinates);
  void getDualVertexPositions(DualVertexData<T>& dualVertexPosition);

  void getEdgeLengths(EdgeData<double>& edgeLength);
  void getDualEdgeLengths(DualEdgeData<double>& dualEdgeLength);
  void getEdgeCotanWeights(EdgeData<double>& edgeCotanWeight,
                           bool faceAreaWeighted = false);

  void getFaceAreas(FaceData<double>& faceArea);
  void getDualFaceAreas(DualFaceData<double>& dualFaceArea);
  void getFaceNormals(FaceData<Vector3>& faceNormal);
  void getFaceBarycenters(FaceData<T>& faceBarycenter);

  void getHalfedgeVectors(HalfedgeData<T>& halfedgeVector);
  void getHalfedgeAngles(HalfedgeData<double>& halfedgeAngle);
  void getCornerAngles(CornerData<double>& cornerAngle);
  void getHalfedgeCotans(HalfedgeData<double>& halfedgeCotan);

  void getAngularCoordinates(HalfedgeData<double>& angularCoordinates);

  void normalize();

  // members
  HalfedgeMesh& mesh;
  HalfedgeDual dualMesh;

 protected:
  VertexData<T>& p;  // convenience reference to "this"
};

template <typename T, typename G>
geometrycentral::SparseMatrix<T> cotanMatrix(
    Geometry<G>* geometry, VertexData<size_t> index,
    bool faceAreaWeighted = false);  // usual symmetric cotan matrix (note this
                                     // is **not** the full Laplacian---for that
                                     // you need a mass matrix (below))
template <typename T, typename G>
geometrycentral::SparseMatrix<T> cotanMatrix(Geometry<G>* geometry,
                                             CornerData<size_t> index,
                                             size_t nC,
                                             bool faceAreaWeighted = false);

template <typename T, typename G>
geometrycentral::SparseMatrix<T> vertexMassMatrix(
    Geometry<G>* geometry,
    VertexData<size_t> index);  // diagonal lumped mass/finite volume mass
                                // matrix, where the geometry of dual cells is
                                // determined by Geometry::dualType
template <typename T, typename G>
geometrycentral::SparseMatrix<T> faceMassMatrix(
    Geometry<G>* geometry,
    FaceData<size_t> index);  // diagonal mass matrix containing face areas

}  // namespace geometrycentral

#include "geometrycentral/geometry.ipp"
