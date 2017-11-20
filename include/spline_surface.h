// spline_surface.h
//
// A SplineSurface approximates a given polygonal mesh by a collection of
// continuous spline patches.  These patches can be evaluated at a given point
// to obtain various geometric quantities (position, normal, tangents, etc.).
// Several different spline patch schemes are available.
//
// The basic object is a SplineSurface, which is constructed from an existing
// Geometry instance via
//
//    SplineSurface<> surface( &geometry );
//
// From here, one can query the SplineSurface in two different ways:
//
//    1. by requesting the patch associated with a given Face, or
//    2. by directly evaluating the surface at a given point.
//
// To get a patch one can call the method SplineSurface::getPatch.  For example:
//
//    FacePtr f = mesh.face(0); // grab the first face in the mesh
//    SplinePatchPtr<> patch = surface.getPatch( f ); // get the spline patch
//    associated w/ f
//
// From here, one can evaluate any point within the patch.  Points are specified
// using a PatchParam instance, which stores either the parametric (u,v)
// coordinates
// for a quadrilateral patch, or the barycentric (u,v,w) coordinates for a
// triangular
// patch.  For instance, to evaluate the patch at its center, one might continue
// the
// code above as
//
//    PatchParam p;
//    if( f.degree() == 3 )
//    {
//       p = PatchParam( 1./3., 1./3., 1./3. );
//    }
//    else if( f.degree() == 4 )
//    {
//       p = PatchParam p( 1./2., 1.2. );
//    }
//    center = patch->evaluate( p );
//

#pragma once

#include <geometry.h>
#include <memory>
#include <array>

// A PatchParam specifies a point on a spline patch via its parametric
// coordinates.  Currently two domain geometries are supported: triangular and
// quadrilateral.  The interpretation of a patch parameter depends on the type
// of patch it is passed to: for triangular patches, all three coordinates are
// used as barycentric coordinates; for a quadrilateral patch, only the first
// two coordinates are used as standard parametric "uv" coordinates.
class PatchParam : public Vector3 {
 public:
  PatchParam();

  // Tri (barycentric)
  PatchParam(const Vector3& p);
  PatchParam(double a, double b, double c);

  // Quad (parametric)
  PatchParam(const Vector2& p);
  PatchParam(double u, double v);
};

// A SplinePatch is an abstract base class defining the interface for
// spline patches used in various patching schemes.  A SplinePatch can
// be evaluated at any point (specified by a PatchParam) to obtain the
// interpolated value.  In the special-but-common case where the patch
// interpolates points in 3D, one can also evaluate geometry quantities
// such as tangents, normals, etc.
//
// Sometimes it is useful to determine the particular patch type, e.g.,
// in order to draw a patch.  Some convenience methods are defined for
// doing type inference.  For example, to determine whether a given
// patch is a Bézier patch, one can write
//
//    SplinePatchPtr<> patch = myObject.getPatch();
//    BezierPatchPtr<> bPatch = toBezierPatch( patch );
//    if( bPatch != nullptr )
//    {
//       std::cout << "This patch is a Bézier patch." << std::endl;
//    }
//
//
template <typename T = Vector3>
class SplinePatch {
 public:
  virtual T evaluate(const PatchParam& point) = 0;
  virtual T evaluate(double u, double v);
  virtual T evaluate(double a, double b, double c);

  // Differential geometry (can be evaluated only when T=Vector3)
  virtual Vector3 position(const PatchParam& point) = 0;
  virtual Vector3 tangentU(const PatchParam& point) = 0;
  virtual Vector3 tangentV(const PatchParam& point) = 0;
  virtual Vector3 normal(const PatchParam& point) = 0;
  // virtual double meanCurvature(const PatchParam& point) = 0;
  // virtual double GaussianCurvature(const PatchParam& point) = 0;
  // virtual Vector3 pushForward(const PatchParam& point, const Vector2&
  // tangent) = 0;
  // virtual void principalDirections(const PatchParam& point, Vector2& X1,
  // Vector2& X2) = 0;
  // virtual void principalCurvatures(const PatchParam& point, double& kappa1,
  // double& kappa2) = 0;
  // virtual Vector2 parallelTransport(const PatchParam& sourcePoint, const
  // PatchParam& targetPoint, Vector2& tangent) = 0;
  // virtual PatchParam exponentialMap(const PatchParam& sourcePoint, const
  // Vector2& tangent) = 0;
};
template <typename T = Vector3>
using SplinePatchPtr = std::shared_ptr<SplinePatch<T>>;

template <typename T>
using array4x4 = std::array<std::array<T, 4>, 4>;

template <typename T = Vector3>
class BezierPatch : public SplinePatch<T> {
 public:
  virtual T evaluate(const PatchParam& point) override;
  virtual T evaluate(double u, double v) override;

  virtual Vector3 position(const PatchParam& point) override;
  virtual Vector3 tangentU(const PatchParam& point) override;
  virtual Vector3 tangentV(const PatchParam& point) override;
  virtual Vector3 normal(const PatchParam& point) override;
  // virtual double meanCurvature(const PatchParam& point) override;
  // virtual double GaussianCurvature(const PatchParam& point) override;
  // virtual Vector3 pushForward(const PatchParam& point, Vector2 tangent)
  // override;
  // virtual void principalDirections(const PatchParam& point, Vector2& X1,
  // Vector2& X2) override;
  // virtual void principalCurvatures(const PatchParam& point, double& kappa1,
  // double& kappa2) override;
  // virtual Vector2 parallelTransport(const PatchParam& sourcePoint, const
  // PatchParam& targetPoint, Vector2& tangent) override;
  // virtual PatchParam exponentialMap(const PatchParam& sourcePoint, Vector2&
  // tangent) override;

  array4x4<T> controlPoints;
};
template <typename T = Vector3>
using BezierPatchPtr = std::shared_ptr<BezierPatch<T>>;
template <typename T = Vector3>
BezierPatchPtr<T> toBezierPatch(SplinePatchPtr<T> patch);

template <typename T = Vector3>
class GregoryPatchTri : public SplinePatch<T> {
 public:
  virtual T evaluate(const PatchParam& point) override;
  virtual T evaluate(double a, double b, double c) override;

  virtual Vector3 position(const PatchParam& point) override;
  virtual Vector3 tangentU(const PatchParam& point) override;
  virtual Vector3 tangentV(const PatchParam& point) override;
  virtual Vector3 normal(const PatchParam& point) override;
  // virtual double meanCurvature(const PatchParam& point) override;
  // virtual double GaussianCurvature(const PatchParam& point) override;
  // virtual Vector3 pushForward(const PatchParam& point, Vector2 tangent)
  // override;
  // virtual void principalDirections(const PatchParam& point, Vector2& X1,
  // Vector2& X2) override;
  // virtual void principalCurvatures(const PatchParam& point, double& kappa1,
  // double& kappa2) override;
  // virtual Vector2 parallelTransport(const PatchParam& sourcePoint, const
  // PatchParam& targetPoint, Vector2& tangent) override;
  // virtual PatchParam exponentialMap(const PatchParam& sourcePoint, Vector2&
  // tangent) override;

  T p0, p1, p2;                    // Corner points
  T e0p, e1m, e1p, e2m, e2p, e0m;  // Edge points
  T f0p, f1p, f2p, f1m, f2m, f0m;  // Face points
};
template <typename T = Vector3>
using GregoryPatchTriPtr = std::shared_ptr<GregoryPatchTri<T>>;
template <typename T = Vector3>
GregoryPatchTriPtr<T> toGregoryPatchTri(SplinePatchPtr<T> patch);

template <typename T = Vector3>
class GregoryPatchQuad : public SplinePatch<T> {
 public:
  virtual T evaluate(const PatchParam& coords) override;
  virtual T evaluate(double u, double v) override;

  virtual Vector3 position(const PatchParam& point) override;
  virtual Vector3 tangentU(const PatchParam& point) override;
  virtual Vector3 tangentV(const PatchParam& point) override;
  virtual Vector3 normal(const PatchParam& point) override;
  // virtual double meanCurvature(const PatchParam& point) override;
  // virtual double GaussianCurvature(const PatchParam& point) override;
  // virtual Vector3 pushForward(const PatchParam& point, Vector2 tangent)
  // override;
  // virtual void principalDirections(const PatchParam& point, Vector2& X1,
  // Vector2& X2) override;
  // virtual void principalCurvatures(const PatchParam& point, double& kappa1,
  // double& kappa2) override;
  // virtual Vector2 parallelTransport(const PatchParam& sourcePoint, const
  // PatchParam& targetPoint, Vector2& tangent) override;
  // virtual PatchParam exponentialMap(const PatchParam& sourcePoint, Vector2&
  // tangent) override;

  T p0, p1, p2, p3;                          // Corner points
  T e0p, e1m, e1p, e2m, e2p, e3m, e3p, e0m;  // Edge points
  T f0p, f0m, f1p, f1m, f2p, f2m, f3p, f3m;  // Face points
};
template <typename T = Vector3>
using GregoryPatchQuadPtr = std::shared_ptr<GregoryPatchQuad<T>>;
template <typename T = Vector3>
GregoryPatchQuadPtr<T> toGregoryPatchQuad(SplinePatchPtr<T> patch);

class PatchPoint {
 public:
  FacePtr face;
  PatchParam coords;
};

class SplineTangent {
 public:
  FacePtr face;
  PatchParam coords;
  Vector2 vector;
};

enum class SplinePatchScheme {
  Bezier,
  Gregory,
  Pm,    // (not yet implemented)
  NURBS  // (not yet implemented)
};

template <typename T = Vector3>
class SplineSurface {
 public:
  SplineSurface(Geometry<T>* geometry);

  SplinePatchScheme scheme = SplinePatchScheme::Gregory;

  T evaluate(const PatchPoint& p);
  SplinePatchPtr<T> getPatch(FacePtr f);

  // Differential geometry (can be evaluated only when T=Vector3)
  Vector3 position(const PatchPoint& point);
  Vector3 tangentU(const PatchPoint& point);
  Vector3 tangentV(const PatchPoint& point);
  Vector3 normal(const PatchPoint& point);
  double meanCurvature(const PatchPoint& point);
  double GaussianCurvature(const PatchPoint& point);
  Vector3 pushForward(SplineTangent& tangent);
  void principalDirections(const PatchPoint& point, Vector2& X1, Vector2& X2);
  void principalCurvatures(const PatchPoint& point, double& kappa1,
                           double& kappa2);
  SplineTangent parallelTransport(const PatchPoint& sourcePoint,
                                  const PatchPoint& targetPoint,
                                  Vector2& tangent);
  PatchPoint exponentialMap(const PatchPoint& sourcePoint, Vector2& tangent);

  HalfedgeMesh& getMesh();
  Geometry<T>& getGeometry();

 protected:
  // Bézier
  SplinePatchPtr<T> buildBezierPatch(FacePtr f);
  SplinePatchPtr<T> buildRegularBezierPatch(FacePtr f);
  SplinePatchPtr<T> buildIrregularBezierPatch(FacePtr f);
  T BezierControlPointVertex(VertexPtr v);
  T BezierControlPointFace(HalfedgePtr h);
  void getBezierControlPointsEdge(HalfedgePtr h, T& c0, T& c1);

  // Gregory
  SplinePatchPtr<T> buildGregoryPatch(FacePtr f);
  SplinePatchPtr<T> buildTriangularGregoryPatch(FacePtr f);
  SplinePatchPtr<T> buildQuadrilateralGregoryPatch(FacePtr f);
  T CatmullClarkLimitPosition(VertexPtr v);
  T GregoryEdge(HalfedgePtr h, const T& p, bool switched);
  T GregoryFace(HalfedgePtr h, const T& p0, const T& e0p, const T& e1m,
                bool switched);

  bool isRegular(FacePtr f);
  HalfedgeMesh& mesh;
  Geometry<T>& geometry;
};

// Implementation
#include "spline_surface.ipp"
