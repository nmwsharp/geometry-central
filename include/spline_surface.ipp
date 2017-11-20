#pragma once

#include <utilities.h>
#include <spline_surface.h>
#include <typeinfo>
#include <string>

// Patch parameter =============================================================

inline PatchParam::PatchParam(const Vector3& p) : Vector3{p.x, p.y, p.z} {}

inline PatchParam::PatchParam(double a, double b, double c)
    : Vector3{a, b, c} {}

inline PatchParam::PatchParam(const Vector2& p) : Vector3{p.x, p.y, 0.} {}

inline PatchParam::PatchParam(double u, double v) : Vector3{u, v, 0.} {}

// Spline patches ==============================================================

template <typename T>
T SplinePatch<T>::evaluate(double u, double v) {
  throw std::invalid_argument("UV evaluation not allowed on a patch of type" +
                              typeNameString(this));
}

template <typename T>
T SplinePatch<T>::evaluate(double a, double b, double c) {
  throw std::invalid_argument(
      "Barycentric evaluation not allowed on a patch of type" +
      typeNameString(this));
}

// Bezier patches ==============================================================

template <typename T>
T BezierPatch<T>::evaluate(const PatchParam& point) {
  return evaluate(point.x, point.y);
}

template <typename T>
T BezierPatch<T>::evaluate(double u, double v) {
  double a = 1. - u;
  double b = 1. - v;

  // Evaluate cubic Bernstein polynomials in u
  double Bu0 = a * a * a;
  double Bu1 = 3. * u * a * a;
  double Bu2 = 3. * u * u * a;
  double Bu3 = u * u * u;

  // Evaluate cubic Bernstein polynomials in v
  double Bv0 = b * b * b;
  double Bv1 = 3. * v * b * b;
  double Bv2 = 3. * v * v * b;
  double Bv3 = v * v * v;

  // Evaluate bicubic interpolant at query point
  array4x4<T>& c = controlPoints;
  return Bv0 * (Bu0 * c[0][0] + Bu1 * c[1][0] + Bu2 * c[2][0] + Bu3 * c[3][0]) +
         Bv1 * (Bu0 * c[0][1] + Bu1 * c[1][1] + Bu2 * c[2][1] + Bu3 * c[3][1]) +
         Bv2 * (Bu0 * c[0][2] + Bu1 * c[1][2] + Bu2 * c[2][2] + Bu3 * c[3][2]) +
         Bv3 * (Bu0 * c[0][3] + Bu1 * c[1][3] + Bu2 * c[2][3] + Bu3 * c[3][3]);
}

template <typename T>
BezierPatchPtr<T> toBezierPatch(SplinePatchPtr<T> patch) {
  return std::dynamic_pointer_cast<BezierPatch<T>>(patch);
}

// GregoryPatchTri =============================================================

template <typename T>
T GregoryPatchTri<T>::evaluate(const PatchParam& point) {
  return evaluate(point.x, point.y, point.z);
}

template <typename T>
T GregoryPatchTri<T>::evaluate(double u, double v, double w) {
  // Evaluate face points, handling degeneracy at corners
  T F0 = (v + w == 0.) ? f0p : (w * f0m + v * f0p) / (v + w);
  T F1 = (w + u == 0.) ? f1p : (u * f1m + w * f1p) / (w + u);
  T F2 = (u + v == 0.) ? f2p : (v * f2m + u * f2p) / (u + v);

  // Evaluate triangular Gregory interpolant at query point
  return u * u * u * p0 + v * v * v * p1 + w * w * w * p2 +
         3. * u * v * (u + v) * (u * e0p + v * e1m) +
         3. * v * w * (v + w) * (v * e1p + w * e2m) +
         3. * w * u * (w + u) * (w * e2p + u * e0m) +
         12. * u * v * w * (u * F0 + v * F1 + w * F2);
}

template <typename T>
GregoryPatchTriPtr<T> toGregoryPatchTri(SplinePatchPtr<T> patch) {
  return std::dynamic_pointer_cast<GregoryPatchTri<T>>(patch);
}

// GregoryPatchQuad ============================================================

template <typename T>
T GregoryPatchQuad<T>::evaluate(const PatchParam& point) {
  return evaluate(point.x, point.y);
}

template <typename T>
T GregoryPatchQuad<T>::evaluate(double u, double v) {
  double a = 1. - u;
  double b = 1. - v;

  // Evaluate cubic Bernstein polynomials in u
  double Bu0 = a * a * a;
  double Bu1 = 3. * u * a * a;
  double Bu2 = 3. * u * u * a;
  double Bu3 = u * u * u;

  // Evaluate cubic Bernstein polynomials in v
  double Bv0 = b * b * b;
  double Bv1 = 3. * v * b * b;
  double Bv2 = 3. * v * v * b;
  double Bv3 = v * v * v;

  // Evaluate face points
  T F0 = (u + v == 0.) ? f0p : (u * f0p + v * f0m) / (u + v);
  T F1 = (a + v == 0.) ? f1p : (a * f1m + v * f1p) / (a + v);
  T F2 = (a + b == 0.) ? f2p : (a * f2p + b * f2m) / (a + b);
  T F3 = (u + b == 0.) ? f3p : (u * f3m + b * f3p) / (u + b);

  // Evaluate quadrilateral Gregory interpolant at query point
  return Bu0 * (Bv0 * p0 + Bv1 * e0m + Bv2 * e3p + Bv3 * p3) +
         Bu1 * (Bv0 * e0p + Bv1 * F0 + Bv2 * F3 + Bv3 * e3m) +
         Bu2 * (Bv0 * e1m + Bv1 * F1 + Bv2 * F2 + Bv3 * e2p) +
         Bu3 * (Bv0 * p1 + Bv1 * e1p + Bv2 * e2m + Bv3 * p2);
}

template <typename T>
GregoryPatchQuadPtr<T> toGregoryPatchQuad(SplinePatchPtr<T> patch) {
  return std::dynamic_pointer_cast<GregoryPatchQuad<T>>(patch);
}

// SplineSurface ===============================================================

template <typename T>
SplineSurface<T>::SplineSurface(Geometry<T>* geometry_)
    : mesh(*(geometry_->getMesh())), geometry(*geometry_) {}

template <typename T>
T SplineSurface<T>::evaluate(const PatchPoint& p) {
  SplinePatchPtr<T> patch = getPatch(p.face);

  return patch->evaluate(p.coords);
}

template <typename T>
SplinePatchPtr<T> SplineSurface<T>::getPatch(FacePtr f) {
  SplinePatchPtr<T> patch(nullptr);
  unsigned int degree = f.degree();

  switch (scheme) {
    case SplinePatchScheme::Bezier:
      if (degree == 4) {
        patch = buildBezierPatch(f);
      }
      break;
    case SplinePatchScheme::Gregory:
      if (degree == 3 || degree == 4) {
        patch = buildGregoryPatch(f);
      }
      break;
    case SplinePatchScheme::Pm:
      throw std::invalid_argument("Pm patch interpolation not yet implemented");
      // patch = buildPmPatch(f);
      break;
    case SplinePatchScheme::NURBS:
      throw std::invalid_argument("NURBS interpolation not yet implemented");
      // patch = buildNURBSPatch(f);
      break;
    default:
      throw std::invalid_argument(
          "No valid spline interpolation scheme specified");
      break;
  }

  return patch;
}

template <typename T>
SplinePatchPtr<T> SplineSurface<T>::buildBezierPatch(FacePtr f) {
  if (f.degree() != 4) {
    throw std::invalid_argument(
        "Bézier interpolation valid only for quadrilaterals");
  }

  SplinePatchPtr<T> patch;

  if (isRegular(f)) {
    patch = buildRegularBezierPatch(f);
  } else {
    patch = buildIrregularBezierPatch(f);
  }

  return patch;
}

template <typename T>
SplinePatchPtr<T> SplineSurface<T>::buildGregoryPatch(FacePtr f) {
  SplinePatchPtr<T> patch;

  if (f.degree() == 3) {
    patch = buildTriangularGregoryPatch(f);

  } else if (f.degree() == 4) {
    if (isRegular(f)) {
      patch = buildRegularBezierPatch(f);
    } else {
      patch = buildQuadrilateralGregoryPatch(f);
    }

  } else {
    throw std::invalid_argument(
        "Gregory interpolation valid only for triangles and quads");
  }

  return patch;
}

template <typename T>
SplinePatchPtr<T> SplineSurface<T>::buildRegularBezierPatch(FacePtr f) {
  // Assuming that f is contained in a regular patch (i.e., all four of its
  // corners are
  // degree 4), this method gathers the 16 vertex coordinates from the
  // surrounding 4x4
  // neighborhood.  By convention we label points as below, letting the halfedge
  // associated
  // with f go from c11 to c21.
  //
  // p03 --- p13 --- p23 --- p33
  //  |       |       |       |
  //  |       |       |       |
  // p02 --- p12 --- p22 --- p32
  //  |       |   f   |       |
  //  |       | ====> |       |
  // p01 --- p11 --- p21 --- p31
  //  |       |       |       |
  //  |       |       |       |
  // p00 --- p10 --- p20 --- p30
  //
  // The algorithm proceeds in two stages:
  //     I. Traverse the mesh, collecting vertex coordinates in an irregular
  //     order.
  //    II. Compute the corresponding Bézier coefficients
  // In particular, we first walk around the middle face (f), then around the
  // outer boundary.

  BezierPatch<T>* patch = new BezierPatch<T>();
  array4x4<T>& c = patch->controlPoints;

  // Get a halfedge of the middle face
  HalfedgePtr h = f.halfedge();

  // Walk around middle face, visiting p11, p21, p22, p12 (in that order)
  T p11 = geometry[h.vertex()];
  h = h.next();
  T p21 = geometry[h.vertex()];
  h = h.next();
  T p22 = geometry[h.vertex()];
  h = h.next();
  T p12 = geometry[h.vertex()];

  // Walk along left side, visiting p02, p01, p00
  h = h.twin().next().next();
  T p02 = geometry[h.vertex()];
  h = h.next();
  T p01 = geometry[h.vertex()];
  h = h.twin().next().next();
  T p00 = geometry[h.vertex()];

  // Walk along bottom, visiting p10, p20, p30
  h = h.next();
  T p10 = geometry[h.vertex()];
  h = h.twin().next().next();
  T p20 = geometry[h.vertex()];
  h = h.twin().next().next();
  T p30 = geometry[h.vertex()];

  // Walk along right side, visiting p31, 32, 33
  h = h.next();
  T p31 = geometry[h.vertex()];
  h = h.twin().next().next();
  T p32 = geometry[h.vertex()];
  h = h.twin().next().next();
  T p33 = geometry[h.vertex()];

  // Walk along top, visiting p23, p13, p03
  h = h.next();
  T p23 = geometry[h.vertex()];
  h = h.twin().next().next();
  T p13 = geometry[h.vertex()];
  h = h.twin().next().next();
  T p03 = geometry[h.vertex()];

  // Convert to Bézier form

  c[0][0] = (p02 + 4. * p12 + p22 + 4. * p01 + 16. * p11 + 4. * p21 + p00 +
             4. * p10 + p20) /
            36.;

  c[1][0] = (2. * p12 + p22 + 8. * p11 + 4. * p21 + 2. * p10 + p20) / 18.;

  c[2][0] = (p12 + 2. * p22 + 4. * p11 + 8. * p21 + p10 + 2. * p20) / 18.;

  c[3][0] = (p12 + 4. * p22 + p32 + 4. * p11 + 16. * p21 + 4. * p31 + p10 +
             4. * p20 + p30) /
            36.;

  // ---

  c[0][1] =
      (1. * p02 + 4. * p12 + 1. * p22 + 2. * p01 + 8. * p11 + 2. * p21) / 18.;

  c[1][1] = (2. * p12 + 1. * p22 + 4. * p11 + 2. * p21) / 9.;

  c[2][1] = (1. * p12 + 2. * p22 + 2. * p11 + 4. * p21) / 9.;

  c[3][1] =
      (1. * p12 + 4. * p22 + 1. * p32 + 2. * p11 + 8. * p21 + 2. * p31) / 18.;

  // ---

  c[0][2] =
      (2. * p02 + 8. * p12 + 2. * p22 + 1. * p01 + 4. * p11 + 1. * p21) / 18.;

  c[1][2] = (4. * p12 + 2. * p22 + 2. * p11 + 1. * p21) / 9.;

  c[2][2] = (2. * p12 + 4. * p22 + 1. * p11 + 2. * p21) / 9.;

  c[3][2] =
      (2. * p12 + 8. * p22 + 2. * p32 + 1. * p11 + 4. * p21 + 1. * p31) / 18.;

  // ---

  c[0][3] = (p03 + 4. * p13 + p23 + 4. * p02 + 16. * p12 + 4. * p22 + p01 +
             4. * p11 + p21) /
            36.;

  c[1][3] = (2. * p13 + p23 + 8. * p12 + 4. * p22 + 2. * p11 + p21) / 18.;

  c[2][3] = (p13 + 2. * p23 + 4. * p12 + 8. * p22 + p11 + 2. * p21) / 18.;

  c[3][3] = (p13 + 4. * p23 + p33 + 4. * p12 + 16. * p22 + 4. * p32 + p11 +
             4. * p21 + p31) /
            36.;

  return SplinePatchPtr<T>(patch);
}

template <typename T>
SplinePatchPtr<T> SplineSurface<T>::buildIrregularBezierPatch(FacePtr f) {
  BezierPatch<T>* patch = new BezierPatch<T>();
  array4x4<T>& c = patch->controlPoints;

  HalfedgePtr h11 = f.halfedge();
  HalfedgePtr h21 = h11.next();
  HalfedgePtr h22 = h21.next();
  HalfedgePtr h12 = h22.next();

  // Corner points
  c[0][0] = BezierControlPointVertex(h11.vertex());
  c[3][0] = BezierControlPointVertex(h21.vertex());
  c[0][3] = BezierControlPointVertex(h12.vertex());
  c[3][3] = BezierControlPointVertex(h22.vertex());

  // Face points
  c[1][1] = BezierControlPointFace(h11);
  c[2][1] = BezierControlPointFace(h21);
  c[1][2] = BezierControlPointFace(h12);
  c[2][2] = BezierControlPointFace(h22);

  // Edge points
  getBezierControlPointsEdge(h11, c[1][0], c[2][0]);
  getBezierControlPointsEdge(h12, c[0][2], c[0][1]);
  getBezierControlPointsEdge(h21, c[3][1], c[3][2]);
  getBezierControlPointsEdge(h22, c[2][3], c[1][3]);

  return SplinePatchPtr<T>(patch);
}

template <typename T>
T SplineSurface<T>::BezierControlPointVertex(VertexPtr v) {
  // TODO can probably be replaced with CatmullClarkLimitPosition()

  double n = v.degree();
  double weight = n * n;
  T sum = weight * geometry[v];

  for (HalfedgePtr h : v.outgoingHalfedges()) {
    sum += 4. * geometry[h.next().vertex()];
    weight += 4.;

    sum += 1. * geometry[h.next().next().vertex()];
    weight += 1.;
  }

  return sum / weight;
}

template <typename T>
T SplineSurface<T>::BezierControlPointFace(HalfedgePtr h) {
  double n = h.vertex().degree();

  T p00 = geometry[h.vertex()];
  h = h.next();
  T p10 = geometry[h.vertex()];
  h = h.next();
  T p11 = geometry[h.vertex()];
  h = h.next();
  T p01 = geometry[h.vertex()];

  return (n * p00 + 2. * p01 + 1. * p11 + 2. * p10) / (n + 5.);
}

template <typename T>
void SplineSurface<T>::getBezierControlPointsEdge(HalfedgePtr h, T& c0, T& c1) {
  double k0 = 2. * h.vertex().degree();
  T p11 = geometry[h.vertex()];
  h = h.next();
  double k1 = 2. * h.vertex().degree();
  T p21 = geometry[h.vertex()];
  h = h.next();
  T p22 = geometry[h.vertex()];
  h = h.next();
  T p12 = geometry[h.vertex()];
  h = h.next().twin().next().next();
  T p10 = geometry[h.vertex()];
  h = h.next();
  T p20 = geometry[h.vertex()];

  c0 = (2. * p12 + 1. * p22 + k0 * p11 + 4. * p21 + 2. * p10 + 1. * p20) /
       (k0 + 10.);

  c1 = (1. * p12 + 2. * p22 + 4. * p11 + k1 * p21 + 1. * p10 + 2. * p20) /
       (k1 + 10.);
}

template <typename T>
SplinePatchPtr<T> SplineSurface<T>::buildTriangularGregoryPatch(FacePtr f) {
  GregoryPatchTri<T>* patch = new GregoryPatchTri<T>();

  HalfedgePtr h0 = f.halfedge();
  HalfedgePtr h1 = h0.next();
  HalfedgePtr h2 = h1.next();

  // Corner points
  patch->p0 = CatmullClarkLimitPosition(h0.vertex());
  patch->p1 = CatmullClarkLimitPosition(h1.vertex());
  patch->p2 = CatmullClarkLimitPosition(h2.vertex());

  // Edge points
  patch->e0p = GregoryEdge(h0, patch->p0, false);
  patch->e1m = GregoryEdge(h0, patch->p1, true);
  patch->e1p = GregoryEdge(h1, patch->p1, false);
  patch->e2m = GregoryEdge(h1, patch->p2, true);
  patch->e2p = GregoryEdge(h2, patch->p2, false);
  patch->e0m = GregoryEdge(h2, patch->p0, true);

  // Face points
  patch->f0p = GregoryFace(h0, patch->p0, patch->e0p, patch->e1m, false);
  patch->f1p = GregoryFace(h1, patch->p1, patch->e1p, patch->e2m, false);
  patch->f2p = GregoryFace(h2, patch->p2, patch->e2p, patch->e0m, false);
  patch->f1m = GregoryFace(h0, patch->p1, patch->e1m, patch->e0p, true);
  patch->f2m = GregoryFace(h1, patch->p2, patch->e2m, patch->e1p, true);
  patch->f0m = GregoryFace(h2, patch->p0, patch->e0m, patch->e2p, true);

  return SplinePatchPtr<T>(patch);
}

template <typename T>
SplinePatchPtr<T> SplineSurface<T>::buildQuadrilateralGregoryPatch(FacePtr f) {
  GregoryPatchQuad<T>* patch = new GregoryPatchQuad<T>();

  HalfedgePtr h0 = f.halfedge();
  HalfedgePtr h1 = h0.next();
  HalfedgePtr h2 = h1.next();
  HalfedgePtr h3 = h2.next();

  // Corner points
  patch->p0 = CatmullClarkLimitPosition(h0.vertex());
  patch->p1 = CatmullClarkLimitPosition(h1.vertex());
  patch->p2 = CatmullClarkLimitPosition(h2.vertex());
  patch->p3 = CatmullClarkLimitPosition(h3.vertex());

  // Edge Points
  patch->e0p = GregoryEdge(h0, patch->p0, false);
  patch->e1m = GregoryEdge(h0, patch->p1, true);
  patch->e1p = GregoryEdge(h1, patch->p1, false);
  patch->e2m = GregoryEdge(h1, patch->p2, true);
  patch->e2p = GregoryEdge(h2, patch->p2, false);
  patch->e3m = GregoryEdge(h2, patch->p3, true);
  patch->e3p = GregoryEdge(h3, patch->p3, false);
  patch->e0m = GregoryEdge(h3, patch->p0, true);

  // Face Points
  patch->f0p = GregoryFace(h0, patch->p0, patch->e0p, patch->e1m, false);
  patch->f0m = GregoryFace(h3, patch->p0, patch->e0m, patch->e3p, true);
  patch->f1p = GregoryFace(h1, patch->p1, patch->e1p, patch->e2m, false);
  patch->f1m = GregoryFace(h0, patch->p1, patch->e1m, patch->e0p, true);
  patch->f2p = GregoryFace(h2, patch->p2, patch->e2p, patch->e3m, false);
  patch->f2m = GregoryFace(h1, patch->p2, patch->e2m, patch->e1p, true);
  patch->f3p = GregoryFace(h3, patch->p3, patch->e3p, patch->e0m, false);
  patch->f3m = GregoryFace(h2, patch->p3, patch->e3m, patch->e2p, true);

  return SplinePatchPtr<T>(patch);
}

template <typename T>
T SplineSurface<T>::CatmullClarkLimitPosition(VertexPtr v) {
  T p = geometry[v];
  T sum = T::zero();
  double n = 0.;

  for (HalfedgePtr h : v.outgoingHalfedges()) {
    T m = geometry.midpoint(h.edge());
    T c = geometry.barycenter(h.face());
    sum += m + c;
    n += 1.;
  }

  return ((n - 3.) / (n + 5.)) * p + (4. / (n * (n + 5.))) * sum;
}

template <typename T>
T SplineSurface<T>::GregoryEdge(HalfedgePtr h, const T& p, bool switched) {
  if (switched) {
    h = h.twin();
  }

  // Get degree of tail vertex
  const double n = h.vertex().degree();

  // Compute constants
  const double sigma = pow(4. + cos(M_PI / n), -.5);
  const double lambda =
      (1. / 16.) * (5. + cos(2. * M_PI / n) +
                    cos(M_PI / n) * sqrt(18. + 2. * cos(2. * M_PI / n)));

  T q = T::zero();
  double i = 0.;
  HalfedgePtr hi = h;
  do {
    T mi, ci;

    mi = geometry.midpoint(hi.edge());
    if (!switched) {
      ci = geometry.barycenter(hi.face());
    } else {
      ci = geometry.barycenter(hi.twin().face());
    }

    q += (1. - sigma * cos(M_PI / n)) * cos(2. * M_PI * i / n) * mi +
         2. * sigma * cos((2 * M_PI * i + M_PI) / n) * ci;

    hi = hi.twin().next();

    if (!switched) {
      i -= 1.;
    } else {
      i += 1.;
    }
  } while (hi != h);

  q *= (2. / n);

  return p + (2. / 3.) * lambda * q;
}

template <typename T>
T SplineSurface<T>::GregoryFace(HalfedgePtr h, const T& p0, const T& e0p,
                                const T& e1m, bool switched) {
  // Get degree of element being interpolated
  double d = h.face().degree() == 4 ? 3. : 4.;

  if (switched) {
    h = h.twin();
  }

  // Get vertex degrees at both halfedge endpoints
  double n0 = h.vertex().degree();
  double n1 = h.twin().vertex().degree();

  // Get midpoints of next and previous edges around h's vertex
  T mip1 = geometry.midpoint(h.prev().edge());
  T mim1 = geometry.midpoint(h.twin().next().edge());
  if (switched) std::swap(mip1, mim1);

  // Get the face barycenters on either side of h
  T ci = geometry.barycenter(h.face());
  T cim1 = geometry.barycenter(h.twin().face());
  if (switched) std::swap(ci, cim1);

  // Compute constants
  double c0 = cos(2. * M_PI / n0);
  double c1 = cos(2. * M_PI / n1);

  // Compute face control point
  T r0p = (1. / 3.) * (mip1 - mim1) + (2. / 3.) * (ci - cim1);
  T f0p = (1. / d) * (c1 * p0 + (d - 2. * c0 - c1) * e0p + 2. * c0 * e1m + r0p);

  return f0p;
}

template <typename T>
bool SplineSurface<T>::isRegular(FacePtr f) {
  for (VertexPtr v : f.adjacentVertices()) {
    if (v.degree() != 4) {
      return false;
    }

    for (FacePtr f : v.adjacentFaces()) {
      if (f.degree() != 4) {
        return false;
      }
    }
  }

  return true;
}

template <typename T>
HalfedgeMesh& SplineSurface<T>::getMesh() {
  return mesh;
}

template <typename T>
Geometry<T>& SplineSurface<T>::getGeometry() {
  return geometry;
}
