#pragma once

#include "geometrycentral/surface/edge_length_geometry.h"

namespace geometrycentral {
namespace surface {

// Face areas
inline double EdgeLengthGeometry::faceArea(Face f) const {
  // WARNING: Logic duplicated between cached and immediate version
  Halfedge he = f.halfedge();
  double a = edgeLengths[he.edge()];
  he = he.next();
  double b = edgeLengths[he.edge()];
  he = he.next();
  double c = edgeLengths[he.edge()];

  GC_SAFETY_ASSERT(he.next() == f.halfedge(), "faces must be triangular");

  // Heron's formula
  double s = (a + b + c) / 2.0;
  double arg = s * (s - a) * (s - b) * (s - c);
  arg = std::fmax(0., arg);
  double area = std::sqrt(arg);
  return area;
}

// Vertex dual areas
inline double EdgeLengthGeometry::vertexDualArea(Vertex v) const {
  // WARNING: Logic duplicated between cached and immediate version
   double area = 0.;
   for( Face f : v.adjacentFaces() ) {
      area += faceArea(f);
   }
   return area/3.;
}

// Corner angles
inline double EdgeLengthGeometry::cornerAngle(Corner c) const {
  // WARNING: Logic duplicated between cached and immediate version
  Halfedge heA = c.halfedge();
  Halfedge heOpp = heA.next();
  Halfedge heB = heOpp.next();

  GC_SAFETY_ASSERT(heB.next() == heA, "faces must be triangular");

  double lOpp = edgeLengths[heOpp.edge()];
  double lA = edgeLengths[heA.edge()];
  double lB = edgeLengths[heB.edge()];

  double q = (lA * lA + lB * lB - lOpp * lOpp) / (2. * lA * lB);
  q = clamp(q, -1.0, 1.0);
  double angle = std::acos(q);

  return angle;
}


inline double EdgeLengthGeometry::halfedgeCotanWeight(Halfedge heI) const {
  // WARNING: Logic duplicated between cached and immediate version
  double cotSum = 0.;

  if (heI.isInterior()) {
    Halfedge he = heI;
    double l_ij = edgeLengths[he.edge()];
    he = he.next();
    double l_jk = edgeLengths[he.edge()];
    he = he.next();
    double l_ki = edgeLengths[he.edge()];
    he = he.next();
    GC_SAFETY_ASSERT(he == heI, "faces must be triangular");
    double area = faceArea(he.face());
    double cotValue = (-l_ij * l_ij + l_jk * l_jk + l_ki * l_ki) / (4. * area);
    return cotValue / 2;
  } else {
    return 0.;
  }
}

inline double EdgeLengthGeometry::vertexGaussianCurvature(Vertex v) const {
  // WARNING: Logic duplicated between cached and immediate version

   // the triangles neighboring any boundary vertex can be flattened into
   // the plane without any stretching/distortion; hence, a boundary
   // vertex has no Gaussian curvature
   if( v.isBoundary() ) return 0.;

   double gaussianCurvature = 2.*PI;
   for (Corner c : v.adjacentCorners() ) {
      gaussianCurvature -= cornerAngle(c);
   }
   return gaussianCurvature;
}

inline double EdgeLengthGeometry::edgeCotanWeight(Edge e) const {
  double sum = 0;
  for (Halfedge he : e.adjacentInteriorHalfedges()) {
    sum += halfedgeCotanWeight(he);
  }
  return sum;
}

inline double EdgeLengthGeometry::faceCircumradius(Face f) const {
  Halfedge he = f.halfedge();
  double a = edgeLengths[he.edge()];
  he = he.next();
  double b = edgeLengths[he.edge()];
  he = he.next();
  double c = edgeLengths[he.edge()];

  double A = faceArea(f);
  return a * b * c / (4. * A);
}


} // namespace surface
} // namespace geometrycentral
