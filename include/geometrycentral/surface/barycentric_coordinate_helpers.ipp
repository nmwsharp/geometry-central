#pragma once


namespace geometrycentral {
namespace surface {

inline double displacementLength2(Vector3 displacement, Vector3 triangleLengths) {

  double distanceSq = 0.;
  for (int i = 0; i < 3; i++) {
    distanceSq += -triangleLengths[i] * triangleLengths[i] * displacement[(i + 1) % 3] * displacement[(i + 2) % 3];
  }

  return distanceSq;
}


inline double displacementLength(Vector3 displacement, Vector3 triangleLengths) {
  double length2 = std::fmax(displacementLength2(displacement, triangleLengths), 0.);
  return std::sqrt(length2);
}

inline Vector2 barycentricDisplacementToCartesian(const std::array<Vector2, 3>& vertCoords, Vector3 baryVec) {
  // Note: happens to be the same as baryCoordsToFaceCoords()
  return vertCoords[0] * baryVec.x + vertCoords[1] * baryVec.y + vertCoords[2] * baryVec.z;
}

inline Vector3 cartesianVectorToBarycentric(const std::array<Vector2, 3>& vertCoords, Vector2 faceVec, bool verbose) {


  // Build matrix for linear transform problem
  // (last constraint comes from chosing the displacement vector with sum = 0)
  Eigen::Matrix3d A;
  Eigen::Vector3d rhs;
  const std::array<Vector2, 3>& c = vertCoords; // short name
  A << c[0].x, c[1].x, c[2].x, c[0].y, c[1].y, c[2].y, 1., 1., 1.;
  rhs << faceVec.x, faceVec.y, 0.;

  // Solve
  Eigen::Vector3d result = A.colPivHouseholderQr().solve(rhs);
  Vector3 resultBary{result(0), result(1), result(2)};

  resultBary = normalizeBarycentricDisplacement(resultBary);

  return resultBary;
}

// Allows not-normalized input
inline Vector3 normalizeBarycentric(Vector3 baryCoords) {
  double s = sum(baryCoords);
  if (s == 0.) {
    return Vector3::constant(1. / 3.);
  }
  return baryCoords / s;
}

inline Vector3 normalizeBarycentricDisplacement(Vector3 baryVec) {
  double s = sum(baryVec);
  return baryVec - Vector3::constant(s / 3.);
}

// Allows not-normalized input
inline Vector3 projectInsideTriangle(Vector3 baryCoords) {
  for (int i = 0; i < 3; i++) {
    baryCoords[i] = std::fmax(baryCoords[i], 0.);
  }
  return normalizeBarycentric(baryCoords);
}

inline bool isInsideTriangle(Vector3 baryCoords) { return baryCoords.x >= 0 && baryCoords.y >= 0 && baryCoords.z >= 0; }

// The index of the halfedge in a triangular face,
inline int halfedgeIndexInTriangle(Halfedge he) {
  Halfedge heTest = he.face().halfedge();
  if (he == heTest) return 0;
  heTest = heTest.next();
  if (he == heTest) return 1;
  heTest = heTest.next();
  if (he == heTest) return 2;

  throw std::runtime_error("called halfedgeIndexInTriangle on non-triangular face");
  return -7777777;
}


// Given barycentric coordinates defined by treating refHe.vertex() as the `x` coordinate, permute to the canonical
// coordinate ordering for the face, which has face.halfedge().vertex() as the `x` coordinate.
inline Vector3 permuteBarycentricToCanonical(Vector3 baryCoords, Halfedge refHe) {
  int heInd = halfedgeIndexInTriangle(refHe);
  Vector3 result;
  for (int i = 0; i < 3; i++) {
    result[(i + heInd) % 3] = baryCoords[i];
  }
  return result;
}

inline Vector3 permuteBarycentricFromCanonical(Vector3 baryCoords, Halfedge refHe) {
  int heInd = halfedgeIndexInTriangle(refHe);
  Vector3 result;
  for (int i = 0; i < 3; i++) {
    result[i] = baryCoords[(i + heInd) % 3];
  }
  return result;
}


} // namespace surface
} // namespace geometrycentral
