#include "geometrycentral/utilities/elementary_geometry.h"


#include "Eigen/Dense"

namespace geometrycentral {

bool inCircleTest(Vector2 pA, Vector2 pB, Vector2 pC, Vector2 pTest) {

  Eigen::Matrix4d A;
  // clang-format off
  A <<  pA.x, pA.y, norm2(pA), 1.,
        pB.x, pB.y, norm2(pB), 1.,
        pC.x, pC.y, norm2(pC), 1.,
        pTest.x, pTest.y, norm2(pTest), 1.;
  // clang-format on
  return A.determinant() > 0.;
}

}
