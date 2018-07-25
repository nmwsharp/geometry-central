#pragma once

#include <iostream>

#include "geometrycentral/halfedge_mesh.h"
#include "geometrycentral/unit_vector3.h"
#include "geometrycentral/vector2.h"
#include "geometrycentral/vector3.h"

namespace geometrycentral {


class IntrinsicGeometry {

public:
  IntrinsicGeometry(HalfedgeMesh& mesh_, EdgeData<double>& edgeLengths);

  // members
  HalfedgeMesh& mesh;

protected:

private:
  EdgeData<double> edgeLengths;
};

} // namespace geometrycentral
