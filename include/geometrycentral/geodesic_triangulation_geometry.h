#pragma once

#include "geometrycentral/dependent_quantity.h"
#include "geometrycentral/halfedge_mesh.h"
#include "geometrycentral/intrinsic_geometry.h"
#include "geometrycentral/unit_vector3.h"
#include "geometrycentral/vector2.h"
#include "geometrycentral/vector3.h"

#include <Eigen/SparseCore>

#include <iostream>

namespace geometrycentral {


class GeodesicTriangulationGeometry : public IntrinsicGeometry {

public:
  GeodesicTriangulationGeometry(HalfedgeMesh* mesh_, EdgeData<double>& edgeLengths);
  GeodesicTriangulationGeometry(HalfedgeMesh* mesh_, VertexData<Vector3>& vertexPositions);


protected:
  std::vector<DependentQuantity*> allQuantities;

  // === Internal interface for all quantities

  DependentQuantity edgeLengthsQ;
  virtual void computeEdgeLengths() override;


private:
  EdgeData<double> geodesicEdgeLengths;
};

} // namespace geometrycentral
