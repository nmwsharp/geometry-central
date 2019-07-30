#pragma once

#include "geometrycentral/surface/geometry.h"
#include "geometrycentral/utilities/utilities.h"

#include <cmath>
#include <utility>
#include <vector>


namespace geometrycentral {
namespace surface {

VertexData<double> FMMDistance(Geometry<Euclidean>* geometry,
                               const std::vector<std::pair<Vertex, double>>& initialDistances);

VertexData<double> FMMDistance(HalfedgeMesh* mesh, const std::vector<std::pair<Vertex, double>>& initialDistances,
                               const EdgeData<double>& edgeLengths, const HalfedgeData<double>& oppAngles);

double eikonalDistanceSubroutine(double a, double b, double theta, double dA, double dB);

} // namespace surface
} // namespace geometrycentral
