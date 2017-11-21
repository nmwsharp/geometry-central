#pragma once

#include <cmath>
#include <utility>
#include <vector>

#include "geometry.h"
#include "utilities.h"

// TODO: Split obtuse triangles instead of being wrong.
namespace geometrycentral {

VertexData<double> FMMDistance(
    Geometry<Euclidean>* geometry,
    const std::vector<std::pair<VertexPtr, double>>& initialDistances);

VertexData<double> FMMDistance(
    HalfedgeMesh* mesh,
    const std::vector<std::pair<VertexPtr, double>>& initialDistances,
    const EdgeData<double>& edgeLengths, const HalfedgeData<double>& oppAngles);

double eikonalDistanceSubroutine(double a, double b, double theta, double dA,
                                 double dB);

}  // namespace geometrycentral