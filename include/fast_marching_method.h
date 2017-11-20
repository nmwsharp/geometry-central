#pragma once

#include <vector>
#include <cmath>
#include <utility>

#include "utilities.h"
#include "geometry.h"

// TODO: Split obtuse triangles instead of being wrong.
namespace GC {

VertexData<double> FMMDistance(
    Geometry<Euclidean>* geometry,
    const std::vector<std::pair<VertexPtr, double>>& initialDistances);

VertexData<double> FMMDistance(
    HalfedgeMesh* mesh,
    const std::vector<std::pair<VertexPtr, double>>& initialDistances,
    const EdgeData<double>& edgeLengths, const HalfedgeData<double>& oppAngles);

double eikonalDistanceSubroutine(double a, double b, double theta, double dA,
                                 double dB);
}