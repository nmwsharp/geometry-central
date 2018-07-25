#include "geometrycentral/intrinsic_geometry.h"
#include <fstream>
#include <limits>

namespace geometrycentral {

IntrinsicGeometry::IntrinsicGeometry(HalfedgeMesh& mesh_, EdgeData<double>& edgeLengths_)
    : mesh(mesh_), edgeLengths(edgeLengths_)

{}

} // namespace geometrycentral
