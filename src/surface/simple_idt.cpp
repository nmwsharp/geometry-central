#include "geometrycentral/surface/simple_idt.h"

#include "geometrycentral/utilities/elementary_geometry.h"

#include <deque>

namespace geometrycentral {
namespace surface {


void flipToDelaunay(HalfedgeMesh& mesh, EdgeData<double>& edgeLengths, double delaunayEPS) {

  // TODO all of these helpers are duplicated from signpost_intrinsic_triangulation

  auto layoutDiamond = [&](Halfedge iHe) {
    // Gather index values
    Halfedge iHeA0 = iHe;
    Halfedge iHeA1 = iHeA0.next();
    Halfedge iHeA2 = iHeA1.next();
    Halfedge iHeB0 = iHe.twin();
    Halfedge iHeB1 = iHeB0.next();
    Halfedge iHeB2 = iHeB1.next();

    // Gather length values
    double l01 = edgeLengths[iHeA1.edge()];
    double l12 = edgeLengths[iHeA2.edge()];
    double l23 = edgeLengths[iHeB1.edge()];
    double l30 = edgeLengths[iHeB2.edge()];
    double l02 = edgeLengths[iHeA0.edge()];

    // Lay out the vertices of the diamond
    Vector2 p3{0., 0.};
    Vector2 p0{l30, 0.};
    Vector2 p2 = layoutTriangleVertex(p3, p0, l02, l23); // involves more arithmetic than strictly necessary
    Vector2 p1 = layoutTriangleVertex(p2, p0, l01, l12);

    return std::array<Vector2, 4>{p0, p1, p2, p3};
  };

  auto area = [&](Face f) {
    Halfedge he = f.halfedge();
    double a = edgeLengths[he.edge()];
    he = he.next();
    double b = edgeLengths[he.edge()];
    he = he.next();
    double c = edgeLengths[he.edge()];
    return triangleArea(a, b, c);
  };

  auto halfedgeCotanWeight = [&](Halfedge heI) {
    if (heI.isInterior()) {
      Halfedge he = heI;
      double l_ij = edgeLengths[he.edge()];
      he = he.next();
      double l_jk = edgeLengths[he.edge()];
      he = he.next();
      double l_ki = edgeLengths[he.edge()];
      he = he.next();
      GC_SAFETY_ASSERT(he == heI, "faces mush be triangular");
      double areaV = area(he.face());
      double cotValue = (-l_ij * l_ij + l_jk * l_jk + l_ki * l_ki) / (4. * areaV);
      return cotValue / 2;
    } else {
      return 0.;
    }
  };

  auto edgeCotanWeight = [&](Edge e) {
    return halfedgeCotanWeight(e.halfedge()) + halfedgeCotanWeight(e.halfedge().twin());
  };


  auto flipEdgeIfNotDelaunay = [&](Edge e) {
    // Can't flip
    if (e.isBoundary()) return false;

    // Don't want to flip
    double cWeight = edgeCotanWeight(e);
    if (cWeight > -delaunayEPS) return false;

    // Get geometric data
    Halfedge he = e.halfedge();
    std::array<Vector2, 4> layoutPositions = layoutDiamond(he);

    // Combinatorial flip
    bool flipped = mesh.flip(e);

    // Should always be possible, something unusual is going on if we end up here
    if (!flipped) {
      return false;
    }

    // Compute the new edge length
    double newLength = (layoutPositions[1] - layoutPositions[3]).norm();

    // If we're going to create a non-finite edge length, abort the flip
    // (only happens if you're in a bad numerical place)
    if (!std::isfinite(newLength)) {
      mesh.flip(e);
      return false;
    }

    // Assign the new edge lengths
    edgeLengths[e] = newLength;

    return true;
  };

  std::deque<Edge> edgesToCheck;
  EdgeData<char> inQueue(mesh, true);
  for (Edge e : mesh.edges()) {
    edgesToCheck.push_back(e);
  }

  size_t nFlips = 0;
  while (!edgesToCheck.empty()) {

    // Get the top element from the queue of possibily non-Delaunay edges
    Edge e = edgesToCheck.front();
    edgesToCheck.pop_front();
    inQueue[e] = false;

    bool wasFlipped = flipEdgeIfNotDelaunay(e);

    if (!wasFlipped) continue;

    // Handle the aftermath of a flip
    nFlips++;

    // Add neighbors to queue, as they may need flipping now
    Halfedge he = e.halfedge();
    Halfedge heN = he.next();
    Halfedge heT = he.twin();
    Halfedge heTN = heT.next();
    std::vector<Edge> neighEdges = {heN.edge(), heN.next().edge(), heTN.edge(), heTN.next().edge()};
    for (Edge nE : neighEdges) {
      if (!inQueue[nE]) {
        edgesToCheck.push_back(nE);
        inQueue[nE] = true;
      }
    }
  }

  std::cout << "nFlips = " << nFlips << std::endl;
}

} // namespace surface
} // namespace geometrycentral
