#include "geometrycentral/surface/normal_coordinates.h"

#include <memory>

namespace geometrycentral {
namespace surface {

// helper functions
namespace {
template <typename T>
inline T positivePart(const T& t) {
  return fmax(t, 0.);
}

template <typename T>
inline T negativePart(const T& t) {
  return fmin(t, 0.);
}
} // namespace


//==============================================================================
//                          Constructors
//==============================================================================
NormalCoordinates::NormalCoordinates(ManifoldSurfaceMesh& mesh_) : mesh(mesh_) {
  edgeCoords = EdgeData<int>(mesh, 0);
  roundabouts = HalfedgeData<int>(mesh, 0);
  roundaboutDegrees = VertexData<int>(mesh, 0);
}

void NormalCoordinates::setCurvesFromEdges(ManifoldSurfaceMesh& mesh) {

  for (Edge e : mesh.edges()) {
    edgeCoords[e] = -1;
  }

  for (Vertex v : mesh.vertices()) {
    size_t iHe = 0;
    size_t D = v.degree();
    roundaboutDegrees[v] = D;

    // Explicitly loop over halfedges in counterclockwise order
    Halfedge he = v.halfedge();
    do {
      roundabouts[he] = iHe;

      // If we have reached a boundary, we're done
      // (Assumes that for boundary vertices, vertex.halfedge() is the one
      // on the interior of the mesh pointing along the boundary)
      // TODO: verify this
      if (!he.isInterior()) break;

      iHe = (iHe + 1) % D;
      he = he.next().next().twin();
    } while (he != v.halfedge());
  }
}

//==============================================================================
//                          Accessors
//==============================================================================

int NormalCoordinates::operator[](Edge e) const { return edgeCoords[e]; }
int NormalCoordinates::operator[](Corner c) const { return cornerCoord(c); }
int& NormalCoordinates::operator[](Edge e) { return edgeCoords[e]; }

//==============================================================================
//                          General Routines
//==============================================================================
// Call immediately before mesh->flip(e);
std::tuple<int, size_t, size_t> NormalCoordinates::computeFlippedData(Edge e) {
  /*       Flip rotates edge clockwise
   *          k                   k
   *        / |\                / ||\
   *       /    \              /  |  \
   *     eki     ejk          /   |   \
   *     /        \         e2    |    e1
   *   |/          \       |/     |     \
   *  i --- eij ---> j    i \     v     /| j
   *    \          /|        \    |    /
   *     \        /           \   |   /
   *    eil      elj          e3  |  e4
   *       \    /               \ | /
   *        \| /                _\|/
   *         l                    l
   */
  Halfedge he = e.halfedge();
  Edge ejk = he.next().edge();
  Edge eki = he.next().next().edge();
  Edge eil = he.twin().next().edge();
  Edge elj = he.twin().next().next().edge();

  // Gather intersection numbers
  int nij = edgeCoords[e];
  int njk = edgeCoords[ejk];
  int nki = edgeCoords[eki];
  int nil = edgeCoords[eil];
  int nlj = edgeCoords[elj];

  int nkl = flipNormalCoordinates(nij, njk, nki, nil, nlj);

  // Compute flipped roundabouts
  size_t rki = roundabouts[he.next().next()];
  size_t rlj = roundabouts[he.twin().next().next()];
  size_t dk = roundaboutDegrees[he.next().next().vertex()];
  size_t dl = roundaboutDegrees[he.twin().next().next().vertex()];

  size_t rkl, rlk;
  std::tie(rkl, rlk) = flipRoundabouts(nij, njk, nki, nil, nlj, nkl, rki, rlj, dk, dl);

  //== After flip, set:
  // edgeCoords[e]          = nkl;
  // roundabouts[he]        = rkl;
  // roundabouts[he.twin()] = rlk;

  return std::make_tuple(nkl, rkl, rlk);
}

void NormalCoordinates::applyFlippedData(Edge e, const std::tuple<int, size_t, size_t>& update) {
  // update = nkl, rkl, rlk

  Halfedge he = e.halfedge();
  edgeCoords[e] = std::get<0>(update);
  roundabouts[he] = std::get<1>(update);
  roundabouts[he.twin()] = std::get<2>(update);
}


// Call immediately before mesh->insertVertex(f);
// Input are the desired normal coordinates of the edges connecting face f's
// vertices to the new vertex, in the order given by f.adjacentVertices()
// TODO: validate that these inputs are possible
std::array<int, 3> NormalCoordinates::computeVertexInsertionData(Face f, const std::array<int, 3>& newCrossingCounts) {
  return {newCrossingCounts[0], newCrossingCounts[2], newCrossingCounts[1]};
}

std::array<int, 3> NormalCoordinates::computeVertexInsertionDataGeodesic(IntrinsicGeometryInterface& geo, Face f,
                                                                         Vector3 location) {

  // Populate the crossing locations for the edges of the triangle
  size_t iHe = 0;
  std::array<std::vector<double>, 3> boundaryCrossings;
  for (Halfedge he : f.adjacentHalfedges()) {
    boundaryCrossings[iHe] = generateGeodesicCrossingLocations(geo, he);
    iHe++;
  }

  // TODO apply some sanity policies, like that the crossings should be
  // correctly ordered

  std::array<int, 3> counts = computeVertexInsertionCrossingCounts(location, boundaryCrossings);

  return computeVertexInsertionData(f, counts);
}

// Call after mesh->insertVertex(f) to update normal coordinates;
void NormalCoordinates::applyVertexInsertionData(Vertex newVertex, const std::array<int, 3>& update) {
  size_t iE = 0;
  GC_SAFETY_ASSERT(newVertex.degree() == 3, "Inserting a vertex into a triangle should produce a degree "
                                            "3 vertex, but " +
                                                std::to_string(newVertex) + " has degree " +
                                                std::to_string(newVertex.degree()));
  for (Edge e : newVertex.adjacentEdges()) {
    edgeCoords[e] = update[iE++];
  }

  for (Corner c : newVertex.adjacentCorners()) {
    GC_SAFETY_ASSERT(strictDegree(c) == 0, "inserted vertices cannot touch curves");
  }
}

std::array<int, 4> NormalCoordinates::computeInteriorEdgeSplitDataGeodesic(IntrinsicGeometryInterface& geo, Edge e,
                                                                           double location) {

  /*       Flip rotates edge clockwise
   *          k
   *        / |\
   *       /    \
   *     eki     ejk
   *     /        \
   *   |/          \
   *  i --- eij ---> j
   *    \          /|
   *     \        /
   *    eil      elj
   *       \    /
   *        \| /
   *         l
   */

  if (edgeCoords[e] > 0) {
    std::vector<double> eCrossings = generateGeodesicCrossingLocations(geo, e.halfedge());

    int n2 = edgeCoords[e];
    int n4 = 0;

    for (double tCross : eCrossings) {
      if (location > tCross) {
        n2--;
        n4++;
      }
    }

    Corner ck = e.halfedge().next().next().corner();
    Corner cl = e.halfedge().twin().next().next().corner();

    int n1 = 0;
    n1 += positivePart(strictCornerCoord(e.halfedge().twin().corner()) - positivePart(n2));
    n1 += positivePart(strictCornerCoord(e.halfedge().twin().next().corner()) - positivePart(n4));
    n1 += strictCornerCoord(cl);
    n1 += strictDegree(e.halfedge().twin().corner());
    n1 += strictDegree(e.halfedge().twin().next().corner());

    int n3 = 0;
    n3 += positivePart(strictCornerCoord(e.halfedge().next().corner()) - positivePart(n2));
    n3 += positivePart(strictCornerCoord(e.halfedge().corner()) - positivePart(n4));
    n3 += strictCornerCoord(ck);
    n3 += strictDegree(e.halfedge().corner());
    n3 += strictDegree(e.halfedge().next().corner());

    return {n1, n2, n3, n4};
  } else {
    int nil = edgeCoords[e.halfedge().twin().next().edge()];
    int nlj = edgeCoords[e.halfedge().twin().next().next().edge()];
    int njk = edgeCoords[e.halfedge().next().edge()];
    int nki = edgeCoords[e.halfedge().next().next().edge()];

    int n1 = fmax(nil, fmax(nlj, 0));
    int n2 = edgeCoords[e];
    int n3 = fmax(njk, fmax(nki, 0));
    int n4 = edgeCoords[e];
    return {n1, n2, n3, n4};
  }
}

std::array<int, 4> NormalCoordinates::computeInteriorEdgeSplitDataCombinatorial(IntrinsicGeometryInterface& geo, Edge e,
                                                                                size_t iSeg) {
  if (edgeCoords[e] > 0) {
    int n2 = edgeCoords[e] - iSeg;
    int n4 = iSeg;

    Corner ck = e.halfedge().next().next().corner();
    Corner cl = e.halfedge().twin().next().next().corner();

    int n1 = 0;
    n1 += positivePart(strictCornerCoord(e.halfedge().twin().corner()) - positivePart(n2));
    n1 += positivePart(strictCornerCoord(e.halfedge().twin().next().corner()) - positivePart(n4));
    n1 += strictCornerCoord(cl);
    n1 += strictDegree(e.halfedge().twin().corner());
    n1 += strictDegree(e.halfedge().twin().next().corner());

    int n3 = 0;
    n3 += positivePart(strictCornerCoord(e.halfedge().next().corner()) - positivePart(n2));
    n3 += positivePart(strictCornerCoord(e.halfedge().corner()) - positivePart(n4));
    n3 += strictCornerCoord(ck);
    n3 += strictDegree(e.halfedge().corner());
    n3 += strictDegree(e.halfedge().next().corner());

    return {n1, n2, n3, n4};
  } else {
    int nil = edgeCoords[e.halfedge().twin().next().edge()];
    int nlj = edgeCoords[e.halfedge().twin().next().next().edge()];
    int njk = edgeCoords[e.halfedge().next().edge()];
    int nki = edgeCoords[e.halfedge().next().next().edge()];

    int n1 = fmax(nil, fmax(nlj, 0));
    int n2 = edgeCoords[e];
    int n3 = fmax(njk, fmax(nki, 0));
    int n4 = edgeCoords[e];
    return {n1, n2, n3, n4};
  }
}

std::array<int, 3> NormalCoordinates::computeBoundaryEdgeSplitDataGeodesic(IntrinsicGeometryInterface& geo, Edge e,
                                                                           double location) {

  GC_SAFETY_ASSERT(e.isBoundary(), "can only call computeBoundaryEdgeSplitDataGeodesic on a "
                                   "boundary edge");

  Halfedge heI = e.halfedge();
  if (!heI.isInterior()) {
    heI = heI.twin();
    location = 1 - location;
  }
  GC_SAFETY_ASSERT(heI.isInterior(), "One of the halfedges must be interior");

  /*       Flip rotates edge clockwise
   *          k
   *        / |\
   *       /    \
   *     eki     ejk
   *     /        \
   *   |/          \
   *  i --- eij ---> j
   */

  if (edgeCoords[e] > 0) {
    std::vector<double> eCrossings = generateGeodesicCrossingLocations(geo, heI);

    int n2 = edgeCoords[e];
    int n4 = 0;

    for (double tCross : eCrossings) {
      if (location > tCross) {
        n2--;
        n4++;
      }
    }

    Corner ck = heI.next().next().corner();

    int n3 = 0;
    n3 += positivePart(strictCornerCoord(heI.next().corner()) - positivePart(n2));
    n3 += positivePart(strictCornerCoord(heI.corner()) - positivePart(n4));
    n3 += strictCornerCoord(ck);
    n3 += strictDegree(heI.corner());
    n3 += strictDegree(heI.next().corner());

    return {n2, n3, n4};
  } else {
    int njk = edgeCoords[e.halfedge().next().edge()];
    int nki = edgeCoords[e.halfedge().next().next().edge()];

    int n2 = edgeCoords[e];
    int n3 = fmax(njk, fmax(nki, 0));
    int n4 = edgeCoords[e];
    return {n2, n3, n4};
  }
}


// As usual, index is be 0-indexed
NormalCoordinatesCurve NormalCoordinates::topologicalTrace(Halfedge he, int index) const {
  Halfedge startHe = he;
  int startIndex = index;

  NormalCoordinatesCurve path;

  do {
    int n = edgeCoords[he.edge()];

    // GC_SAFETY_ASSERT(0 <= index && index < n,
    //"invalid tracing index : " + std::to_string(index) +
    //" of " + std::to_string(n));
    // remove >= 0 case which always holds
    GC_SAFETY_ASSERT(index < n, "invalid tracing index : " + std::to_string(index) + " of " + std::to_string(n));

    path.crossings.push_back(std::make_pair(index, he));

  } while (!stepTopologicalCurve(he, index) && !(he == startHe && index == startIndex) && !he.edge().isBoundary());

  // If we made a complete loop, record the starting point (also the ending
  // point) again
  if ((he == startHe && index == startIndex) || he.edge().isBoundary())
    path.crossings.push_back(std::make_pair(index, he));

  return path;
}

// As usual, index is be 0-indexed
NormalCoordinatesCurve NormalCoordinates::topologicalTrace(Corner c, int index) const {
  Halfedge he = c.halfedge().next();
  GC_SAFETY_ASSERT(strictDegree(c) > 0, "Tried to trace paths out of a corner that no paths come out of");

  // Take the positive part of the adjacent normal coordinate. When we're
  // doing this tracing, we don't care about arcs parallel to c.halfedge()
  int offset = positivePart(edgeCoords[c.halfedge().edge()]);

  // Return early if curve immediately hits the boundary
  if (he.edge().isBoundary()) {
    return NormalCoordinatesCurve{{std::make_pair(index + offset, he)}};
  }

  return topologicalTrace(he, index + offset);
}


// Return curve that intersects he w/ positive orientation (i.e. he is in the
// list of halfedges crossed by the curve, rather than he.twin())
// Recall also that curves always intersect left-pointing halfedges
std::tuple<NormalCoordinatesCurve, int> NormalCoordinates::topologicalTraceBidirectional(Halfedge he, int index) const {

  // doesn't make sense to call for coincident edges (or if there are no
  // crossings)
  GC_SAFETY_ASSERT(edgeCoords[he.edge()] > 0, "should not be coincident, or have no crossings");


  // == Trace forward
  NormalCoordinatesCurve forwardTrace = topologicalTrace(he, index);

  // Case where the curve is a cycle
  // TODO this may be wrong
  if (forwardTrace.crossings.size() > 1 &&
      forwardTrace.crossings.front().second == forwardTrace.crossings.back().second) {
    throw std::runtime_error("encountered a loop while tracing normal coordinates");
  }

  // == Trace backward
  NormalCoordinatesCurve backwardTrace = topologicalTrace(he.twin(), edgeCoords[he.edge()] - index - 1);

  // == Combine

  // Reverse and flip the backward trace list
  std::reverse(backwardTrace.crossings.begin(), backwardTrace.crossings.end());
  for (std::pair<int, Halfedge>& p : backwardTrace.crossings) {
    p.first = edgeCoords[p.second.edge()] - p.first - 1;
    p.second = p.second.twin();
  }

  // Concatenate forward list after second
  int origIndex = backwardTrace.crossings.size() - 1;
  for (size_t i = 1; i < forwardTrace.crossings.size(); i++) {
    backwardTrace.crossings.push_back(forwardTrace.crossings[i]);
  }

  return std::tuple<NormalCoordinatesCurve, int>{backwardTrace, origIndex};
}

// HACK: represents arcs parallel to a mesh edge with a single pair {-n, he}
// where n is the number of arcs parallel to he.edge()
// TODO: Doesn't handle boundary properly
std::vector<NormalCoordinatesCurve> NormalCoordinates::topologicalTraceAllCurves() const {
  std::vector<NormalCoordinatesCurve> curves;

  //== Trace paths
  for (Corner c : mesh.corners()) {
    for (size_t iL = 0; iL < strictDegree(c); ++iL) {
      curves.push_back(topologicalTrace(c, iL));
      // cout << curves.back() << myendl;
    }
  }

  //== Trace remaining closed loops
  // Record whether we've extraced a curve going through each intersection
  EdgeData<std::vector<char>> visited(mesh);
  for (Edge e : mesh.edges()) visited[e] = std::vector<char>(positivePart(edgeCoords[e]), false);

  auto visitCurve = [&](const NormalCoordinatesCurve& curve) {
    for (const auto& p : curve.crossings) {
      int iC = std::get<0>(p);
      Halfedge he = std::get<1>(p);
      Edge e = he.edge();
      GC_SAFETY_ASSERT(0 <= iC && iC < edgeCoords[e], "traced invalid intersection");

      // Flip index if halfedge points backwards
      if (he != e.halfedge()) {
        iC = edgeCoords[e] - iC - 1;
      }

      visited[e][iC] = true;
    }
  };

  for (const auto& curve : curves) {
    visitCurve(curve);
  }

  // Trace out remaining loops. Start with boundary edges
  for (BoundaryLoop b : mesh.boundaryLoops()) {
    for (Edge e : b.adjacentEdges()) {
      for (size_t iC = 0; (int)iC < edgeCoords[e]; ++iC) {
        if (!visited[e][iC]) {
          auto curve = topologicalTrace(e.halfedge().twin(), edgeCoords[e] - iC - 1);
          visitCurve(curve);
          curves.push_back(curve);
        }
      }
    }
  }
  for (Edge e : mesh.edges()) {
    for (size_t iC = 0; (int)iC < edgeCoords[e]; ++iC) {
      if (!visited[e][iC]) {
        auto curve = topologicalTrace(e.halfedge(), iC);
        visitCurve(curve);
        curves.push_back(curve);
      }
    }
  }

  //==Add in paths parallel to mesh edges
  for (Edge e : mesh.edges()) {
    if (edgeCoords[e] < 0) {
      curves.push_back(NormalCoordinatesCurve{{std::make_pair(edgeCoords[e], e.halfedge())}});
    }
  }

  return curves;
}

std::vector<std::vector<SurfacePoint>> NormalCoordinates::generateAnyGeometry() const {
  std::vector<NormalCoordinatesCurve> curves = topologicalTraceAllCurves();

  std::vector<std::vector<SurfacePoint>> embeddedCurves;
  for (const NormalCoordinatesCurve& curve : curves) {
    if (std::get<0>(curve.crossings[0]) >= 0) {
      std::vector<SurfacePoint> embeddedCurve;

      // NormalCoordinatesCurve of size 1 (just immediately goes from a corner into the
      // boundary)
      if (curve.crossings.size() == 1) {

        Halfedge heStart = std::get<1>(curve.crossings[0]);
        SurfacePoint vStart(heStart.next().next().vertex());
        embeddedCurve.push_back(vStart);

        int index = std::get<0>(curve.crossings[0]);
        Halfedge he = std::get<1>(curve.crossings[0]);
        GC_SAFETY_ASSERT(index >= 0, "Only paths along a single mesh edge are "
                                     "allowed to have negative indices");
        GC_SAFETY_ASSERT(edgeCoords[he.edge()] > 0, "Can't cross an edge with nonpositive edgeCoord");

        double t = (index + 1) / (edgeCoords[he.edge()] + 1.);
        embeddedCurve.push_back(SurfacePoint(he, t));
      }

      // In this case, the path crosses mesh edges
      if (curve.crossings[0] == curve.crossings[curve.crossings.size() - 1]) {
        // In this case, the curve is a loop
        for (const std::pair<int, Halfedge>& pt : curve.crossings) {
          int index = std::get<0>(pt);
          Halfedge he = std::get<1>(pt);
          GC_SAFETY_ASSERT(index >= 0, "Only paths along a single mesh edge are "
                                       "allowed to have negative indices");
          GC_SAFETY_ASSERT(edgeCoords[he.edge()] > 0, "Can't cross an edge with nonpositive edgeCoord");

          double t = (index + 1) / (edgeCoords[he.edge()] + 1.);
          embeddedCurve.push_back(SurfacePoint(he, t));
        }

      } else {
        // In this case, the curve is a path

        // If the curve starts on a boundary edge, it doesn't get a
        // starting vertex
        Halfedge heStart = std::get<1>(curve.crossings[0]);
        if (!heStart.edge().isBoundary()) {
          SurfacePoint vStart(heStart.next().next().vertex());
          embeddedCurve.push_back(vStart);
        }

        for (const std::pair<int, Halfedge>& pt : curve.crossings) {
          int index = std::get<0>(pt);
          Halfedge he = std::get<1>(pt);
          GC_SAFETY_ASSERT(index >= 0, "Only paths along a single mesh edge are "
                                       "allowed to have negative indices");
          GC_SAFETY_ASSERT(edgeCoords[he.edge()] > 0, "Can't cross an edge with nonpositive edgeCoord");

          double t = (index + 1) / (edgeCoords[he.edge()] + 1.);
          embeddedCurve.push_back(SurfacePoint(he, t));
        }

        // If the curve ends on a boundary edge, it doesn't get an
        // ending vertex
        Halfedge heEnd = std::get<1>(curve.crossings[curve.crossings.size() - 1]);
        if (!heEnd.edge().isBoundary()) {
          SurfacePoint vEnd(heEnd.twin().next().next().vertex());
          embeddedCurve.push_back(vEnd);
        }
      }

      embeddedCurves.push_back(embeddedCurve);
    } else {
      // In this case, the path goes along a single mesh edge
      GC_SAFETY_ASSERT(curve.crossings.size() == 1, "Only paths along a single mesh edge are allowed to "
                                                    "have negative indices");

      Halfedge he = std::get<1>(curve.crossings[0]);
      SurfacePoint vStart(he.vertex());
      SurfacePoint vEnd(he.next().vertex());

      std::vector<SurfacePoint> embeddedCurve{vStart, vEnd};
      embeddedCurves.push_back(embeddedCurve);
    }
  }
  return embeddedCurves;
}

//==============================================================================
//                            Geodesic Routines
//==============================================================================
std::vector<std::vector<SurfacePoint>>
NormalCoordinates::generateGeodesicGeometry(IntrinsicGeometryInterface& geo) const {

  std::vector<NormalCoordinatesCurve> curves = topologicalTraceAllCurves();

  return geometrycentral::surface::generateGeodesicGeometry(mesh, geo, curves);
}


double NormalCoordinates::generateGeodesicCrossingLocation(IntrinsicGeometryInterface& geo, Halfedge he,
                                                           int ind) const {

  // Get the topological crossings for the curve
  NormalCoordinatesCurve crossings;
  int centerCrossInd;
  std::tie(crossings, centerCrossInd) = topologicalTraceBidirectional(he, ind);

  // Generate the geometry for those crossings by laying out
  std::vector<SurfacePoint> crossingLocations = generateSingleGeodesicGeometry(mesh, geo, crossings);

  // Get the crossing
  SurfacePoint& thisCross = crossingLocations[centerCrossInd + 1];

  GC_SAFETY_ASSERT(thisCross.type == SurfacePointType::Edge, "crossing should be an edge point");


  // Flip if wrong halfedge orientation
  if (he.edge().halfedge() == he) {
    return thisCross.tEdge;
  } else {
    return 1. - thisCross.tEdge;
  }
}

std::vector<double> NormalCoordinates::generateGeodesicCrossingLocations(IntrinsicGeometryInterface& geo,
                                                                         Halfedge he) const {
  std::vector<double> crossings;
  for (int ind = 0; ind < edgeCoords[he.edge()]; ind++) {
    double tCross = generateGeodesicCrossingLocation(geo, he, ind);
    crossings.push_back(tCross);
  }
  return crossings;
}

//==============================================================================
//                                Helpers
//==============================================================================

// Compute a corner's coordinate from the edge coordinates in its triangle
int NormalCoordinates::cornerCoord(Corner c) const {
  int nki = edgeCoords[c.halfedge().edge()];
  int nij = edgeCoords[c.halfedge().next().edge()];
  int njk = edgeCoords[c.halfedge().next().next().edge()];

  return geometrycentral::surface::cornerCoord(nij, njk, nki);
}

// Counts how many arcs come out of this corner and exit at the opposite
// edge (i.e. ignores arcs parallel to triangle edges)
size_t NormalCoordinates::strictDegree(Corner c) const {
  int nki = edgeCoords[c.halfedge().edge()];
  int nij = edgeCoords[c.halfedge().next().edge()];
  int njk = edgeCoords[c.halfedge().next().next().edge()];

  return geometrycentral::surface::strictDegree(nij, njk, nki);
}

// Counts how many arcs clip off this corner, ignoring arcs that emanate
size_t NormalCoordinates::strictCornerCoord(Corner c) const {
  int nki = edgeCoords[c.halfedge().edge()];
  int nij = edgeCoords[c.halfedge().next().edge()];
  int njk = edgeCoords[c.halfedge().next().next().edge()];

  return geometrycentral::surface::strictCornerCoord(nij, njk, nki);
}

// He is the twin of the halfedge in a face that the curve starts from,
// index is the index of the curve, with 0 being closest to the source
// vertex of he.
// If the curve terminates at a vertex, returns true
// Otherwise returns false and udpates he and index to be the next crossing
// As usual, index is 0-indexed
bool NormalCoordinates::stepTopologicalCurve(Halfedge& he, int& index) const {
  he = he.twin();
  int cikj = strictCornerCoord(he.next().corner());
  int cjik = strictCornerCoord(he.corner());
  int nij = positivePart(edgeCoords[he.edge()]);

  if (index < cikj) {
    he = he.next();
    // index = index;
    return false;
  } else if (index >= nij - cjik) {
    int nkj = positivePart(edgeCoords[he.next().next().edge()]);
    he = he.next().next();
    index = nkj + index - nij;
    return false;
  } else {
    return true;
  }
}


bool NormalCoordinates::triangleInequalityViolation(Face f, Halfedge& violatingHe) const {
  Halfedge ha = f.halfedge();
  Halfedge hb = ha.next();
  Halfedge hc = hb.next();

  GC_SAFETY_ASSERT(ha == hc.next(), "Attempting to use normal coordinates on a non-triangular face");

  size_t a = positivePart(edgeCoords[ha.edge()]);
  size_t b = positivePart(edgeCoords[hb.edge()]);
  size_t c = positivePart(edgeCoords[hc.edge()]);

  if (a > b + c) {
    violatingHe = ha;
    return true;
  } else if (b > c + a) {
    violatingHe = hb;
    return true;
  } else if (c > a + b) {
    violatingHe = hc;
    return true;
  } else {
    return false;
  }
}


bool NormalCoordinates::isEncircledByLoopCurve(Vertex v) const {
  for (Corner c : v.adjacentCorners()) {
    if (cornerCoord(c) <= 0) return false;
  }
  return true;
}

bool NormalCoordinates::isHookedByCurve(Vertex v) const {
  bool seenEnterExit = false;
  for (Corner c : v.adjacentCorners()) {
    if (cornerCoord(c) <= 0) {
      if (seenEnterExit) return false;
      seenEnterExit = true;
    }
  }
  return seenEnterExit;
}

void NormalCoordinates::setRoundaboutFromPrevRoundabout(Halfedge he) {
  if (!he.isInterior()) {
    // Any non-interior halfedge must be the the last one, and must be
    // shared if we're storing a triangulation
    roundabouts[he] = roundaboutDegrees[he.vertex()] - 1;
  } else {
    Halfedge hePrev = he.twin().next();
    size_t deltaPrev = -negativePart(edgeCoords[hePrev.edge()]);

    size_t rPrev = roundabouts[hePrev];
    size_t em = strictDegree(hePrev.corner());
    size_t dV = roundaboutDegrees[he.vertex()];

    size_t rNew = dV == 0 ? 0 : (rPrev + deltaPrev + em) % dV;

    roundabouts[he] = rNew;
  }
}

/*       Flip rotates edge clockwise
 *          k                   k
 *        / |\                / ||\
 *       /    \              /  |  \
 *     eki     ejk          /   |   \
 *     /        \         e2    |    e1
 *   |/          \       |/     |     \
 *  i --- eij ---> j    i \     v     /| j
 *    \          /|        \    |    /
 *     \        /           \   |   /
 *    eil      elj          e3  |  e4
 *       \    /               \ | /
 *        \| /                _\|/
 *         l                    l
 */
int flipNormalCoordinates(int nij, int njk, int nki, int nil, int nlj) {
  // For any edge except nij, we can forget about the negative part, since
  // arcs along those edges cannot intersect edge kl
  // Once we do this, we can use the flip formula from CEPS
  // TODO: can we simplify the formula now that we allow negative n?

  njk = positivePart(njk);
  nki = positivePart(nki);
  nil = positivePart(nil);
  nlj = positivePart(nlj);


  // For nij, we can decompose it into its positive and negative parts (only
  // one of which can be nonzero)
  int arcsAlongIJ = -negativePart(nij);
  nij = positivePart(nij);

  // Compute flipped edgeCoord
  int Eilj = positivePart(nlj - nij - nil);
  int Ejil = positivePart(nil - nlj - nij);
  int Elji = positivePart(nij - nil - nlj);

  int Eijk = positivePart(njk - nki - nij);
  int Ejki = positivePart(nki - nij - njk);
  int Ekij = positivePart(nij - njk - nki);

  double Cilj = -(negativePart(nlj - nij - nil) + Ejil + Elji) / 2.;
  double Cjil = -(negativePart(nil - nlj - nij) + Eilj + Elji) / 2.;
  double Clji = -(negativePart(nij - nil - nlj) + Eilj + Ejil) / 2.;

  double Cijk = -(negativePart(njk - nki - nij) + Ejki + Ekij) / 2.;
  double Cjki = -(negativePart(nki - nij - njk) + Eijk + Ekij) / 2.;
  double Ckij = -(negativePart(nij - njk - nki) + Eijk + Ejki) / 2.;

  int doubledAnswer = 2 * Clji + 2 * Ckij + abs(Cjil - Cjki) + abs(Cilj - Cijk) - Elji - Ekij + 2 * Eilj + 2 * Eijk +
                      2 * Ejil + 2 * Ejki;

  int answer = std::round(doubledAnswer / 2) + arcsAlongIJ;
  return answer;
}


// See NormalCoordinates::computeFlippedData for an explanation of all these
// variable names;
std::pair<size_t, size_t> flipRoundabouts(int nij, int njk, int nki, int nil, int nlj, int nkl, size_t rki, size_t rlj,
                                          size_t dk, size_t dl) {
  size_t deltaki = -negativePart(nki);
  size_t deltalj = -negativePart(nlj);

  size_t ekil = strictDegree(nil, nki, nkl);
  size_t eljk = strictDegree(njk, nlj, nkl);

  size_t rkl = dk == 0 ? 0 : (rki + deltaki + ekil) % dk;
  size_t rlk = dl == 0 ? 0 : (rlj + deltalj + eljk) % dl;

  return std::make_pair(rkl, rlk);
}

// Count the number of arcs strictly leaving corner k (i.e. not parallel to jk
// or ki)
size_t strictDegree(int nij, int njk, int nki) {
  return positivePart(positivePart(nij) - positivePart(njk) - positivePart(nki));
}

// Counts the number of arcs clipping off corner k. If arcs leave corner k, then
// they get counted as -1/2 each
int cornerCoord(int nij, int njk, int nki) {
  return (njk + nki - nij - strictDegree(njk, nki, nij) - strictDegree(nki, nij, njk)) / 2;
}

// Counts the number of arcs clipping off corner k, returns 0 if arcs leave
// corner k
size_t strictCornerCoord(int nij, int njk, int nki) { return positivePart(cornerCoord(nij, njk, nki)); }

// Given trace info (like the output from topotrace, generate geodesic
// geometry)
std::vector<std::vector<SurfacePoint>>
generateGeodesicGeometry(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geo,
                         const std::vector<NormalCoordinatesCurve>& traceCounts) {

  std::vector<std::vector<SurfacePoint>> embeddedCurves;
  for (const NormalCoordinatesCurve& curve : traceCounts) {
    if (std::get<0>(curve.crossings[0]) >= 0) {
      // In this case, the path crosses mesh edges
      embeddedCurves.push_back(generateSingleGeodesicGeometry(mesh, geo, curve));
    } else {
      // In this case, the path goes along a single mesh edge
      GC_SAFETY_ASSERT(curve.crossings.size() == 1, "Only paths along a single mesh edge are allowed to "
                                                    "have negative indices");

      Halfedge he = std::get<1>(curve.crossings[0]);
      SurfacePoint vStart(he.vertex());
      SurfacePoint vEnd(he.next().vertex());

      std::vector<SurfacePoint> embeddedCurve{vStart, vEnd};
      embeddedCurves.push_back(embeddedCurve);
    }
  }
  return embeddedCurves;
}

std::vector<SurfacePoint> generateSingleGeodesicGeometry(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geo,
                                                         const NormalCoordinatesCurve& curve) {

  std::vector<SurfacePoint> toRet;
  std::vector<std::pair<SurfacePoint, double>> fullGeo = generateFullSingleGeodesicGeometry(mesh, geo, curve);
  for (const std::pair<SurfacePoint, double>& pt : fullGeo) {
    toRet.push_back(std::get<0>(pt));
  }
  return toRet;
}

// Record barycentric coordinate of each SurfacePoint along the curve
// Barycentric coordinates are 0 at the src of and edge and 1 at the dst
std::vector<std::pair<SurfacePoint, double>> generateFullSingleGeodesicGeometry(ManifoldSurfaceMesh& mesh,
                                                                                IntrinsicGeometryInterface& geo,
                                                                                const NormalCoordinatesCurve& curve,
                                                                                bool verbose) {

  // == Compute a geodesic path by laying the triangle strip that the path
  // passes through out in the plane
  // Do the layout by first computing the orientation of each halfedge in the
  // triangle strip, and then laying out faces using those directions
  // Can I get away without requiring these?
  // geo.requireEdgeLengths();

  // Stolen from edge_length_geometry.ipp
  // TODO: refactor better
  auto cornerAngle = [&](Corner c) -> double {
    Halfedge heA = c.halfedge();
    Halfedge heOpp = heA.next();
    Halfedge heB = heOpp.next();

    GC_SAFETY_ASSERT(heB.next() == heA, "faces must be triangular");

    double lOpp = geo.edgeLengths[heOpp.edge()];
    double lA = geo.edgeLengths[heA.edge()];
    double lB = geo.edgeLengths[heB.edge()];

    double q = (lA * lA + lB * lB - lOpp * lOpp) / (2. * lA * lB);
    q = clamp(q, -1.0, 1.0);
    double angle = std::acos(q);

    return angle;
  };

  // Helper functions to lay out faces from different sets of input data

  auto layOutFirstFace = [&](Halfedge he) {
    double len = geo.edgeLengths[he.edge()];
    Vector2 pos1 = Vector2{0, 0};
    Vector2 pos2 = Vector2{len, 0};
    len = geo.edgeLengths[he.next().next().edge()];
    double angle = cornerAngle(he.corner());
    Vector2 pos0 = Vector2{len * cos(angle), len * sin(angle)};

    return std::make_tuple(pos0, pos1, pos2);
  };

  auto nextAngle = [&](Halfedge he, Halfedge prevHe, double prevAngle) {
    double angle = cornerAngle(prevHe.twin().next().corner());
    double theta = prevAngle - angle;

    if (prevHe.twin().next() == he) {
      return std::make_tuple(theta, theta);
    } else if (prevHe.twin().next().next() == he) {
      double nextAngle = cornerAngle(prevHe.twin().next().next().corner());
      return std::make_tuple(theta, theta + M_PI - nextAngle);
    } else {
      throw std::runtime_error("Path connectivity error in "
                               "Tracing.cpp:straightenEuclidean:nextAngle");
    }
  };

  auto layOutNextFaceFromAngle = [&](Halfedge he, Halfedge prevHe, const std::array<Vector2, 3>& prevPositions,
                                     double theta) {
    double len = geo.edgeLengths[prevHe.twin().next().edge()];

    Vector2 newPos = prevPositions[1] + Vector2{len * cos(theta), len * sin(theta)};

    if (prevHe.twin().next() == he) {
      return std::tuple<Vector2, Vector2, Vector2>{prevPositions[2], prevPositions[1], newPos};
    } else if (prevHe.twin().next().next() == he) {
      return std::tuple<Vector2, Vector2, Vector2>{prevPositions[1], newPos, prevPositions[2]};
    } else {
      throw std::runtime_error("Path connectivity error in "
                               "Tracing.cpp:straightenEuclidean:layOutNextFaceFromAngle");
    }
  };

  // TODO: replace with placeTrianglePoint or whatever it's called
  auto layOutLastFace = [&](Halfedge he, Halfedge prevHe, const std::array<Vector2, 3>& prevPositions) {
    double angle = cornerAngle(prevHe.twin().next().corner());
    double len = geo.edgeLengths[prevHe.twin().next().edge()];
    Vector2 disp = prevPositions[2] - prevPositions[1];
    double theta = atan2(disp.y, disp.x) - angle;

    Vector2 newPos = prevPositions[1] + Vector2{len * cos(theta), len * sin(theta)};


    // TODO: it's kind of wasteful to return a whole tuple, but it makes
    // my life easier for now
    if (prevHe.twin().next() == he) {
      return std::tuple<Vector2, Vector2, Vector2>{prevPositions[2], prevPositions[1], newPos};
    } else if (prevHe.twin().next().next() == he) {
      return std::tuple<Vector2, Vector2, Vector2>{prevPositions[1], newPos, prevPositions[2]};
    } else {
      exit(1);
    }
  };

  std::vector<std::array<Vector2, 2>> edgePositions;
  std::vector<double> angles;

  Vector2 base, src, dst;

  Halfedge prevHe = std::get<1>(curve.crossings[0]);

  std::tie(base, src, dst) = layOutFirstFace(prevHe);
  edgePositions.push_back({src, dst});

  Vector2 start = base;

  double heAngle = 0, refAngle = 0;
  for (size_t iC = 1; iC < curve.crossings.size(); ++iC) {
    Halfedge he = std::get<1>(curve.crossings[iC]);
    std::tie(refAngle, heAngle) = nextAngle(he, prevHe, heAngle);

    angles.push_back(refAngle);
    prevHe = he;
  }

  prevHe = std::get<1>(curve.crossings[0]);
  for (size_t iC = 1; iC < curve.crossings.size(); ++iC) {
    Halfedge he = std::get<1>(curve.crossings[iC]);
    std::tie(base, src, dst) = layOutNextFaceFromAngle(he, prevHe, {base, src, dst}, angles[iC - 1]);

    edgePositions.push_back({src, dst});
    prevHe = he;
  }

  Halfedge almostLastHe = std::get<1>(curve.crossings[curve.crossings.size() - 1]);
  // lastHe points to the ending vertex
  Halfedge lastHe = almostLastHe.twin().next();
  std::tie(base, src, dst) = layOutLastFace(lastHe, almostLastHe, {base, src, dst});
  Vector2 end = dst;

  Halfedge heStart = std::get<1>(curve.crossings[0]);
  SurfacePoint vStart(heStart.next().next().vertex());
  Halfedge heEnd = std::get<1>(curve.crossings[curve.crossings.size() - 1]);
  SurfacePoint vEnd(heEnd.twin().next().next().vertex());

  if (verbose) {
    std::cout << "start: " << start << std::endl;
    for (const auto& positions : edgePositions) {
      std::cout << "  " << positions[0] << "   |   " << positions[1] << std::endl;
    }
    std::cout << "end: " << start << std::endl;
  }

  std::vector<std::pair<SurfacePoint, double>> embeddedCurve{std::make_pair(vStart, 0)};

  for (size_t iC = 0; iC < curve.crossings.size(); ++iC) {
    Halfedge he = std::get<1>(curve.crossings[iC]);

    auto intersection = segmentSegmentIntersection(edgePositions[iC][0], edgePositions[iC][1], start, end);

    double tA = intersection.tA;
    double tB = intersection.tB;

    if (verbose) {
      std::cout << tA << ", " << tB << std::endl;
    }

    tA = clamp(tA, 0., 1.);
    tB = clamp(tB, 0., 1.);

    embeddedCurve.push_back(std::make_pair(SurfacePoint{he, tA}, tB));
  }
  if (verbose) {
    std::cout << std::endl;
  }

  embeddedCurve.push_back(std::make_pair(vEnd, 1));

  return embeddedCurve;
}

// Compute the new normal coordinates for an inserted vertex v in face f given
// that v's barycentric coordinates and the barycentric coordinates of the
// crossings along f's edges
std::array<int, 3> computeVertexInsertionCrossingCounts(Vector3 bary,
                                                        const std::array<std::vector<double>, 3>& boundaryCrossings) {

  // Pick some arbitrary positions for the triangle's 3 vertices
  Vector2 vi{0, 0};
  Vector2 vj{1, 0};
  Vector2 vk{0, 1};

  std::array<Vector2, 3> positions{vi, vj, vk};

  size_t nij = boundaryCrossings[0].size();
  size_t njk = boundaryCrossings[1].size();
  size_t nki = boundaryCrossings[2].size();

  // ns[i] is the normal coordinate of the edge opposite vertex i
  std::array<size_t, 3> ns{njk, nki, nij};

  int ci = cornerCoord(njk, nki, nij);
  int cj = cornerCoord(nki, nij, njk);
  int ck = cornerCoord(nij, njk, nki);

  std::array<int, 3> cs{ci, cj, ck};

  Vector2 newPt = bary.x * vi + bary.y * vj + bary.z * vk;

  // Get the barycentric coordinate of the iC'th crossing along the edge
  // from vertex iSrc to vertex iDst
  auto crossingCoord = [&](size_t iSrc, size_t iDst, size_t iC) -> double {
    size_t iEdge = 0;

    if (iSrc == 0 && iDst == 1) {
      iEdge = 0;
    } else if (iSrc == 1 && iDst == 2) {
      iEdge = 1;
    } else if (iSrc == 2 && iDst == 0) {
      iEdge = 2;
    } else if (iSrc == 1 && iDst == 0) {
      iEdge = 0;
      iC = nij - iC - 1;
    } else if (iSrc == 2 && iDst == 1) {
      iEdge = 1;
      iC = njk - iC - 1;
    } else {
      iEdge = 2;
      iC = nki - iC - 1;
    }

    return boundaryCrossings[iEdge][iC];
  };

  auto next = [](size_t iV) -> size_t { return (iV + 1) % 3; };
  auto prev = [](size_t iV) -> size_t { return (iV + 2) % 3; };

  // Get the line (represented by the 2D coordinates of its endpoints)
  // that is the iLine'th curve cutting off corner iV
  auto cornerLine = [&](size_t iV, size_t iLine) -> std::array<Vector2, 2> {
    size_t pV = prev(iV);
    size_t nV = next(iV);

    double tPrev = crossingCoord(iV, pV, iLine);
    double tNext = crossingCoord(iV, nV, iLine);

    Vector2 src = (1 - tPrev) * positions[pV] + tPrev * positions[iV];
    Vector2 dst = (1 - tNext) * positions[iV] + tNext * positions[nV];

    return {src, dst};
  };

  // Return true if pt is to the left of the oriented line, false
  // otherwise
  auto outside = [&](const std::array<Vector2, 2>& line, Vector2 pt) -> bool {
    // HACK: if line degenerates to a point, the cross product computation
    // gets very unstable. So if line is too short, we always say that pt is
    // outside the line

    Vector2 lineDir = line[1] - line[0];
    if (lineDir.norm2() < 1e-6) {
      return true;
    }


    Vector2 ptDir = pt - line[0];

    return cross(ptDir, lineDir) <= 0;
  };

  // Count how many curves the line from newPt to corner iCorner crosses
  auto cornerSection = [&](const size_t iCorner) -> int {
    int crossings = 0;
    while (crossings < cs[iCorner] && outside(cornerLine(iCorner, crossings), newPt)) {
      crossings++;
    }
    return crossings;
  };

  if (ci >= 0 && cj >= 0 && ck >= 0) {
    // Triforce
    int iCrossings = cornerSection(0);
    int jCrossings = cornerSection(1);
    int kCrossings = cornerSection(2);

    int iSlack = ci - iCrossings;
    int jSlack = cj - jCrossings;
    int kSlack = ck - kCrossings;

    // Only one vertex can have slack (i.e. crossings > c); on bad meshes,
    // floating point errors can (very rarely) cause multiple vertices to
    // have slack. In this case, pick the biggest
    if (iSlack >= jSlack && iSlack >= kSlack) {
      // This case also happens when all slacks are 0, ie iCrossings = ci,
      // for all ijk. In that case, iSlack is 0, so this is fine
      jCrossings = iSlack + cj;
      kCrossings = iSlack + ck;
    } else if (jSlack >= kSlack && jSlack >= iSlack) {
      kCrossings = jSlack + ck;
      iCrossings = jSlack + ci;
    } else if (kSlack >= iSlack && kSlack >= jSlack) {
      iCrossings = kSlack + ci;
      jCrossings = kSlack + cj;
    }

    return {iCrossings, jCrossings, kCrossings};
  } else {
    // Fan

    // Compute the normal coordinates, assuming fan vertex comes first
    auto fanCrossings = [&](size_t vFan, size_t vNext, size_t vPrev) -> std::tuple<int, int, int> {
      int fanCrossings = 0;
      int nextCrossings = cornerSection(vNext);
      int prevCrossings = cornerSection(vPrev);

      // Number of curves emanating from of vFan
      int nEmanating = strictDegree(ns[vFan], ns[vNext], ns[vPrev]);

      if (prevCrossings < cs[vPrev]) {
        int slack = cs[vPrev] - prevCrossings;
        nextCrossings += slack + nEmanating;
        fanCrossings += slack;
      } else if (nextCrossings < cs[vNext]) {
        int slack = cs[vNext] - nextCrossings;
        fanCrossings += slack;
        prevCrossings += slack + nEmanating;
      } else {

        int crossingIndex = 0;
        while (crossingIndex < nEmanating) {
          double tCross = boundaryCrossings[(vFan + 1) % 3][cs[vNext] + crossingIndex];

          Vector2 pCross = (1 - tCross) * positions[vNext] + tCross * positions[vPrev];
          if (outside({positions[vFan], pCross}, newPt)) {
            crossingIndex++;
          } else {
            break;
          }
        }

        nextCrossings += crossingIndex;
        prevCrossings += nEmanating - crossingIndex;
      }

      return std::tuple<int, int, int>(fanCrossings, nextCrossings, prevCrossings);
    };

    int iCrossings, jCrossings, kCrossings;
    if (ci < 0) {
      std::tie(iCrossings, jCrossings, kCrossings) = fanCrossings(0, 1, 2);
    } else if (cj < 0) {
      std::tie(jCrossings, kCrossings, iCrossings) = fanCrossings(1, 2, 0);
    } else {
      std::tie(kCrossings, iCrossings, jCrossings) = fanCrossings(2, 0, 1);
    }
    return {iCrossings, jCrossings, kCrossings};
  }
}

} // namespace surface
} // namespace geometrycentral
