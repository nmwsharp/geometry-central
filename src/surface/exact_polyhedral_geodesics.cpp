#include "geometrycentral/surface/exact_polyhedral_geodesics.h"

#include <queue>
#include <tuple>

using std::cout;
using std::endl;

namespace geometrycentral {
namespace surface {

ExactPolyhedralGeodesics::ExactPolyhedralGeodesics(EdgeLengthGeometry* geom_) : geom(geom_) {
  mesh = geom->mesh;
  clear();
}

void ExactPolyhedralGeodesics::addSource(Vertex v) { srcVerts.push_back(v); }

void ExactPolyhedralGeodesics::clear() {
  while (not winQ.empty()) winQ.pop();
  while (not pseudoSrcQ.empty()) pseudoSrcQ.pop();

  splitInfos = HalfedgeData<SplitInfo>(*mesh);
  vertInfos = VertexData<VertInfo>(*mesh);

  for (Halfedge he : mesh->halfedges()) {
    splitInfos[he].dist = std::numeric_limits<double>::infinity();
    splitInfos[he].x = std::numeric_limits<double>::infinity();
    splitInfos[he].pseudoSrc = Vertex();
    splitInfos[he].src = Vertex();
    splitInfos[he].level = -1;
  }

  for (Vertex v : mesh->vertices()) {
    vertInfos[v].birthTime = -1;
    vertInfos[v].dist = std::numeric_limits<double>::infinity();
    vertInfos[v].enterHalfedge = Halfedge();
    vertInfos[v].isSource = false;
    vertInfos[v].pseudoSrc = Vertex();
    vertInfos[v].src = Vertex();
  }

  srcVerts.clear();

  storedWindows.clear();
  keptFaces.clear();
  allWindows = HalfedgeData<std::list<Window>>(*mesh);
  keptAllWindows = false;

  numOfWinGen = 0;
  maxWinQSize = 0;
  maxPseudoQSize = 0;
  totalCalcVertNum = 0;
  geodesicRadius = std::numeric_limits<double>::infinity();
  geodesicRadiusReached = false;
}

// TODO: allow arbitrary surface points as input
void ExactPolyhedralGeodesics::initialize() {
  // initialize
  for (auto v : srcVerts) {
    for (Halfedge he : v.outgoingHalfedges()) {
      Halfedge oppHE = he.next();

      Window win;
      win.halfedge = oppHE;
      win.b0 = 0.;
      win.b1 = geom->edgeLengths[oppHE.edge()];
      win.d0 = geom->edgeLengths[he.edge()];
      win.d1 = geom->edgeLengths[oppHE.next().edge()];
      win.pseudoSrcDist = 0.;
      win.computeMinDist();
      win.src = v;
      win.pseudoSrc = v;
      win.pseudoSrcBirthTime = 0;
      win.level = 0;
      winQ.push(win);
      numOfWinGen += 1;

      Vertex oppV = he.twin().vertex();
      if (geom->edgeLengths[he.edge()] < vertInfos[oppV].dist) {
        vertInfos[oppV].birthTime = 0;
        vertInfos[oppV].dist = geom->edgeLengths[he.edge()];
        vertInfos[oppV].enterHalfedge = he.twin();
        vertInfos[oppV].src = v;
        vertInfos[oppV].pseudoSrc = v;

        // only need to create a new pseudoWin if the vertex is hyperbolic
        if (geom->vertexAngleDefects[oppV] > 0.) continue;

        PseudoWindow pseudoWin;
        pseudoWin.v = oppV;
        pseudoWin.dist = geom->edgeLengths[he.edge()];
        pseudoWin.src = v;
        pseudoWin.pseudoSrc = v;
        pseudoWin.pseudoSrcBirthTime = vertInfos[oppV].birthTime;
        pseudoWin.level = 0;
        pseudoSrcQ.push(pseudoWin);
      }
    }
    vertInfos[v].birthTime = 0;
    vertInfos[v].dist = 0.;
    vertInfos[v].enterHalfedge = Halfedge();
    vertInfos[v].isSource = true;
    vertInfos[v].src = v;
    vertInfos[v].pseudoSrc = v;
  }
}

bool ExactPolyhedralGeodesics::isValidWindow(const Window& win, bool isLeftChild) {
  if (win.b1 <= win.b0) return false;
  // apply ICH's filter
  Vertex v1 = win.halfedge.vertex();
  Vertex v2 = win.halfedge.twin().vertex();
  Vertex v3 = win.halfedge.next().twin().vertex();
  double l0 = geom->edgeLengths[win.halfedge.edge()];
  double l1 = geom->edgeLengths[win.halfedge.next().edge()];
  double l2 = geom->edgeLengths[win.halfedge.next().next().edge()];

  Vector2 p1{0.0, 0.0}, p2{l0, 0.0}, p3;
  p3.x = (l2 * l2 + l0 * l0 - l1 * l1) / (2.0 * l0);
  p3.y = sqrt(fabs(l2 * l2 - p3.x * p3.x));

  Vector2 A{win.b0, 0.0}, B{win.b1, 0.0};
  Vector2 src2D = win.flattenedSrc();
  if (win.pseudoSrcDist + norm(src2D - B) > vertInfos[v1].dist + win.b1 &&
      (win.pseudoSrcDist + norm(src2D - B)) / (vertInfos[v1].dist + win.b1) - 1.0 > 0.0) {
    return false;
  }
  if (win.pseudoSrcDist + norm(src2D - A) > vertInfos[v2].dist + l0 - win.b0 &&
      (win.pseudoSrcDist + norm(src2D - A)) / (vertInfos[v2].dist + l0 - win.b0) - 1.0 > 0.0) {
    return false;
  }
  if (isLeftChild) {
    if (win.pseudoSrcDist + norm(src2D - A) > vertInfos[v3].dist + norm(p3 - A) &&
        (win.pseudoSrcDist + norm(src2D - A)) / (vertInfos[v3].dist + norm(p3 - A)) - 1.0 > 0.0) {
      return false;
    }
  } else {
    if (win.pseudoSrcDist + norm(src2D - B) > vertInfos[v3].dist + norm(p3 - B) &&
        (win.pseudoSrcDist + norm(src2D - B)) / (vertInfos[v3].dist + norm(p3 - B)) - 1.0 > REL_ERR) {
      return false;
    }
  }
  return true;
}

void ExactPolyhedralGeodesics::buildWindow(const Window& pWin, Halfedge& he, double t0, double t1, const Vector2& v0,
                                           const Vector2& v1, Window& win) {
  Vector2 src2D = pWin.flattenedSrc();
  win.halfedge = he;
  win.b0 = (1 - t0) * geom->edgeLengths[he.edge()];
  win.b1 = (1 - t1) * geom->edgeLengths[he.edge()];
  win.d0 = norm(src2D - (t0 * v0 + (1 - t0) * v1));
  win.d1 = norm(src2D - (t1 * v0 + (1 - t1) * v1));
  win.pseudoSrcDist = pWin.pseudoSrcDist;
  win.computeMinDist();
  win.src = pWin.src;
  win.pseudoSrc = pWin.pseudoSrc;
  win.pseudoSrcBirthTime = pWin.pseudoSrcBirthTime;
  win.level = pWin.level + 1;
}

double ExactPolyhedralGeodesics::intersect(const Vector2& v0, const Vector2& v1, const Vector2& p0, const Vector2& p1) {
  double a00 = p0.x - p1.x, a01 = v1.x - v0.x, b0 = v1.x - p1.x;
  double a10 = p0.y - p1.y, a11 = v1.y - v0.y, b1 = v1.y - p1.y;
  return (b0 * a11 - b1 * a01) / (a00 * a11 - a10 * a01);
}

void ExactPolyhedralGeodesics::propogateWindow(const Window& win) {
  Halfedge he0 = win.halfedge.twin();
  if (he0 == Halfedge()) return;

  Halfedge he1 = he0.next();
  Halfedge he2 = he1.next();

  Vertex oppV = he1.twin().vertex();
  Vector2 src2D = win.flattenedSrc();
  Vector2 left{win.b0, 0.}, right{win.b1, 0.};

  double l0 = geom->edgeLengths[he0.edge()];
  double l1 = geom->edgeLengths[he1.edge()];
  double l2 = geom->edgeLengths[he2.edge()];

  Vector2 v0{0., 0.}, v1{l0, 0.}, v2;
  v2.x = (l1 * l1 + l0 * l0 - l2 * l2) / (2. * l0);
  v2.y = -sqrt(fabs(l1 * l1 - v2.x * v2.x));

  double interX = v2.x - v2.y * (v2.x - src2D.x) / (v2.y - src2D.y);
  Window leftChildWin, rightChildWin;
  bool hasLeftChild = true, hasRightChild = true;

  // only generate right window
  if (interX <= left.x) {
    hasLeftChild = false;
    double t0 = intersect(src2D, left, v2, v1);
    double t1 = intersect(src2D, right, v2, v1);

    buildWindow(win, he2, t0, t1, v2, v1, rightChildWin);
    if (not isValidWindow(rightChildWin, false)) hasRightChild = false;
  }

  // only generate left window
  else if (interX >= right.x) {
    hasRightChild = false;
    double t0 = intersect(src2D, left, v0, v2);
    double t1 = intersect(src2D, right, v0, v2);
    buildWindow(win, he1, t0, t1, v0, v2, leftChildWin);
    if (not isValidWindow(leftChildWin, true)) hasLeftChild = false;
  }

  // generate both left and right windows
  else {
    double directDist = norm(v2 - src2D);
    // if( directDist + win.pseudoSrcDist > splitInfos[he0].dist
    //     && (directDist + win.pseudoSrcDist)/splitInfos[he0].dist - 1. > REL_ERR )
    // {
    //   hasLeftChild = splitInfos[he0].x < interX;
    //   hasRightChild = not hasLeftChild;
    // }
    // else {
    {
      if (directDist + win.pseudoSrcDist < splitInfos[he0].dist) {
        splitInfos[he0].dist = directDist + win.pseudoSrcDist;
        splitInfos[he0].pseudoSrc = win.pseudoSrc;
        splitInfos[he0].src = win.src;
        splitInfos[he0].level = win.level;
        splitInfos[he0].x = l0 - interX;
      }

      if (directDist + win.pseudoSrcDist < vertInfos[oppV].dist) {
        if (vertInfos[oppV].dist == std::numeric_limits<double>::infinity()) {
          ++totalCalcVertNum;
        }
        if (directDist + win.pseudoSrcDist >= geodesicRadius) {
          geodesicRadiusReached = true;
        }

        vertInfos[oppV].birthTime++;
        vertInfos[oppV].dist = directDist + win.pseudoSrcDist;
        vertInfos[oppV].enterHalfedge = he0;
        vertInfos[oppV].src = win.src;
        vertInfos[oppV].pseudoSrc = win.pseudoSrc;

        if (geom->vertexAngleDefects[oppV] < 0.) {
          PseudoWindow pseudoWin;
          pseudoWin.v = oppV;
          pseudoWin.dist = vertInfos[oppV].dist;
          pseudoWin.src = win.src;
          pseudoWin.pseudoSrc = win.pseudoSrc;
          pseudoWin.pseudoSrcBirthTime = vertInfos[oppV].birthTime;
          pseudoWin.level = win.level + 1;
          pseudoSrcQ.push(pseudoWin);
        }
      }
    }
    if (hasLeftChild) {
      double t0 = intersect(src2D, left, v0, v2);
      buildWindow(win, he1, t0, 0., v0, v2, leftChildWin);
      if (not isValidWindow(leftChildWin, true)) hasLeftChild = false;
    }
    if (hasRightChild) {
      double t1 = intersect(src2D, right, v2, v1);
      buildWindow(win, he2, 1., t1, v2, v1, rightChildWin);
      if (not isValidWindow(rightChildWin, false)) hasRightChild = false;
    }
  }
  if (hasLeftChild) {
    numOfWinGen += 1;
    winQ.push(leftChildWin);
  }
  if (hasRightChild) {
    numOfWinGen += 1;
    winQ.push(rightChildWin);
  }
}

void ExactPolyhedralGeodesics::generateSubWinsForPseudoSrc(const PseudoWindow& pseudoWin) {
  Halfedge startHe, endHe;

  if (vertInfos[pseudoWin.v].enterHalfedge == Halfedge() && vertInfos[pseudoWin.v].birthTime != -1) {
    startHe = pseudoWin.v.halfedge();
    endHe = startHe;
  } else if (vertInfos[pseudoWin.v].enterHalfedge.vertex() == pseudoWin.v)
    generateSubWinsForPseudoSrcFromPseudoSrc(pseudoWin, startHe, endHe);
  else if (vertInfos[pseudoWin.v].enterHalfedge.next().twin().vertex() == pseudoWin.v)
    generateSubWinsForPseudoSrcFromWindow(pseudoWin, startHe, endHe);
  else
    assert(false);

  // generate windows
  do {
    Window win;
    win.halfedge = startHe.next();
    win.b0 = 0.;
    win.b1 = geom->edgeLengths[win.halfedge.edge()];
    win.d0 = geom->edgeLengths[startHe.edge()];
    win.d1 = geom->edgeLengths[win.halfedge.next().edge()];
    win.pseudoSrcDist = pseudoWin.dist;
    win.computeMinDist();
    win.src = pseudoWin.src;
    win.pseudoSrc = pseudoWin.v;
    win.pseudoSrcBirthTime = pseudoWin.pseudoSrcBirthTime;
    win.level = pseudoWin.level + 1;
    winQ.push(win);
    numOfWinGen += 1;

    startHe = startHe.next().next().twin();
  } while (startHe != endHe);

  Vertex vx = pseudoWin.v;
  for (Halfedge he : vx.outgoingHalfedges()) {
    Vertex oppV = he.twin().vertex();

    // if( geom->vertexAngleDefects[oppV] > 0. ) continue;
    if (vertInfos[oppV].dist < pseudoWin.dist + geom->edgeLengths[he.edge()]) continue;

    if (vertInfos[oppV].dist == std::numeric_limits<double>::infinity()) {
      totalCalcVertNum += 1;
    }
    if (pseudoWin.dist + geom->edgeLengths[he.edge()] >= geodesicRadius) {
      geodesicRadiusReached = true;
    }

    vertInfos[oppV].dist = pseudoWin.dist + geom->edgeLengths[he.edge()];
    vertInfos[oppV].birthTime += 1;
    vertInfos[oppV].enterHalfedge = he.twin();
    vertInfos[oppV].src = pseudoWin.src;
    vertInfos[oppV].pseudoSrc = pseudoWin.pseudoSrc;

    PseudoWindow childPseudoWin;
    childPseudoWin.v = oppV;
    childPseudoWin.dist = vertInfos[oppV].dist;
    childPseudoWin.src = pseudoWin.src;
    childPseudoWin.pseudoSrc = pseudoWin.v;
    childPseudoWin.pseudoSrcBirthTime = vertInfos[oppV].birthTime;
    childPseudoWin.level = pseudoWin.level;
    pseudoSrcQ.push(childPseudoWin);
  }
}

void ExactPolyhedralGeodesics::generateSubWinsForPseudoSrcFromWindow(const PseudoWindow& pseudoWin, Halfedge& startHe,
                                                                     Halfedge& endHe) {
  Halfedge he0 = vertInfos[pseudoWin.v].enterHalfedge;
  Halfedge he1 = he0.next();
  Halfedge he2 = he1.next();

  double l0 = geom->edgeLengths[he0.edge()];
  double l1 = geom->edgeLengths[he1.edge()];
  double l2 = geom->edgeLengths[he2.edge()];

  Vertex pseudoSrc = pseudoWin.src;
  Vector2 enterPt;
  enterPt.x = l0 - splitInfos[he0].x;
  enterPt.y = 0.;

  Vector2 v0{0., 0.}, v1{l0, 0.}, v2;
  v2.x = (l1 * l1 + l0 * l0 - l2 * l2) / (2. * l0);
  v2.y = -sqrt(fabs(l1 * l1 - v2.x * v2.x));

  double angleFromLeft = dot(enterPt - v2, v0 - v2) / norm(enterPt - v2) / l1;
  double angleFromRight = dot(enterPt - v2, v1 - v2) / norm(enterPt - v2) / l2;
  if (angleFromLeft > 1.)
    angleFromLeft = 1.;
  else if (angleFromLeft < -1.)
    angleFromLeft = -1.;
  if (angleFromRight > 1.)
    angleFromRight = 1.;
  else if (angleFromRight < -1.)
    angleFromRight = -1.;
  angleFromLeft = acos(angleFromLeft);
  angleFromRight = acos(angleFromRight);

  startHe = Halfedge();
  endHe = Halfedge();
  Halfedge currHe = he1.twin();
  while (angleFromLeft < M_PI && currHe != Halfedge()) {
    Halfedge oppHe = currHe.next();
    Halfedge nextHe = oppHe.next();
    double L0 = geom->edgeLengths[currHe.edge()];
    double L1 = geom->edgeLengths[nextHe.edge()];
    double L2 = geom->edgeLengths[oppHe.edge()];
    double currAngle = (L0 * L0 + L1 * L1 - L2 * L2) / (2. * L0 * L1);
    if (currAngle > 1.)
      currAngle = 1.;
    else if (currAngle < -1.)
      currAngle = -1.;
    currAngle = acos(currAngle);
    angleFromLeft += currAngle;
    currHe = nextHe.twin();
  }
  if (currHe != Halfedge()) startHe = currHe.twin().next();

  currHe = he2.twin();
  while (angleFromRight < M_PI && currHe != Halfedge()) {
    Halfedge nextHe = currHe.next();
    Halfedge oppHe = nextHe.next();
    double L0 = geom->edgeLengths[currHe.edge()];
    double L1 = geom->edgeLengths[nextHe.edge()];
    double L2 = geom->edgeLengths[oppHe.edge()];
    double currAngle = (L0 * L0 + L1 * L1 - L2 * L2) / (2. * L0 * L1);
    if (currAngle > 1.)
      currAngle = 1.;
    else if (currAngle < -1.)
      currAngle = -1.;
    currAngle = acos(currAngle);
    angleFromRight += currAngle;
    currHe = nextHe.twin();
  }
  if (currHe != Halfedge()) endHe = currHe.twin().next().next().twin();
}

void ExactPolyhedralGeodesics::generateSubWinsForPseudoSrcFromPseudoSrc(const PseudoWindow& pseudoWin,
                                                                        Halfedge& startHe, Halfedge& endHe) {
  Vertex pseudoSrc = pseudoWin.v;
  double angleFromLeft = 0., angleFromRight = 0.;
  startHe = Halfedge();
  endHe = Halfedge();
  Halfedge currHe = vertInfos[pseudoWin.v].enterHalfedge;
  while (angleFromLeft < M_PI && currHe != Halfedge()) {
    Halfedge oppHe = currHe.next();
    Halfedge nextHe = oppHe.next();
    double L0 = geom->edgeLengths[currHe.edge()];
    double L1 = geom->edgeLengths[nextHe.edge()];
    double L2 = geom->edgeLengths[oppHe.edge()];
    double currAngle = (L0 * L0 + L1 * L1 - L2 * L2) / (2. * L0 * L1);
    if (currAngle > 1.)
      currAngle = 1.;
    else if (currAngle < -1.)
      currAngle = -1.;
    currAngle = acos(currAngle);
    angleFromLeft += currAngle;
    currHe = nextHe.twin();
  }
  if (currHe != Halfedge()) startHe = currHe.twin().next();

  currHe = vertInfos[pseudoWin.v].enterHalfedge.twin();
  while (angleFromRight < M_PI && currHe != Halfedge()) {
    Halfedge nextHe = currHe.next();
    Halfedge oppHe = nextHe.next();
    double L0 = geom->edgeLengths[currHe.edge()];
    double L1 = geom->edgeLengths[nextHe.edge()];
    double L2 = geom->edgeLengths[oppHe.edge()];
    double currAngle = (L0 * L0 + L1 * L1 - L2 * L2) / (2. * L0 * L1);
    if (currAngle > 1.)
      currAngle = 1.;
    else if (currAngle < -1.)
      currAngle = -1.;
    currAngle = acos(currAngle);
    angleFromRight += currAngle;
    currHe = nextHe.twin();
  }
  if (currHe != Halfedge()) endHe = currHe.twin().next().next().twin();
}

VertexData<double> ExactPolyhedralGeodesics::computeDistance() {
  geom->requireEdgeLengths();
  geom->requireVertexAngleDefects();
  geom->requireVertexIndices();

  initialize();

  while (not winQ.empty() || not pseudoSrcQ.empty()) {
    maxWinQSize = fmax(maxWinQSize, winQ.size());
    maxPseudoQSize = fmax(maxPseudoQSize, pseudoSrcQ.size());

    while (not winQ.empty() && winQ.top().pseudoSrc != Vertex() &&
           winQ.top().pseudoSrcBirthTime != vertInfos[winQ.top().pseudoSrc].birthTime)
      winQ.pop();

    while (not pseudoSrcQ.empty() && winQ.top().pseudoSrc != Vertex() &&
           (int)pseudoSrcQ.top().pseudoSrcBirthTime != vertInfos[pseudoSrcQ.top().v].birthTime)
      pseudoSrcQ.pop();

    if (not winQ.empty() && (pseudoSrcQ.empty() || winQ.top().minDist < pseudoSrcQ.top().dist)) {
      Window win = winQ.top();
      winQ.pop();
      if (win.level > (int)mesh->nFaces()) continue;

      Halfedge twin = win.halfedge.twin();
      propogateWindow(win);
    } else if (not pseudoSrcQ.empty() && (winQ.empty() || winQ.top().minDist >= pseudoSrcQ.top().dist)) {
      PseudoWindow pseudoWin = pseudoSrcQ.top();
      pseudoSrcQ.pop();
      if (pseudoWin.level >= (unsigned)mesh->nFaces()) continue;
      generateSubWinsForPseudoSrc(pseudoWin);
    }

    if (geodesicRadiusReached) break;
  }

  VertexData<double> dists(*mesh);
  for (Vertex v : mesh->vertices()) {
    dists[v] = vertInfos[v].dist;
  }
  return dists;
}

} // namespace surface
} // namespace geometrycentral
