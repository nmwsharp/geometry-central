#pragma once

#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/utilities/utilities.h"

#include <unordered_map>
#include <vector>


namespace geometrycentral {
namespace surface {


// Return all vertices and their geodesic distance within the given radius threshold. Uses sparse data structures, so is
// efficient for small queries on large meshes. However, performs greedy triangle unfolding to compute distance, so has
// exponential worst-case complexity. Should be used for small queries only!
std::unordered_map<Vertex, double> vertexGeodesicDistanceWithinRadius(IntrinsicGeometryInterface& geom,
                                                                      Vertex centerVert, double radius);


/*
const double REL_ERR = 1e-8;

struct Window {
  Halfedge halfedge;
  double b0, b1;
  double d0, d1;

  double pseudoSrcDist, minDist;
  Vertex src, pseudoSrc;
  int pseudoSrcBirthTime;
  int level;

  bool operator<(const Window& right) const { return minDist > right.minDist; }
  void computeMinDist() {
    double wLen = b1 - b0;
    double xProj = (d0 * d0 + wLen * wLen - d1 * d1) / (2. * wLen);
    if (xProj < 0.)
      minDist = d0 + pseudoSrcDist;
    else if (xProj > wLen)
      minDist = d1 + pseudoSrcDist;
    else
      minDist = sqrt(fabs(d0 * d0 - xProj * xProj)) + pseudoSrcDist;
  }
  Vector2 flattenedSrc() const {
    Vector2 src2D;
    double wLen = b1 - b0;
    src2D.x = (d0 * d0 + wLen * wLen - d1 * d1) / (2. * wLen);
    src2D.y = sqrt(fabs(d0 * d0 - src2D.x * src2D.x));
    src2D.x += b0;
    return src2D;
  }
};

struct PseudoWindow {
  Vertex v;
  double dist;
  Vertex src, pseudoSrc;
  unsigned pseudoSrcBirthTime;
  unsigned level;

  bool operator<(const PseudoWindow& right) const { return dist > right.dist; }
};

struct SplitInfo {
  double dist;
  Vertex pseudoSrc, src;
  unsigned level;
  double x;

  SplitInfo() {
    dist = std::numeric_limits<double>::infinity();
    pseudoSrc = Vertex();
    src = Vertex();
    level = -1;
    x = std::numeric_limits<double>::infinity();
  }
};

struct VertInfo {
  int birthTime;
  double dist;
  bool isSource;
  Halfedge enterHalfedge;
  Vertex pseudoSrc, src;

  VertInfo() {
    birthTime = -1;
    dist = std::numeric_limits<double>::infinity();
    isSource = false;
    enterHalfedge = Halfedge();
    pseudoSrc = Vertex();
    src = Vertex();
  }
};

struct GeodesicKeyPoint {
  bool isVertex;
  unsigned id;
  double pos;
};

class ExactPolyhedralGeodesics {

public:
  ExactPolyhedralGeodesics(EdgeLengthGeometry* geom_);

  void addSource(Vertex v);
  VertexData<double> computeDistance(Vertex v);
  VertexData<double> computeDistance();

private:
  HalfedgeMesh* mesh;
  EdgeLengthGeometry* geom;

  std::vector<Vertex> srcVerts;
  HalfedgeData<SplitInfo> splitInfos;
  VertexData<VertInfo> vertInfos;
  std::priority_queue<Window> winQ;
  std::priority_queue<PseudoWindow> pseudoSrcQ;

  std::vector<Window> storedWindows;
  std::vector<Face> keptFaces;
  HalfedgeData<std::list<Window>> allWindows;
  bool keptAllWindows = false;

  // stats
  int numOfWinGen;
  int maxWinQSize, maxPseudoQSize;
  int totalCalcVertNum;
  double geodesicRadius;
  bool geodesicRadiusReached;

  void clear();
  void initialize();
  void propogateWindow(const Window& win);
  void generateSubWinsForPseudoSrc(const PseudoWindow& pseudoWin);
  void generateSubWinsForPseudoSrcFromWindow(const PseudoWindow& pseudoWin, Halfedge& startHe, Halfedge& endHe);
  void generateSubWinsForPseudoSrcFromPseudoSrc(const PseudoWindow& pseudoWin, Halfedge& startHe, Halfedge& endHe);
  void buildWindow(const Window& pWin, Halfedge& he, double t0, double t1, const Vector2& v0, const Vector2& v1,
                   Window& win);
  bool isValidWindow(const Window& win, bool isLeftChild);
  double intersect(const Vector2& v0, const Vector2& v1, const Vector2& p0, const Vector2& p1);
};

*/

} // namespace surface
} // namespace geometrycentral
