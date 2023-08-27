// Copyright (C) 2008 Danil Kirsanov, MIT License
// (Modified to work in geometry-central. Original code can be found here: https://code.google.com/p/geodesic/)

#include "geometrycentral/surface/exact_geodesic_helpers.h"

namespace geometrycentral {
namespace surface {
double Interval::signal(double x) const {
  assert(x >= 0.0 && x <= m_edge_length);

  if (m_d == GEODESIC_INF) {
    return GEODESIC_INF;
  } else {
    double dx = x - m_pseudo_x;
    if (m_pseudo_y == 0.0) {
      return m_d + std::abs(dx);
    } else {
      return m_d + sqrt(dx * dx + m_pseudo_y * m_pseudo_y);
    }
  }
}

double Interval::max_distance(double end) const {
  if (m_d == GEODESIC_INF) {
    return GEODESIC_INF;
  } else {
    double a = std::abs(m_start - m_pseudo_x);
    double b = std::abs(end - m_pseudo_x);

    return a > b ? m_d + sqrt(a * a + m_pseudo_y * m_pseudo_y) : m_d + sqrt(b * b + m_pseudo_y * m_pseudo_y);
  }
}

void Interval::compute_min_distance(double stop) {
  assert(stop > m_start);

  if (m_d == GEODESIC_INF) {
    m_min = GEODESIC_INF;
  } else if (m_start > m_pseudo_x) {
    m_min = signal(m_start);
  } else if (stop < m_pseudo_x) {
    m_min = signal(stop);
  } else {
    assert(m_pseudo_y <= 0);
    m_min = m_d - m_pseudo_y;
  }
}

bool Interval::operator()(interval_pointer const x, interval_pointer const y) const {
  if (x->min() != y->min()) {
    return x->min() < y->min();
  } else if (x->start() != y->start()) {
    return x->start() < y->start();
  } else {
    return x->edge() < y->edge();
  }
}

double Interval::stop() const { return m_next ? m_next->start() : m_edge_length; }

double Interval::hypotenuse(double a, double b) const { return sqrt(a * a + b * b); }

void Interval::find_closest_point(double const rs, double const hs, double& r, double& d_out) const {
  if (m_d == GEODESIC_INF) {
    r = GEODESIC_INF;
    d_out = GEODESIC_INF;
    return;
  }

  double hc = -m_pseudo_y;
  double rc = m_pseudo_x;
  double end = stop();

  double local_epsilon = SMALLEST_INTERVAL_RATIO * m_edge_length;
  if (std::abs(hs + hc) < local_epsilon) {
    if (rs <= m_start) {
      r = m_start;
      d_out = signal(m_start) + std::abs(rs - m_start);
    } else if (rs >= end) {
      r = end;
      d_out = signal(end) + fabs(end - rs);
    } else {
      r = rs;
      d_out = signal(rs);
    }
  } else {
    double ri = (rs * hc + hs * rc) / (hs + hc);

    if (ri < m_start) {
      r = m_start;
      d_out = signal(m_start) + hypotenuse(m_start - rs, hs);
    } else if (ri > end) {
      r = end;
      d_out = signal(end) + hypotenuse(end - rs, hs);
    } else {
      r = ri;
      d_out = m_d + hypotenuse(rc - rs, hc + hs);
    }
  }
}

void Interval::initialize(IntrinsicGeometryInterface& geom, Edge edge, double edge_length, SurfacePoint* source,
                          unsigned source_index) {
  m_next = nullptr;
  // m_geodesic_previous = nullptr;
  m_direction = DirectionType::UNDEFINED_DIRECTION;
  m_edge = edge;
  m_edge_length = edge_length;
  m_source_index = source_index;

  m_start = 0.0;
  // m_stop = m_edge_length;
  if (!source) {
    m_d = GEODESIC_INF;
    m_min = GEODESIC_INF;
    return;
  }
  m_d = 0;

  if (source->type == SurfacePointType::Vertex) {
    if (source->vertex == edge.halfedge().tailVertex()) {
      m_pseudo_x = 0.0;
      m_pseudo_y = 0.0;
      m_min = 0.0;
      return;
    } else if (source->vertex == edge.halfedge().tipVertex()) {
      m_pseudo_x = stop();
      m_pseudo_y = 0.0;
      m_min = 0.0;
      return;
    }
  }

  std::tie(m_pseudo_x, m_pseudo_y) = exactgeodesic::compute_local_coordinates(geom, edge, *source);
  // edge->local_coordinates(source, m_pseudo_x, m_pseudo_y);
  m_pseudo_y = -m_pseudo_y;

  compute_min_distance(stop());
}

IntervalList::IntervalList() { m_first = nullptr; }
void IntervalList::clear() { m_first = nullptr; };

void IntervalList::initialize(Edge e) {
  m_edge = e;
  m_first = nullptr;
};

interval_pointer IntervalList::covering_interval(double offset) {
  assert(m_first == nullptr || (offset >= 0.0 && offset <= m_first->edge_length()));

  interval_pointer p = m_first;
  while (p && p->stop() < offset) {
    p = p->next();
  }

  return p; // && p->start() <= offset ? p : nullptr;
};

const_interval_pointer IntervalList::covering_interval(double offset) const {
  assert(m_first == nullptr || (offset >= 0.0 && offset <= m_first->edge_length()));

  const_interval_pointer p = m_first;
  while (p && p->stop() < offset) {
    p = p->next();
  }

  return p; // && p->start() <= offset ? p : nullptr;
};

void IntervalList::find_closest_point(IntrinsicGeometryInterface& geom, const SurfacePoint point, double& offset,
                                      double& distance, const_interval_pointer& interval, bool verbose) const {
  const_interval_pointer p = m_first;
  distance = GEODESIC_INF;
  interval = nullptr;

  double x, y;
  std::tie(x, y) = exactgeodesic::compute_local_coordinates(geom, m_edge, point);
  if (verbose) {
    std::cout << "\t local point coordinates are " << x << ", " << y << std::endl;
  }

  while (p) {
    if (p->min() < GEODESIC_INF) {
      double o, d;
      p->find_closest_point(x, y, o, d);
      if (d < distance) {
        distance = d;
        offset = o;
        interval = p;
        // if (verbose) {
        //     std::cout << "\t\t new min interval!" << std::endl;
        // }
      }
    }
    p = p->next();
  }
};

unsigned IntervalList::number_of_intervals() const {
  interval_pointer p = m_first;
  unsigned count = 0;
  while (p) {
    ++count;
    p = p->next();
  }
  return count;
}

interval_pointer IntervalList::last() {
  interval_pointer p = m_first;
  if (p) {
    while (p->next()) {
      p = p->next();
    }
  }
  return p;
}

double IntervalList::signal(double x) const {
  const_interval_pointer interval = covering_interval(x);

  return interval ? interval->signal(x) : GEODESIC_INF;
}
interval_pointer& IntervalList::first() { return m_first; };
Edge& IntervalList::edge() { return m_edge; };

unsigned SurfacePointWithIndex::index() const { return m_index; }
void SurfacePointWithIndex::initialize(const SurfacePoint& p, unsigned index) {
  type = p.type;
  vertex = p.vertex;
  edge = p.edge;
  tEdge = p.tEdge;
  face = p.face;
  faceCoords = p.faceCoords;
  m_index = index;
}

bool SurfacePointWithIndex::operator()(const SurfacePointWithIndex* x, const SurfacePointWithIndex* y) const {
  return compare(*x, *y);
}

bool SurfacePointWithIndex::operator()(const SurfacePointWithIndex& x, const SurfacePointWithIndex* y) const {
  return compare(x, *y);
}

bool SurfacePointWithIndex::operator()(const SurfacePointWithIndex* x, const SurfacePointWithIndex& y) const {
  return compare(*x, y);
}

bool SurfacePointWithIndex::compare(const SurfacePointWithIndex& x, const SurfacePointWithIndex& y) const {
  if (x.type != y.type) {
    return x.type < y.type;
  } else {
    switch (x.type) {
    case SurfacePointType::Vertex:
      return x.vertex.getIndex() < y.vertex.getIndex();
    case SurfacePointType::Edge:
      return x.edge.getIndex() < y.edge.getIndex();
    case SurfacePointType::Face:
      return x.face.getIndex() < y.face.getIndex();
    }

    GC_SAFETY_ASSERT(false, "this should be unreachable");
    return false;
  }
}

SortedSources::sorted_iterator_pair SortedSources::sources(const SurfacePoint& mesh_element) {
  m_search_dummy = SurfacePointWithIndex(mesh_element);

  return equal_range(m_sorted.begin(), m_sorted.end(), &m_search_dummy, m_compare_less);
}

SortedSources::const_sorted_iterator_pair SortedSources::sources(const SurfacePoint& mesh_element) const {
  SurfacePointWithIndex temp_search_dummy = SurfacePointWithIndex(mesh_element);

  return equal_range(m_sorted.begin(), m_sorted.end(), temp_search_dummy, m_compare_less);
}

void SortedSources::initialize(const std::vector<SurfacePoint>& sources) {
  resize(sources.size());
  m_sorted.resize(sources.size());
  for (unsigned i = 0; i < sources.size(); ++i) {
    SurfacePointWithIndex& p = *(begin() + i);

    p.initialize(sources[i], i);
    m_sorted[i] = &p;
  }

  std::sort(m_sorted.begin(), m_sorted.end(), m_compare_less);
}

SurfacePointWithIndex& SortedSources::operator[](unsigned i) {
  assert(i < size());
  return *(begin() + i);
}
const SurfacePointWithIndex& SortedSources::operator[](unsigned i) const {
  assert(i < size());
  return *(begin() + i);
}


namespace exactgeodesic {
double compute_surface_distance(IntrinsicGeometryInterface& geom, const SurfacePoint& p1, const SurfacePoint& p2) {
  Face fShared = sharedFace(p1, p2);

  GC_SAFETY_ASSERT(fShared != Face(), "compute_surface_distance() err: Point " + std::to_string(p1) +
                                          " not adjacent to " + std::to_string(p2));

  SurfacePoint p1InFace = p1.inFace(fShared);
  SurfacePoint p2InFace = p2.inFace(fShared);

  // lengths[i] is the length of the edge opposite the i'th vertex
  Vector3 triangleLengths{geom.edgeLengths[fShared.halfedge().next().edge()],
                          geom.edgeLengths[fShared.halfedge().next().next().edge()],
                          geom.edgeLengths[fShared.halfedge().edge()]};
  Vector3 displacement = p1InFace.faceCoords - p2InFace.faceCoords;
  return displacementLength(displacement, triangleLengths);
}

unsigned compute_closest_vertices(SurfacePoint p, std::vector<Vertex>* storage) {

  switch (p.type) {
  case SurfacePointType::Vertex: {
    if (storage) {
      storage->push_back(p.vertex);
    }
    return 1;
  }
  case SurfacePointType::Edge: {
    Edge e = p.edge;
    Halfedge he = e.halfedge();

    if (storage) {
      storage->push_back(he.tailVertex());
      storage->push_back(he.tipVertex());
      storage->push_back(he.next().tipVertex());
      if (!e.isBoundary()) {
        storage->push_back(he.twin().next().tipVertex());
      }
    }
    return e.isBoundary() ? 3 : 4;
  }
  case SurfacePointType::Face: {
    if (storage) {
      for (Vertex v : p.face.adjacentVertices()) {
        storage->push_back(v);
      }
    }
    return 3;
  }
  }

  GC_SAFETY_ASSERT(false, "this should not be reachable");
  return 0;
}

std::pair<double, double> compute_local_coordinates(IntrinsicGeometryInterface& geom, Edge e,
                                                    const SurfacePoint& point) {
  SurfacePoint eDummyPt(e, 0.5);
  Face fShared = sharedFace(point, eDummyPt);

  GC_SAFETY_ASSERT(fShared != Face(), "compute_local_coordinates() err: Point " + std::to_string(point) +
                                          " not adjacent to " + std::to_string(e));

  Halfedge he = e.halfedge();

  SurfacePoint pointInFace = point.inFace(fShared);
  SurfacePoint tailInFace = SurfacePoint(he, 0).inFace(fShared);
  SurfacePoint tipInFace = SurfacePoint(he, 1).inFace(fShared);

  // lengths[i] is the length of the edge opposite the i'th vertex
  Vector3 triangleLengths{geom.edgeLengths[fShared.halfedge().next().edge()],
                          geom.edgeLengths[fShared.halfedge().next().next().edge()],
                          geom.edgeLengths[fShared.halfedge().edge()]};
  Vector3 tailDisp = pointInFace.faceCoords - tailInFace.faceCoords;
  double d0 = displacementLength(tailDisp, triangleLengths);

  if (d0 < 1e-50) {
    return std::make_pair(0, 0);
  }

  Vector3 tipDisp = pointInFace.faceCoords - tipInFace.faceCoords;
  double d1 = displacementLength(tipDisp, triangleLengths);
  if (d1 < 1e-50) {
    return std::make_pair(geom.edgeLengths[e], 0);
  }

  double x = geom.edgeLengths[e] / 2.0 + (d0 * d0 - d1 * d1) / (2.0 * geom.edgeLengths[e]);
  double y = sqrt(std::max(0.0, d0 * d0 - x * x));

  return std::make_pair(x, y);
}

Vector2 transformToCoordinateSystem(IntrinsicGeometryInterface& geom, Vector2 v, Face f, SurfacePoint p) {
  Vector2 vTransformed = v;

  geom.requireHalfedgeVectorsInFace();
  geom.requireHalfedgeVectorsInVertex();
  // Transform to coordinate system at point
  switch (p.type) {
  case SurfacePointType::Vertex: {
    // find shared halfedge
    Halfedge heShared;
    for (Halfedge h : f.adjacentHalfedges()) {
      if (h.vertex() == p.vertex) {
        heShared = h;
        break;
      }
    }
    GC_SAFETY_ASSERT(heShared != Halfedge(), "failed to find shared halfedge");

    // Transform from face to halfedge
    vTransformed /= geom.halfedgeVectorsInFace[heShared];

    // Transform from halfedge to vertex
    vTransformed *= geom.halfedgeVectorsInVertex[heShared];
    break;
  }
  case SurfacePointType::Edge: {
    // find shared halfedge
    Halfedge heShared;
    for (Halfedge h : p.edge.adjacentHalfedges()) {
      if (h.face() == f) {
        heShared = h;
        break;
      }
    }

    GC_SAFETY_ASSERT(heShared != Halfedge(), "failed to find shared halfedge");

    // Transform from face to halfedge
    vTransformed /= geom.halfedgeVectorsInFace[heShared];

    // Transform from halfedge to edge
    if (!heShared.orientation()) {
      vTransformed *= -1;
    }

    break;
  }
  case SurfacePointType::Face:
    // do nothing; already using face coordinate system
    break;
  }
  geom.unrequireHalfedgeVectorsInVertex();
  geom.unrequireHalfedgeVectorsInFace();

  return vTransformed;
}

} // namespace exactgeodesic

} // namespace surface
} // namespace geometrycentral
