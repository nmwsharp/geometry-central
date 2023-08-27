// Copyright (C) 2008 Danil Kirsanov, MIT License
// (Modified to work in geometry-central. Original code can be found here: https://code.google.com/p/geodesic/)

#include "geometrycentral/surface/exact_geodesics.h"

namespace geometrycentral {
namespace surface {

VertexData<double> exactGeodesicDistance(SurfaceMesh& mesh, IntrinsicGeometryInterface& geom, Vertex v) {
  GeodesicAlgorithmExact mmp(mesh, geom);
  mmp.propagate({SurfacePoint(v)});
  return mmp.getDistanceFunction();
}

GeodesicAlgorithmExact::GeodesicAlgorithmExact(SurfaceMesh& mesh_, IntrinsicGeometryInterface& geom_)
    : m_max_propagation_distance(1e100), mesh(mesh_), geom(geom_), m_memory_allocator(mesh_.nEdges(), mesh_.nEdges()) {

  geom.requireEdgeLengths();

  m_edge_interval_lists = EdgeData<IntervalList>(mesh);
  for (Edge e : mesh.edges()) {
    m_edge_interval_lists[e].initialize(e);
  }

  // Cache vertex manifold status so we don't have to repeatedly check vertices
  vertexIsManifold = mesh.getVertexManifoldStatus();
  vertexIsBoundary = mesh.getVertexBoundaryStatus();
};

// == Adapters for various input types
void GeodesicAlgorithmExact::propagate(const std::vector<Vertex>& sources, double max_propagation_distance,
                                       const std::vector<Vertex>& stop_points) {
  // Call general version
  std::vector<SurfacePoint> source_surface_points, stop_surface_points;
  for (const Vertex v : sources) source_surface_points.push_back(SurfacePoint(v));
  for (const Vertex v : stop_points) stop_surface_points.push_back(SurfacePoint(v));
  propagate(source_surface_points, max_propagation_distance, stop_surface_points);
}
void GeodesicAlgorithmExact::propagate(const SurfacePoint& source, double max_propagation_distance,
                                       const std::vector<SurfacePoint>& stop_points) {
  // Call general version
  return propagate(std::vector<SurfacePoint>{source}, max_propagation_distance, stop_points);
}
void GeodesicAlgorithmExact::propagate(const Vertex& source, double max_propagation_distance,
                                       const std::vector<Vertex>& stop_points) {
  // Call general version
  std::vector<SurfacePoint> source_surface_points{SurfacePoint(source)};
  std::vector<SurfacePoint> stop_surface_points;
  for (const Vertex v : stop_points) stop_surface_points.push_back(SurfacePoint(v));
  propagate(source_surface_points, max_propagation_distance, stop_surface_points);
}

std::vector<SurfacePoint> GeodesicAlgorithmExact::traceBack(const Vertex& point) const {
  // Call general version
  double ignore;
  return traceBack(SurfacePoint(point), ignore);
}

std::vector<SurfacePoint> GeodesicAlgorithmExact::traceBack(const Vertex& point, double& pathLength) const {
  // Call general version
  return traceBack(SurfacePoint(point), pathLength);
}

std::vector<SurfacePoint> GeodesicAlgorithmExact::traceBack(const SurfacePoint& point) const {
  // Call general version
  double ignore;
  return traceBack(point, ignore);
}

std::pair<unsigned, double> GeodesicAlgorithmExact::closestSource(const Vertex& v) const {
  // Call general version
  return closestSource(SurfacePoint(v));
}

double GeodesicAlgorithmExact::getDistance(const Vertex& v) const {
  // Call general version
  return getDistance(SurfacePoint(v));
}

Vector2 GeodesicAlgorithmExact::getDistanceGradient(const Vertex& v) const {
  // Call general version
  return getDistanceGradient(SurfacePoint(v));
}

Vector2 GeodesicAlgorithmExact::getLog(const Vertex& v) const {
  // Call general version
  return getLog(SurfacePoint(v));
}

void GeodesicAlgorithmExact::best_point_on_the_edge_set(const SurfacePoint& point, std::vector<Edge> const& storage,
                                                        const_interval_pointer& best_interval,
                                                        double& best_total_distance, double& best_interval_position,
                                                        bool verbose) const {

  best_total_distance = 1e100;
  for (Edge e : storage) {
    if (verbose) {
      std::cout << "\tconsidering edge " << e << std::endl;
    }
    const_list_pointer list = interval_list(e);

    double offset;
    double distance;
    const_interval_pointer interval;

    list->find_closest_point(geom, point, offset, distance, interval, verbose);

    if (distance < best_total_distance) {
      best_interval = interval;
      best_total_distance = distance;
      best_interval_position = offset;
    }
  }
}

void GeodesicAlgorithmExact::possible_traceback_edges(const SurfacePoint& point, std::vector<Edge>& storage) const {
  storage.clear();

  switch (point.type) {
  case SurfacePointType::Vertex: {
    Vertex v = point.vertex;
    for (Halfedge he : v.outgoingHalfedges()) {
      if (he.face().isBoundaryLoop()) continue; // don't try to trace along virtual boundary face
      storage.push_back(he.next().edge());
    }
    break;
  }
  case SurfacePointType::Edge: {
    for (Halfedge he : point.edge.adjacentHalfedges()) {
      if (he.isInterior()) {
        storage.push_back(he.next().edge());
        storage.push_back(he.next().next().edge());
      }
    }
    break;
  }
  case SurfacePointType::Face: {
    Face f = point.face;
    for (Edge e : f.adjacentEdges()) {
      storage.push_back(e);
    }
    break;
  }
  }
}


// negative if not visible
long GeodesicAlgorithmExact::visible_from_source(const SurfacePoint& point) const {
  // TODO: interesting branches still untested
  switch (point.type) {
  case SurfacePointType::Vertex: {
    Vertex v = point.vertex;
    for (Halfedge he : v.outgoingHalfedges()) {
      Edge e = he.edge();
      const_list_pointer list = interval_list(e);

      double edgeLength = geom.edgeLengths[e];
      double position = he.orientation() ? 0.0 : edgeLength;

      const_interval_pointer interval = list->covering_interval(position);
      if (interval && interval->visible_from_source()) {
        return (long)interval->source_index();
      }
    }
    return -1;
  }
  case SurfacePointType::Edge: {
    Edge e = point.edge;
    const_list_pointer list = interval_list(e);

    Vertex gcTail = e.halfedge().tailVertex();

    // TODO: extrinsic distance
    double edgeLength = geom.edgeLengths[e];
    double position = std::min(exactgeodesic::compute_surface_distance(geom, point, gcTail), edgeLength);

    const_interval_pointer interval = list->covering_interval(position);
    // assert(interval);
    if (interval && interval->visible_from_source()) {
      return (long)interval->source_index();
    } else {
      return -1;
    }
  }
  case SurfacePointType::Face: {
    return -1;
  }
  }

  GC_SAFETY_ASSERT(false, "This should be unreachable");
  return 0;
}

double GeodesicAlgorithmExact::compute_positive_intersection(double start, double pseudo_x, double pseudo_y,
                                                             double sin_alpha, double cos_alpha) {
  assert(pseudo_y < 0);

  double denominator = sin_alpha * (pseudo_x - start) - cos_alpha * pseudo_y;
  if (denominator < 0.0) {
    return -1.0;
  }

  double numerator = -pseudo_y * start;

  if (numerator < 1e-30) {
    return 0.0;
  }

  if (denominator < 1e-30) {
    return -1.0;
  }

  return numerator / denominator;
}

void GeodesicAlgorithmExact::list_edges_visible_from_source(const SurfacePoint& source,
                                                            std::vector<Edge>& storage) const {

  switch (source.type) {
  case SurfacePointType::Vertex: {
    Vertex v = source.vertex;
    for (Edge e : v.adjacentEdges()) {
      storage.push_back(e);
    }
    return;
  }
  case SurfacePointType::Edge: {
    Edge e = source.edge;
    storage.push_back(e);
    return;
  }
  case SurfacePointType::Face: {
    Face f = source.face;
    for (Edge e : f.adjacentEdges()) {
      storage.push_back(e);
    }
    return;
  }
  }
  GC_SAFETY_ASSERT(false, "This should be unreachable");
}

bool GeodesicAlgorithmExact::erase_from_queue(interval_pointer p) {
  if (p->min() < GEODESIC_INF / 10.0) // && p->min >= queue->begin()->first)
  {
    assert(m_queue.count(p) <= 1); // the set is unique

    IntervalQueue::iterator it = m_queue.find(p);

    if (it != m_queue.end()) {
      m_queue.erase(it);
      return true;
    }
  }

  return false;
}

// intersecting two intervals with up to three intervals in the end
unsigned GeodesicAlgorithmExact::intersect_intervals(interval_pointer zero, IntervalWithStop* one) {
  assert(zero->edge() == one->edge());
  assert(zero->stop() > one->start() && zero->start() < one->stop());
  assert(one->min() < GEODESIC_INF / 10.0);

  double const local_epsilon = SMALLEST_INTERVAL_RATIO * one->edge_length();

  unsigned N = 0;
  if (zero->min() > GEODESIC_INF / 10.0) {
    start[0] = zero->start();
    if (zero->start() < one->start() - local_epsilon) {
      map[0] = MapType::OLD;
      start[1] = one->start();
      map[1] = MapType::NEW;
      N = 2;
    } else {
      map[0] = MapType::NEW;
      N = 1;
    }

    if (zero->stop() > one->stop() + local_epsilon) {
      map[N] = MapType::OLD; // "zero" interval
      start[N++] = one->stop();
    }

    start[N + 1] = zero->stop();
    return N;
  }

  double const local_small_epsilon = 1e-8 * one->edge_length();

  double D = zero->d() - one->d();
  double x0 = zero->pseudo_x();
  double x1 = one->pseudo_x();
  double R0 = x0 * x0 + zero->pseudo_y() * zero->pseudo_y();
  double R1 = x1 * x1 + one->pseudo_y() * one->pseudo_y();

  double inter[2]; // points of intersection
  int Ninter = 0;  // number of the points of the intersection

  if (std::abs(D) < local_epsilon) { // if d1 == d0, equation is linear
    double denom = x1 - x0;
    if (std::abs(denom) > local_small_epsilon) {
      inter[0] = (R1 - R0) / (2. * denom); // one solution
      Ninter = 1;
    }
  } else {
    double D2 = D * D;
    double Q = 0.5 * (R1 - R0 - D2);
    double X = x0 - x1;

    double A = X * X - D2;
    double B = Q * X + D2 * x0;
    double C = Q * Q - D2 * R0;

    if (std::abs(A) < local_small_epsilon) { // if A == 0, linear equation
      if (std::abs(B) > local_small_epsilon) {
        inter[0] = -C / B; // one solution
        Ninter = 1;
      }
    } else {
      double det = B * B - A * C;
      if (det > local_small_epsilon * local_small_epsilon) { // two roots
        det = sqrt(det);
        if (A > 0.0) // make sure that the roots are ordered
        {
          inter[0] = (-B - det) / A;
          inter[1] = (-B + det) / A;
        } else {
          inter[0] = (-B + det) / A;
          inter[1] = (-B - det) / A;
        }

        if (inter[1] - inter[0] > local_small_epsilon) {
          Ninter = 2;
        } else {
          Ninter = 1;
        }
      } else if (det >= 0.0) { // single root
        inter[0] = -B / A;
        Ninter = 1;
      }
    }
  }

  //-----------------------find possible intervals-------------------------------
  // define left and right boundaries of the intersection of the intervals
  double left = std::max(zero->start(), one->start());
  double right = std::min(zero->stop(), one->stop());

  // points of intersection within the (left, right) limits + "left" + "right"
  double good_start[4];
  good_start[0] = left;
  int Ngood_start = 1; // number of the points of the intersection

  for (int i = 0; i < Ninter; ++i) { // for all points of intersection
    double x = inter[i];
    if (x > left + local_epsilon && x < right - local_epsilon) {
      good_start[Ngood_start++] = x;
    }
  }
  good_start[Ngood_start++] = right;

  MapType mid_map[3];
  for (int i = 0; i < Ngood_start - 1; ++i) {
    double mid = (good_start[i] + good_start[i + 1]) * 0.5;
    mid_map[i] = zero->signal(mid) <= one->signal(mid) ? MapType::OLD : MapType::NEW;
  }

  //-----------------------------------output----------------------------------
  N = 0;
  if (zero->start() < left - local_epsilon) { // additional "zero" interval

    if (mid_map[0] == MapType::OLD) { // first interval in the map is already the old one
      good_start[0] = zero->start();
    } else {
      map[N] = MapType::OLD; //"zero" interval
      start[N++] = zero->start();
    }
  }

  for (long i = 0; i < Ngood_start - 1; ++i) { // for all intervals
    MapType current_map = mid_map[i];
    if (N == 0 || map[N - 1] != current_map) {
      map[N] = current_map;
      start[N++] = good_start[i];
    }
  }

  if (zero->stop() > one->stop() + local_epsilon) {
    if (N == 0 || map[N - 1] == MapType::NEW) {
      map[N] = MapType::OLD; //"zero" interval
      start[N++] = one->stop();
    }
  }

  start[0] = zero->start(); // just to make sure that epsilons do not damage anything
  // start[N] = zero->stop();

  return N;
}

void GeodesicAlgorithmExact::initialize_propagation_data() {
  clear();

  IntervalWithStop candidate;
  std::vector<Edge> edges_visible_from_source;
  for (unsigned i = 0; i < m_sources.size(); ++i) // for all edges adjacent to the starting vertex
  {
    SurfacePoint* source = &m_sources[i];

    edges_visible_from_source.clear();
    list_edges_visible_from_source(*source, edges_visible_from_source);

    for (Edge e : edges_visible_from_source) {
      double edge_length = geom.edgeLengths[e];
      candidate.initialize(geom, e, edge_length, source, i);
      candidate.stop() = edge_length;
      candidate.compute_min_distance(candidate.stop());
      candidate.direction() = Interval::DirectionType::FROM_SOURCE;

      update_list_and_queue(interval_list(e), &candidate, 1);
    }
  }
}

static int reps = 0;

// propagation algorithm stops after reaching the certain distance from the
// source
void GeodesicAlgorithmExact::propagate(const std::vector<SurfacePoint>& sources, double max_propagation_distance,
                                       const std::vector<SurfacePoint>& stop_points) {

  set_stop_conditions(stop_points, max_propagation_distance);
  set_sources(sources);
  initialize_propagation_data();

  clock_t start = clock();

  unsigned satisfied_index = 0;

  m_iterations = 0; // for statistics
  m_queue_max_size = 0;

  IntervalWithStop candidates[2];
  geom.requireCornerAngles();

  while (!m_queue.empty()) {
    m_queue_max_size = std::max(m_queue.size(), m_queue_max_size);

    unsigned const check_period = 10;
    if (++m_iterations % check_period == 0) // check if we covered all required vertices
    {
      if (check_stop_conditions(satisfied_index)) {
        break;
      }
    }

    interval_pointer min_interval = *m_queue.begin();
    m_queue.erase(m_queue.begin());
    Edge edge = min_interval->edge();
    Halfedge he = edge.halfedge();
    list_pointer list = interval_list(edge);

    assert(min_interval->d() < GEODESIC_INF);

    bool const first_interval = min_interval->start() == 0.0;
    // bool const last_interval = min_interval->stop() == edge->length();
    bool const last_interval = min_interval->next() == nullptr;

    auto saddleOrBoundary = [&](Vertex v) -> bool {
      geom.requireVertexGaussianCurvatures();
      bool saddle = geom.vertexGaussianCurvatures[v] < 0;
      geom.unrequireVertexGaussianCurvatures();
      return saddle || vertexIsBoundary[v];
    };

    bool const turn_left = saddleOrBoundary(he.tailVertex());
    bool const turn_right = saddleOrBoundary(he.tipVertex());

    auto propagateNonmanifoldVertex = [&](Vertex v, bool isTailVertex) {
      IntervalWithStop candidate;
      double Le = geom.edgeLengths[min_interval->edge()];
      double dx = isTailVertex ? min_interval->pseudo_x() : Le - min_interval->pseudo_x();
      double dy = min_interval->pseudo_y();
      double source_d = min_interval->d() + sqrt(dx * dx + dy * dy);

      for (Halfedge heSpoke : v.outgoingHalfedges()) {
        double L = geom.edgeLengths[heSpoke.edge()];
        candidate.start() = 0;
        candidate.stop() = L;
        candidate.d() = source_d;
        candidate.pseudo_x() = heSpoke.orientation() ? 0 : L;
        candidate.pseudo_y() = 0;

        candidate.next() = nullptr;
        candidate.edge() = heSpoke.edge();
        candidate.halfedge() = heSpoke;
        candidate.edge_length() = L;
        candidate.source_index() = min_interval->source_index();
        candidate.direction() = Interval::DirectionType::FROM_HALFEDGE;
        candidate.min() = candidate.d();

        update_list_and_queue(interval_list(heSpoke.edge()), &candidate, 1);
      }
    };

    if (!vertexIsManifold[he.tailVertex()]) {
      propagateNonmanifoldVertex(he.tailVertex(), true);
    }
    if (!vertexIsManifold[he.tipVertex()]) {
      propagateNonmanifoldVertex(he.tipVertex(), false);
    }

    for (Halfedge neighboring_halfedge : edge.adjacentHalfedges()) {

      //== Check some early exit conditions
      if (!neighboring_halfedge.isInterior()) continue;

      // just in case, always propagate boundary edges
      if (!edge.isBoundary() && min_interval->halfedge() == neighboring_halfedge) continue;

      Face face = neighboring_halfedge.face();

      // if we come from 1, go to 2
      // Edge next_edge = (neighboring_halfedge.orientation() == he.orientation()) ?
      // neighboring_halfedge.next().next().edge() : neighboring_halfedge.next().edge();
      Halfedge next_halfedge = (neighboring_halfedge.orientation() == he.orientation())
                                   ? neighboring_halfedge.next().next()
                                   : neighboring_halfedge.next();
      Edge next_edge = next_halfedge.edge();

      double next_edge_length = geom.edgeLengths[next_edge];
      double corner_angle = (neighboring_halfedge.orientation() == he.orientation())
                                ? geom.cornerAngles[neighboring_halfedge.corner()]
                                : geom.cornerAngles[neighboring_halfedge.next().corner()];

      unsigned num_propagated = compute_propagated_parameters(min_interval->pseudo_x(), min_interval->pseudo_y(),
                                                              min_interval->d(), // parameters of the interval
                                                              min_interval->start(),
                                                              min_interval->stop(), // start/end of the interval
                                                              corner_angle,         // corner angle
                                                              next_edge_length,     // length of the new edge
                                                              first_interval, // if it is the first interval on the edge
                                                              last_interval,  // if it is the last interval on the edge
                                                              turn_left, turn_right, candidates);
      bool propagate_to_right = true;

      if (num_propagated) {
        if (candidates[num_propagated - 1].stop() != next_edge_length) {
          propagate_to_right = false;
        }

        // if the origins coinside, do not invert intervals
        // bool const invert = next_halfedge.edge().halfedge().tailVertex() != he.tailVertex();
        bool const invert = (next_halfedge == neighboring_halfedge.next().next()) ? next_halfedge.orientation()
                                                                                  : !next_halfedge.orientation();

        // invert: do not inverse
        construct_propagated_intervals(invert, next_halfedge, candidates, num_propagated, min_interval);

        update_list_and_queue(interval_list(next_edge), candidates, num_propagated);
      }

      if (propagate_to_right) {
        // propogation to the right edge
        double length = geom.edgeLengths[edge];

        // next_edge = (neighboring_halfedge.orientation() == he.orientation()) ? neighboring_halfedge.next().edge()
        // : neighboring_halfedge.next().next().edge();
        next_halfedge = (neighboring_halfedge.orientation() == he.orientation()) ? neighboring_halfedge.next()
                                                                                 : neighboring_halfedge.next().next();
        next_edge = next_halfedge.edge();

        next_edge_length = geom.edgeLengths[next_edge];
        corner_angle = (neighboring_halfedge.orientation() == he.orientation())
                           ? geom.cornerAngles[neighboring_halfedge.next().corner()]
                           : geom.cornerAngles[neighboring_halfedge.corner()];

        num_propagated = compute_propagated_parameters(length - min_interval->pseudo_x(), min_interval->pseudo_y(),
                                                       min_interval->d(), // parameters of the interval
                                                       length - min_interval->stop(),
                                                       length - min_interval->start(), // start/end of the interval
                                                       corner_angle,                   // corner angle
                                                       next_edge_length,               // length of the new edge
                                                       last_interval, // if it is the first interval on the edge
                                                       first_interval, turn_right, turn_left,
                                                       candidates); // if it is the last interval on the edge

        if (num_propagated) {
          // if the origins coinside, do not invert intervals
          // bool const invert = next_halfedge.edge().halfedge().tailVertex() != he.tipVertex();
          bool const invert = (next_halfedge == neighboring_halfedge.next().next()) ? next_halfedge.orientation()
                                                                                    : !next_halfedge.orientation();

          // invert: do not inverse
          construct_propagated_intervals(invert, next_halfedge, candidates, num_propagated, min_interval);

          update_list_and_queue(interval_list(next_edge), candidates, num_propagated);
        }
      }
    }
  }

  m_propagation_distance_stopped = m_queue.empty() ? GEODESIC_INF : (*m_queue.begin())->min();
  clock_t stop = clock();
  m_time_consumed = (static_cast<double>(stop) - static_cast<double>(start)) / CLOCKS_PER_SEC;

  geom.unrequireCornerAngles();

  /*	for(unsigned i=0; i<m_edge_interval_lists.size(); ++i)
    {
      list_pointer list = &m_edge_interval_lists[i];
      interval_pointer p = list->first();
      assert(p->start() == 0.0);
      while(p->next())
      {
        assert(p->stop() == p->next()->start());
        assert(p->d() < GEODESIC_INF);
        p = p->next();
      }
    }*/
}


bool GeodesicAlgorithmExact::check_stop_conditions(unsigned& index) const {
  double queue_distance = (*m_queue.begin())->min();
  double stop_dist = stop_distance();
  if (stop_dist < GEODESIC_INF && queue_distance < stop_distance()) {
    return false;
  }

  if (m_stop_vertices.empty()) {
    return false;
  }

  while (index < m_stop_vertices.size()) {
    Vertex v = m_stop_vertices[index].first;
    Edge edge = v.halfedge().edge(); // take any edge

    double distance = v.halfedge().orientation() ? interval_list(edge)->signal(0.0)
                                                 : interval_list(edge)->signal(geom.edgeLengths[edge]);

    if (queue_distance < distance + m_stop_vertices[index].second) {
      return false;
    }

    ++index;
  }
  return true;
}


void GeodesicAlgorithmExact::update_list_and_queue(list_pointer list,
                                                   IntervalWithStop* candidates, // up to two candidates
                                                   unsigned num_candidates) {
  assert(num_candidates <= 2);
  // assert(list->first() != nullptr);
  Edge edge = list->edge();
  double edge_length = geom.edgeLengths[edge];
  double const local_epsilon = SMALLEST_INTERVAL_RATIO * edge_length;

  if (list->first() == nullptr) {
    interval_pointer* p = &list->first();
    IntervalWithStop* first;
    IntervalWithStop* second;

    if (num_candidates == 1) {
      first = candidates;
      second = candidates;
      first->compute_min_distance(first->stop());
    } else {
      if (candidates->start() <= (candidates + 1)->start()) {
        first = candidates;
        second = candidates + 1;
      } else {
        first = candidates + 1;
        second = candidates;
      }
      assert(first->stop() == second->start());

      first->compute_min_distance(first->stop());
      second->compute_min_distance(second->stop());
    }

    if (first->start() > 0.0) {
      *p = m_memory_allocator.allocate();
      (*p)->initialize(geom, edge, edge_length);
      p = &(*p)->next();
    }

    *p = m_memory_allocator.allocate();
    // memcpy(*p, first, sizeof(Interval));
    **p = *first;
    m_queue.insert(*p);

    if (num_candidates == 2) {
      p = &(*p)->next();
      *p = m_memory_allocator.allocate();
      // memcpy(*p, second, sizeof(Interval));
      **p = *second;
      m_queue.insert(*p);
    }

    if (second->stop() < edge_length) {
      p = &(*p)->next();
      *p = m_memory_allocator.allocate();
      (*p)->initialize(geom, edge, edge_length);
      (*p)->start() = second->stop();
    } else {
      (*p)->next() = nullptr;
    }
    return;
  }

  bool propagate_flag;

  // for all new intervals
  for (unsigned i = 0; i < num_candidates; ++i) {
    IntervalWithStop* q = &candidates[i];

    interval_pointer previous = nullptr;

    interval_pointer p = list->first();
    assert(p->start() == 0.0);

    while (p != nullptr && p->stop() - local_epsilon < q->start()) {
      p = p->next();
    }

    // go through all old intervals
    while (p != nullptr && p->start() < q->stop() - local_epsilon) {
      unsigned const N = intersect_intervals(p, q); // interset two intervals

      if (N == 1) {
        // if "p" is always better, we do not need to update anything
        if (map[0] == MapType::OLD) {
          if (previous) { // close previous interval and put in into the queue
            previous->next() = p;
            previous->compute_min_distance(p->start());
            m_queue.insert(previous);
            previous = nullptr;
          }

          p = p->next();

        } else if (previous) { // extend previous interval to cover everything; remove p
          previous->next() = p->next();
          erase_from_queue(p);
          m_memory_allocator.deallocate(p);

          p = previous->next();
        } else { // p becomes "previous"
          previous = p;
          interval_pointer next = p->next();
          erase_from_queue(p);

          // memcpy(previous, q, sizeof(Interval));
          *previous = *q;

          previous->start() = start[0];
          previous->next() = next;

          p = next;
        }
        continue;
      }

      // update_flag = true;

      Interval swap(*p); // used for swapping information
      propagate_flag = erase_from_queue(p);

      for (unsigned j = 1; j < N; ++j) { // no memory is needed for the first one
        // create new intervals
        i_new[j] = m_memory_allocator.allocate();
      }

      if (map[0] == MapType::OLD) { // finish previous, if any
        if (previous) {
          previous->next() = p;
          previous->compute_min_distance(previous->stop());
          m_queue.insert(previous);
          previous = nullptr;
        }
        i_new[0] = p;
        p->next() = i_new[1];
        p->start() = start[0];
      } else if (previous) { // extend previous interval to cover everything; remove p
        i_new[0] = previous;
        previous->next() = i_new[1];
        m_memory_allocator.deallocate(p);
        previous = nullptr;
      } else { // p becomes "previous"
        i_new[0] = p;
        // memcpy(p, q, sizeof(Interval));
        *p = *q;

        p->next() = i_new[1];
        p->start() = start[0];
      }

      assert(!previous);

      for (unsigned j = 1; j < N; ++j) {
        interval_pointer current_interval = i_new[j];

        if (map[j] == MapType::OLD) {
          // memcpy(current_interval, &swap, sizeof(Interval));
          *current_interval = swap;
        } else {
          // memcpy(current_interval, q, sizeof(Interval));
          *current_interval = *q;
        }

        if (j == N - 1) {
          current_interval->next() = swap.next();
        } else {
          current_interval->next() = i_new[j + 1];
        }

        current_interval->start() = start[j];
      }

      for (unsigned j = 0; j < N; ++j) { // find "min" and add the intervals to the queue
        if (j == N - 1 && map[j] == MapType::NEW) {
          previous = i_new[j];
        } else {
          interval_pointer current_interval = i_new[j];

          current_interval->compute_min_distance(current_interval->stop()); // compute minimal distance

          if (map[j] == MapType::NEW || (map[j] == MapType::OLD && propagate_flag)) {
            m_queue.insert(current_interval);
          }
        }
      }

      p = swap.next();
    }

    if (previous) { // close previous interval and put in into the queue
      previous->compute_min_distance(previous->stop());
      m_queue.insert(previous);
      previous = nullptr;
    }
  }
}

unsigned
GeodesicAlgorithmExact::compute_propagated_parameters(double pseudo_x, double pseudo_y,
                                                      double d, // parameters of the interval
                                                      double begin,
                                                      double end,          // start/end of the interval
                                                      double alpha,        // corner angle
                                                      double L,            // length of the new edge
                                                      bool first_interval, // if it is the first interval on the edge
                                                      bool last_interval,  // if it is the last interval on the edge
                                                      bool turn_left, bool turn_right, IntervalWithStop* candidates) {
  assert(pseudo_y <= 0.0);
  assert(d < GEODESIC_INF / 10.0);
  assert(begin <= end);
  assert(first_interval ? (begin == 0.0) : true);

  IntervalWithStop* p = candidates;

  // pseudo-source is on the edge
  if (std::abs(pseudo_y) <= 1e-30) {
    if (first_interval && pseudo_x <= 0.0) {
      p->start() = 0.0;
      p->stop() = L;
      p->d() = d - pseudo_x;
      p->pseudo_x() = 0.0;
      p->pseudo_y() = 0.0;
      return 1;
    } else if (last_interval && pseudo_x >= end) {
      p->start() = 0.0;
      p->stop() = L;
      p->d() = d + pseudo_x - end;
      p->pseudo_x() = end * cos(alpha);
      p->pseudo_y() = -end * sin(alpha);
      return 1;
    } else if (pseudo_x >= begin && pseudo_x <= end) {
      p->start() = 0.0;
      p->stop() = L;
      p->d() = d;
      p->pseudo_x() = pseudo_x * cos(alpha);
      p->pseudo_y() = -pseudo_x * sin(alpha);
      return 1;
    } else {
      return 0;
    }
  }

  double sin_alpha = sin(alpha);
  double cos_alpha = cos(alpha);

  // important: for the first_interval, this function returns zero only if the
  // new edge is "visible" from the source.
  // if the new edge can be covered only after turn_over, the value is negative (-1.0)
  double L1 = compute_positive_intersection(begin, pseudo_x, pseudo_y, sin_alpha, cos_alpha);

  if (L1 < 0 || L1 >= L) {
    if (first_interval && turn_left) {
      p->start() = 0.0;
      p->stop() = L;
      p->d() = d + sqrt(pseudo_x * pseudo_x + pseudo_y * pseudo_y);
      p->pseudo_y() = 0.0;
      p->pseudo_x() = 0.0;
      return 1;
    } else {
      return 0;
    }
  }

  double L2 = compute_positive_intersection(end, pseudo_x, pseudo_y, sin_alpha, cos_alpha);

  if (L2 < 0 || L2 >= L) {
    p->start() = L1;
    p->stop() = L;
    p->d() = d;
    p->pseudo_x() = cos_alpha * pseudo_x + sin_alpha * pseudo_y;
    p->pseudo_y() = -sin_alpha * pseudo_x + cos_alpha * pseudo_y;

    return 1;
  }

  p->start() = L1;
  p->stop() = L2;
  p->d() = d;
  p->pseudo_x() = cos_alpha * pseudo_x + sin_alpha * pseudo_y;
  p->pseudo_y() = -sin_alpha * pseudo_x + cos_alpha * pseudo_y;
  assert(p->pseudo_y() <= 0.0);

  if (!(last_interval && turn_right)) {
    return 1;
  } else {
    p = candidates + 1;

    p->start() = L2;
    p->stop() = L;
    double dx = pseudo_x - end;
    p->d() = d + sqrt(dx * dx + pseudo_y * pseudo_y);
    p->pseudo_x() = end * cos_alpha;
    p->pseudo_y() = -end * sin_alpha;

    return 2;
  }
}

void GeodesicAlgorithmExact::construct_propagated_intervals(bool invert, Halfedge halfedge,
                                                            IntervalWithStop* candidates,
                                                            unsigned& num_candidates, // up to two candidates
                                                            interval_pointer source_interval) {
  // constructs iNew from the rest of the data

  Edge edge = halfedge.edge();

  double edge_length = geom.edgeLengths[edge];
  double local_epsilon = SMALLEST_INTERVAL_RATIO * edge_length;

  // kill very small intervals in order to avoid precision problems
  if (num_candidates == 2) {
    double start = std::min(candidates->start(), (candidates + 1)->start());
    double stop = std::max(candidates->stop(), (candidates + 1)->stop());
    if (candidates->stop() - candidates->start() < local_epsilon) // kill interval 0
    {
      *candidates = *(candidates + 1);
      num_candidates = 1;
      candidates->start() = start;
      candidates->stop() = stop;
    } else if ((candidates + 1)->stop() - (candidates + 1)->start() < local_epsilon) {
      num_candidates = 1;
      candidates->start() = start;
      candidates->stop() = stop;
    }
  }

  IntervalWithStop* first;
  IntervalWithStop* second;
  if (num_candidates == 1) {
    first = candidates;
    second = candidates;
  } else {
    if (candidates->start() <= (candidates + 1)->start()) {
      first = candidates;
      second = candidates + 1;
    } else {
      first = candidates + 1;
      second = candidates;
    }
    assert(first->stop() == second->start());
  }

  if (first->start() < local_epsilon) {
    first->start() = 0.0;
  }
  if (edge_length - second->stop() < local_epsilon) {
    second->stop() = edge_length;
  }


  // invert intervals if necessary; fill missing data and set pointers
  // correctly
  Interval::DirectionType direction = Interval::DirectionType::FROM_HALFEDGE;

  if (!invert) // in this case everything is straighforward, we do not have to
               // invert the intervals
  {
    for (unsigned i = 0; i < num_candidates; ++i) {
      IntervalWithStop* p = candidates + i;

      p->next() = (i == num_candidates - 1) ? nullptr : candidates + i + 1;
      p->edge() = edge;
      p->halfedge() = halfedge;
      p->edge_length() = edge_length;
      p->direction() = direction;
      p->source_index() = source_interval->source_index();

      p->min() = 0.0; // it will be changed later on

      assert(p->start() < p->stop());
    }
  } else // now we have to invert the intervals
  {
    for (unsigned i = 0; i < num_candidates; ++i) {
      IntervalWithStop* p = candidates + i;

      p->next() = (i == 0) ? nullptr : candidates + i - 1;
      p->edge() = edge;
      p->halfedge() = halfedge;
      p->edge_length() = edge_length;
      p->direction() = direction;
      p->source_index() = source_interval->source_index();

      double length = edge_length;
      p->pseudo_x() = length - p->pseudo_x();

      double start = length - p->stop();
      p->stop() = length - p->start();
      p->start() = start;

      p->min() = 0;

      assert(p->start() < p->stop());
      assert(p->start() >= 0.0);
      assert(p->stop() <= edge_length);
    }
  }
}


std::pair<unsigned, double> GeodesicAlgorithmExact::closestSource(const SurfacePoint& point) const {

  double best_interval_position;
  unsigned best_source_index;
  double best_source_distance;

  best_first_interval(point, best_source_distance, best_interval_position, best_source_index);

  return std::make_pair(best_source_index, best_source_distance);
}

double GeodesicAlgorithmExact::getDistance(const SurfacePoint& point) const { return closestSource(point).second; }

VertexData<double> GeodesicAlgorithmExact::getDistanceFunction() const {
  VertexData<double> distances(mesh);
  for (Vertex v : mesh.vertices()) {
    distances[v] = closestSource(SurfacePoint(v)).second;
  }
  return distances;
}

Vector2 GeodesicAlgorithmExact::getDistanceGradient(const SurfacePoint& point) const {

  SurfacePoint closestPoint;
  long visibleSource = visible_from_source(point);
  if (visibleSource >= 0) {
    // closest source shares face with point
    // This only handles edge or vertex points which share a face with the source. Face points are handled below
    closestPoint = static_cast<SurfacePoint>(m_sources[visibleSource]);
  } else {
    double best_interval_position;
    unsigned best_source_index;
    double best_source_distance;
    const_interval_pointer best_interval;

    if (point.type == SurfacePointType::Face) {
      best_interval = best_first_interval(point, best_source_distance, best_interval_position, best_source_index);
    } else {
      std::vector<Edge> possible_edges;
      possible_edges.reserve(10);
      possible_traceback_edges(point, possible_edges);

      best_point_on_the_edge_set(point, possible_edges, best_interval, best_source_distance, best_interval_position);
    }

    if (best_interval) {
      // Identity point in interval which lies along the geodesic
      Edge edge = best_interval->edge();
      Halfedge he = edge.halfedge();
      double edge_length = geom.edgeLengths[edge];
      double local_epsilon = SMALLEST_INTERVAL_RATIO * edge_length;

      if (best_interval_position < local_epsilon) {
        closestPoint = SurfacePoint(he.tailVertex());
      } else if (best_interval_position > edge_length - local_epsilon) {
        closestPoint = SurfacePoint(he.tipVertex());
      } else {
        double normalized_position = best_interval_position / edge_length;
        closestPoint = SurfacePoint(edge, normalized_position);
      }
    } else {
      // If this happens, then point should be in the same face as some source
      closestPoint = static_cast<SurfacePoint>(m_sources[best_source_index]);
    }
  }

  // gradient of distance points from closestPoint to point
  Face fShared = sharedFace(point, closestPoint);
  Vector3 baryPoint = point.inFace(fShared).faceCoords;
  Vector3 baryClosestPoint = closestPoint.inFace(fShared).faceCoords;

  auto vertexCoordinatesInTriangle = [](IntrinsicGeometryInterface& geom, Face face) -> std::array<Vector2, 3> {
    geom.requireHalfedgeVectorsInFace();
    std::array<Vector2, 3> result = {Vector2{0., 0.}, geom.halfedgeVectorsInFace[face.halfedge()],
                                     -geom.halfedgeVectorsInFace[face.halfedge().next().next()]};
    geom.unrequireHalfedgeVectorsInFace();
    return result;
  };

  // gradient vector in face fShared
  Vector2 faceDisplacement =
      barycentricDisplacementToCartesian(vertexCoordinatesInTriangle(geom, fShared), baryPoint - baryClosestPoint);

  Vector2 displacement = exactgeodesic::transformToCoordinateSystem(geom, faceDisplacement, fShared, point);

  return displacement.normalize();
}

Vector2 GeodesicAlgorithmExact::getLog(const SurfacePoint& point) const {
  double pathLen;
  std::vector<SurfacePoint> path = traceBack(point, pathLen);

  if (path.size() <= 1) { // degenerate path means point is closest source
    return Vector2::zero();
  }

  int N = path.size() - 1;
  SurfacePoint closestSource = path[N];
  Face fSource = sharedFace(closestSource, path[N - 1]);

  // compute displacement between second-to-last point and closest source
  // (in barycentric coordinates)
  Vector3 logDirBary = path[N - 1].inFace(fSource).faceCoords - path[N].inFace(fSource).faceCoords;

  geom.requireHalfedgeVectorsInFace();
  std::array<Vector2, 3> vertCoordsInSourceFace{Vector2{0., 0.}, geom.halfedgeVectorsInFace[fSource.halfedge()],
                                                -geom.halfedgeVectorsInFace[fSource.halfedge().next().next()]};
  geom.unrequireHalfedgeVectorsInFace();

  // convert displacement to cartesian coordinates in fSource
  Vector2 faceLogDir = barycentricDisplacementToCartesian(vertCoordsInSourceFace, logDirBary);

  // transform displacement into coordinate system of closestSource
  Vector2 logDir = exactgeodesic::transformToCoordinateSystem(geom, faceLogDir, fSource, closestSource);

  // set log magnitude and return
  return logDir.normalize() * pathLen;
}

const_interval_pointer GeodesicAlgorithmExact::best_first_interval(const SurfacePoint& point,
                                                                   double& best_total_distance,
                                                                   double& best_interval_position,
                                                                   unsigned& best_source_index) const {

  const_interval_pointer best_interval = nullptr;
  best_total_distance = GEODESIC_INF;

  switch (point.type) {
  case SurfacePointType::Vertex: {
    Vertex v = point.vertex;
    for (Halfedge he : v.outgoingHalfedges()) {
      Edge e = he.edge();
      const_list_pointer list = interval_list(e);

      // double position = e->v0()->id() == v->id() ? 0.0 : e->length();
      double position = he.orientation() ? 0.0 : geom.edgeLengths[e];

      const_interval_pointer interval = list->covering_interval(position);
      if (interval) {
        double distance = interval->signal(position);

        if (distance < best_total_distance) {
          best_interval = interval;
          best_total_distance = distance;
          best_interval_position = position;
          best_source_index = interval->source_index();
        }
      }
    }
    break;
  }
  case SurfacePointType::Edge: {
    // TODO: untested branch
    Edge e = point.edge;
    const_list_pointer list = interval_list(e);

    Vertex gcTail = e.halfedge().tailVertex();

    best_interval_position = exactgeodesic::compute_surface_distance(geom, point, gcTail);
    best_interval = list->covering_interval(best_interval_position);
    if (best_interval) {
      // assert(best_interval && best_interval->d() < GEODESIC_INF);
      best_total_distance = best_interval->signal(best_interval_position);
      best_source_index = best_interval->source_index();
    }
    break;
  }
  case SurfacePointType::Face: {
    Face f = point.face;
    for (Edge e : f.adjacentEdges()) {
      const_list_pointer list = interval_list(e);

      double offset;
      double distance;
      const_interval_pointer interval;

      list->find_closest_point(geom, point, offset, distance, interval);

      if (interval && distance < best_total_distance) {
        best_interval = interval;
        best_total_distance = distance;
        best_interval_position = offset;
        best_source_index = interval->source_index();
      }
    }

    // check for all sources that might be located inside this face
    SortedSources::const_sorted_iterator_pair local_sources = m_sources.sources(point);
    for (SortedSources::const_sorted_iterator it = local_sources.first; it != local_sources.second; ++it) {
      SurfacePointWithIndex* source = *it;
      double distance = exactgeodesic::compute_surface_distance(geom, point, *source);
      if (distance < best_total_distance) {
        best_interval = nullptr;
        best_total_distance = distance;
        best_interval_position = 0.0;
        best_source_index = source->index();
      }
    }
    break;
  }
  }

  if (best_total_distance > m_propagation_distance_stopped) // result is unreliable
  {
    best_total_distance = GEODESIC_INF;
    return nullptr;
  } else {
    return best_interval;
  }
}

// trace back piecewise-linear path
std::vector<SurfacePoint> GeodesicAlgorithmExact::traceBack(const SurfacePoint& point, double& pathLength) const {
  std::vector<SurfacePoint> path;

  double best_interval_position;
  unsigned source_index = std::numeric_limits<unsigned>::max();
  const_interval_pointer best_interval = best_first_interval(point, pathLength, best_interval_position, source_index);

  // unable to find the right path
  if (pathLength >= GEODESIC_INF / 2.0) {
    return {};
  }

  path.push_back(point);

  if (best_interval) // if we did not hit the face source immediately
  {
    std::vector<Edge> possible_edges;
    possible_edges.reserve(10);

    // while this point is not in the direct visibility of some
    // source (if we are inside the FACE, we obviously hit the source)
    while (visible_from_source(path.back()) < 0) {
      SurfacePoint& q = path.back();

      possible_traceback_edges(q, possible_edges);

      const_interval_pointer interval;
      double total_distance;
      double position;

      best_point_on_the_edge_set(q, possible_edges, interval, total_distance, position);
      assert(total_distance < GEODESIC_INF);
      source_index = interval->source_index();


      Edge edge = interval->edge();
      Halfedge he = edge.halfedge();
      double edge_length = geom.edgeLengths[edge];
      double local_epsilon = SMALLEST_INTERVAL_RATIO * edge_length;

      if (position < local_epsilon) {
        path.push_back(SurfacePoint(he.tailVertex()));
      } else if (position > edge_length - local_epsilon) {
        path.push_back(SurfacePoint(he.tipVertex()));
      } else {
        double normalized_position = position / edge_length;
        path.push_back(SurfacePoint(edge, normalized_position));
      }
    }
  }

  SurfacePoint source = static_cast<SurfacePoint>(m_sources[source_index]);
  if (exactgeodesic::compute_surface_distance(geom, path.back(), source) > 0) {
    path.push_back(source);
  }

  return path;
}

void GeodesicAlgorithmExact::print_statistics() const {
  std::cout << "propagation step took " << m_time_consumed << " seconds " << std::endl;

  unsigned interval_counter = 0;
  for (unsigned i = 0; i < m_edge_interval_lists.size(); ++i) {
    interval_counter += m_edge_interval_lists[i].number_of_intervals();
  }
  double intervals_per_edge = (double)interval_counter / (double)m_edge_interval_lists.size();

  double memory = m_edge_interval_lists.size() * sizeof(IntervalList) + interval_counter * sizeof(Interval);

  std::cout << "uses about " << memory / 1e6 << "Mb of memory" << std::endl;
  std::cout << interval_counter << " total intervals, or " << intervals_per_edge << " intervals per edge" << std::endl;
  std::cout << "maximum interval queue size is " << m_queue_max_size << std::endl;
  std::cout << "number of interval propagations is " << m_iterations << std::endl;
}

IntervalList GeodesicAlgorithmExact::getEdgeIntervals(Edge e) const { return m_edge_interval_lists[e]; }

void GeodesicAlgorithmExact::set_stop_conditions(const std::vector<SurfacePoint>& stop_points, double stop_distance) {
  m_max_propagation_distance = stop_distance;

  if (stop_points.empty()) {
    m_stop_vertices.clear();
    return;
  }

  m_stop_vertices.reserve(stop_points.size());

  // TODO: not tested
  std::vector<Vertex> possible_vertices;
  for (const SurfacePoint& point : stop_points) {
    possible_vertices.clear();
    exactgeodesic::compute_closest_vertices(point, &possible_vertices);

    Vertex closest_vertex;
    double min_distance = 1e100;
    for (Vertex v : possible_vertices) {
      double distance = exactgeodesic::compute_surface_distance(geom, point, SurfacePoint(v));
      if (distance < min_distance) {
        min_distance = distance;
        closest_vertex = v;
      }
    }
    assert(closest_vertex != Vertex());

    m_stop_vertices.push_back(std::make_pair(closest_vertex, min_distance));
  }
}

void GeodesicAlgorithmExact::clear() {
  m_memory_allocator.clear();
  m_queue.clear();
  for (unsigned i = 0; i < m_edge_interval_lists.size(); ++i) {
    m_edge_interval_lists[i].clear();
  }
  m_propagation_distance_stopped = GEODESIC_INF;
};

} // namespace surface
} // namespace geometrycentral
