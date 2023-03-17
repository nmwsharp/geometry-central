// Copyright (C) 2008 Danil Kirsanov, MIT License
// (Modified to work in geometry-central. Original code can be found here: https://code.google.com/p/geodesic/)

#pragma once

#include "geometrycentral/surface/exact_geodesic_helpers.h"
#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/surface_point.h"

#include <assert.h>
#include <cmath>
#include <set>
#include <vector>

namespace geometrycentral {
namespace surface {

// One-off function to compute distance from a vertex
VertexData<double> exactGeodesicDistance(SurfaceMesh& mesh, IntrinsicGeometryInterface& geom, Vertex v);

class GeodesicAlgorithmExact {
public:
  GeodesicAlgorithmExact(SurfaceMesh& mesh_, IntrinsicGeometryInterface& geom_);
  ~GeodesicAlgorithmExact(){};

  // propagation algorithm stops after reaching the certain distance from the
  // source or after ensuring that all the stop_points are covered
  void propagate(const std::vector<SurfacePoint>& sources, double max_propagation_distance = GEODESIC_INF,
                 const std::vector<SurfacePoint>& stop_points = {});
  void propagate(const std::vector<Vertex>& sources, double max_propagation_distance = GEODESIC_INF,
                 const std::vector<Vertex>& stop_points = {});
  void propagate(const SurfacePoint& source, double max_propagation_distance = GEODESIC_INF,
                 const std::vector<SurfacePoint>& stop_points = {});
  void propagate(const Vertex& source, double max_propagation_distance = GEODESIC_INF,
                 const std::vector<Vertex>& stop_points = {});

  // trace back piecewise-linear path
  // the resulting path starts at "point" and ends at the closest source
  // optionally also returns the length of the path via the pathLength argument
  std::vector<SurfacePoint> traceBack(const SurfacePoint& point) const;
  std::vector<SurfacePoint> traceBack(const Vertex& point) const;
  std::vector<SurfacePoint> traceBack(const SurfacePoint& point, double& pathLength) const;
  std::vector<SurfacePoint> traceBack(const Vertex& point, double& pathLength) const;

  // quickly find what source this point belongs to and what is the distance
  // to this source
  std::pair<unsigned, double> closestSource(const SurfacePoint& point) const;
  std::pair<unsigned, double> closestSource(const Vertex& point) const;

  // evaluate distance function at a point
  double getDistance(const SurfacePoint& point) const;
  double getDistance(const Vertex& point) const;

  // evaluate gradient of distance function at a point
  Vector2 getDistanceGradient(const SurfacePoint& point) const;
  Vector2 getDistanceGradient(const Vertex& point) const;

  // evaluate log map at closest source
  Vector2 getLog(const SurfacePoint& point) const;
  Vector2 getLog(const Vertex& point) const;

  // evaluate distance function at all vertices
  VertexData<double> getDistanceFunction() const;

  void print_statistics() const;

  IntervalList getEdgeIntervals(Edge e) const;

protected:
  typedef std::set<interval_pointer, Interval> IntervalQueue;

  void update_list_and_queue(list_pointer list,
                             IntervalWithStop* candidates, // up to two candidates
                             unsigned num_candidates);

  unsigned compute_propagated_parameters(double pseudo_x, double pseudo_y,
                                         double d, // parameters of the interval
                                         double start,
                                         double end,          // start/end of the interval
                                         double alpha,        // corner angle
                                         double L,            // length of the new edge
                                         bool first_interval, // if it is the first interval on the edge
                                         bool last_interval, bool turn_left, bool turn_right,
                                         IntervalWithStop* candidates); // if it is the last interval on the edge

  void construct_propagated_intervals(bool invert, Halfedge halfedge, IntervalWithStop* candidates,
                                      unsigned& num_candidates, interval_pointer source_interval);
  // constructs iNew from the rest of the data

  double compute_positive_intersection(double start, double pseudo_x, double pseudo_y, double sin_alpha,
                                       double cos_alpha); // used in construct_propagated_intervals

  // intersecting two intervals with up to three intervals in the end
  unsigned intersect_intervals(interval_pointer zero, IntervalWithStop* one);

  const_interval_pointer best_first_interval(const SurfacePoint& point, double& best_total_distance,
                                             double& best_interval_position, unsigned& best_source_index) const;

  bool check_stop_conditions(unsigned& index) const;

  void clear();

  list_pointer interval_list(Edge e) { return &m_edge_interval_lists[e]; };
  const_list_pointer interval_list(Edge e) const { return &m_edge_interval_lists[e]; };

  void set_sources(const std::vector<SurfacePoint>& sources) { m_sources.initialize(sources); }

  void initialize_propagation_data();

  // used in initialization
  void list_edges_visible_from_source(const SurfacePoint& source, std::vector<Edge>& storage) const;

  long visible_from_source(const SurfacePoint& point) const; // used in backtracing

  void best_point_on_the_edge_set(const SurfacePoint& point, std::vector<Edge> const& storage,
                                  const_interval_pointer& best_interval, double& best_total_distance,
                                  double& best_interval_position, bool verbose = false) const;

  void possible_traceback_edges(const SurfacePoint& point, std::vector<Edge>& storage) const;

  bool erase_from_queue(interval_pointer p);

  void set_stop_conditions(const std::vector<SurfacePoint>& stop_points, double stop_distance);
  double stop_distance() const { return m_max_propagation_distance; }

  //== Data
  typedef std::pair<Vertex, double> stop_vertex_with_distace_type;
  // algorithm stops propagation after covering certain vertices
  std::vector<stop_vertex_with_distace_type> m_stop_vertices;
  double m_max_propagation_distance; // or reaching the certain distance

  SurfaceMesh& mesh;
  IntrinsicGeometryInterface& geom;

  double m_time_consumed;                // how much time does the propagation step takes
  double m_propagation_distance_stopped; // at what distance (if any) the
                                         // propagation algorithm stopped

  IntervalQueue m_queue; // interval queue

  MemoryAllocator<Interval> m_memory_allocator; // quickly allocate and deallocate intervals
  EdgeData<IntervalList> m_edge_interval_lists; // every edge has its interval data

  enum class MapType { OLD, NEW }; // used for interval intersection
  MapType map[5];
  double start[6];
  interval_pointer i_new[5];

  size_t m_queue_max_size; // used for statistics
  size_t m_iterations;     // used for statistics

  SortedSources m_sources;

  VertexData<bool> vertexIsManifold, vertexIsBoundary;
};

} // namespace surface
} // namespace geometrycentral
