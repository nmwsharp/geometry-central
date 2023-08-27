// Copyright (C) 2008 Danil Kirsanov, MIT License
// (Modified to work in geometry-central. Original code can be found here: https://code.google.com/p/geodesic/)

#pragma once

#include "geometrycentral/surface/barycentric_coordinate_helpers.h"
#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/surface_point.h"

#include <algorithm>
#include <assert.h>
#include <cmath>
#include <vector>

// added by nsharp
#include <memory>
#include <string.h>

namespace geometrycentral {
namespace surface {

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// double const GEODESIC_INF = std::numeric_limits<double>::max();
double const GEODESIC_INF = 1e100;

// in order to avoid numerical problems with "infinitely small" intervals,
// we drop all the intervals smaller than SMALLEST_INTERVAL_RATIO*edge_length
double const SMALLEST_INTERVAL_RATIO = 1e-6;

class Interval;
class IntervalList;
typedef Interval* interval_pointer;
typedef const Interval* const_interval_pointer;
typedef IntervalList* list_pointer;
typedef const IntervalList* const_list_pointer;

// interval of the edge
class Interval {
public:
  Interval(){};
  ~Interval(){};

  enum class DirectionType { FROM_HALFEDGE, FROM_SOURCE, UNDEFINED_DIRECTION };

  // geodesic distance function at point x
  double signal(double x) const;

  double max_distance(double end) const;

  // compute min, given c,d theta, start, end. save value to m_min
  void compute_min_distance(double stop);

  // compare two intervals in the queue
  bool operator()(interval_pointer const x, interval_pointer const y) const;

  // return the endpoint of the interval
  double stop() const;

  double hypotenuse(double a, double b) const;

  // find the point on the interval that is closest to the point (x, y)
  void find_closest_point(double const x, double const y, double& offset, double& distance) const;

  double& start() { return m_start; };
  const double& start() const { return m_start; };
  double& d() { return m_d; };
  double& pseudo_x() { return m_pseudo_x; };
  double& pseudo_y() { return m_pseudo_y; };
  double& min() { return m_min; };
  double min() const { return m_min; };
  interval_pointer& next() { return m_next; };
  const_interval_pointer next() const { return m_next; };
  Edge& edge() { return m_edge; };
  const Edge& edge() const { return m_edge; };
  Halfedge& halfedge() { return m_halfedge; };
  const Halfedge& halfedge() const { return m_halfedge; };
  double& edge_length() { return m_edge_length; };
  DirectionType& direction() { return m_direction; };
  bool visible_from_source() const { return m_direction == DirectionType::FROM_SOURCE; };
  unsigned& source_index() { return m_source_index; };
  unsigned source_index() const { return m_source_index; };

  void initialize(IntrinsicGeometryInterface& geom, Edge edge, double edge_length, SurfacePoint* point = nullptr,
                  unsigned source_index = 0);

protected:
  double m_start;    // initial point of the interval on the edge
  double m_d;        // distance from the source to the pseudo-source
  double m_pseudo_x; // coordinates of the pseudo-source in the local coordinate system
  double m_pseudo_y; // y-coordinate should be always negative
  double m_min;      // minimum distance on the interval

  interval_pointer m_next; // pointer to the next interval in the list
  Edge m_edge;             // edge that the interval belongs to
  Halfedge m_halfedge;     // halfedge indicating which direction the interval comes from (only set if DirectionType is
                           // FROM_HALFEDGE)
  double m_edge_length = -1; // length of m_edge
  unsigned m_source_index;   // the source it belongs to
  DirectionType m_direction; // where the interval is coming from
};

struct IntervalWithStop : public Interval {
public:
  double& stop() { return m_stop; };

public:
  double m_stop;
};

// list of the of intervals of the given edge
class IntervalList {
public:
  IntervalList();
  ~IntervalList(){};

  void clear();
  void initialize(Edge e);

  // returns the interval that covers the offset
  interval_pointer covering_interval(double offset);
  const_interval_pointer covering_interval(double offset) const;

  void find_closest_point(IntrinsicGeometryInterface& geom, const SurfacePoint point, double& offset, double& distance,
                          const_interval_pointer& interval, bool verbose = false) const;

  unsigned number_of_intervals() const;
  interval_pointer last();
  double signal(double x) const;

  interval_pointer& first();
  Edge& edge();

public:
  interval_pointer m_first; // pointer to the first member of the list
  Edge m_edge;              // edge that owns this list
};

class SurfacePointWithIndex : public SurfacePoint {
public:
  SurfacePointWithIndex() : SurfacePoint(){};
  SurfacePointWithIndex(const SurfacePoint& p) : SurfacePoint(p){};
  unsigned index() const;

  void initialize(const SurfacePoint& p, unsigned index);

  // used for sorting
  bool operator()(const SurfacePointWithIndex* x, const SurfacePointWithIndex* y) const;
  bool operator()(const SurfacePointWithIndex& x, const SurfacePointWithIndex* y) const;
  bool operator()(const SurfacePointWithIndex* x, const SurfacePointWithIndex& y) const;
  bool compare(const SurfacePointWithIndex& x, const SurfacePointWithIndex& y) const;

public:
  unsigned m_index;
};

class SortedSources : public std::vector<SurfacePointWithIndex> {
public:
  typedef std::vector<SurfacePointWithIndex*> sorted_vector_type;

public:
  typedef sorted_vector_type::iterator sorted_iterator;
  typedef std::pair<sorted_iterator, sorted_iterator> sorted_iterator_pair;
  typedef sorted_vector_type::const_iterator const_sorted_iterator;
  typedef std::pair<const_sorted_iterator, const_sorted_iterator> const_sorted_iterator_pair;

  sorted_iterator_pair sources(const SurfacePoint& mesh_element);
  const_sorted_iterator_pair sources(const SurfacePoint& mesh_element) const;

  // we initialize the sources by copy
  void initialize(const std::vector<SurfacePoint>& sources);

  SurfacePointWithIndex& operator[](unsigned i);
  const SurfacePointWithIndex& operator[](unsigned i) const;

public:
  sorted_vector_type m_sorted;
  SurfacePointWithIndex m_search_dummy; // used as a search template
  SurfacePointWithIndex m_compare_less; // used as a compare functor
};

// Assorted helper functions
namespace exactgeodesic {
double compute_surface_distance(IntrinsicGeometryInterface& geom, const SurfacePoint& p1, const SurfacePoint& p2);

unsigned compute_closest_vertices(SurfacePoint p, std::vector<Vertex>* storage);

std::pair<double, double> compute_local_coordinates(IntrinsicGeometryInterface& geom, Edge e,
                                                    const SurfacePoint& point);

// maps a tangent vector v from f's coordinate system to p's coordinate system. p must be located on f, or one of its
// vertices or edges
Vector2 transformToCoordinateSystem(IntrinsicGeometryInterface& geom, Vector2 v, Face f, SurfacePoint p);
} // namespace exactgeodesic

//== A fast and simple memory allocator
// quickly allocates and deallocates single elements of a given type
template <class T>
class MemoryAllocator {
public:
  typedef T* pointer;

  MemoryAllocator(unsigned block_size = 1024, unsigned max_number_of_blocks = 1024) {
    reset(block_size, max_number_of_blocks);
  }
  ~MemoryAllocator(){};

  void clear();

  void reset(unsigned block_size, unsigned max_number_of_blocks);

  // allocates single unit of memory
  pointer allocate();

  // allocate n units
  void deallocate(pointer p);

protected:
  std::vector<std::vector<T>> m_storage;
  unsigned m_block_size;           // size of a single block
  unsigned m_max_number_of_blocks; // maximum allowed number of blocks
  unsigned m_current_position;     // first unused element inside the current
                                   // block

  std::vector<pointer> m_deleted; // pointers to deleted elemets
};

} // namespace surface
} // namespace geometrycentral

#include "geometrycentral/surface/exact_geodesic_helpers.ipp"
