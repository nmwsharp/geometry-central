#pragma once

#include "geometrycentral/surface/mutation_manager.h"

namespace geometrycentral {
namespace surface {

// ======================================================
// ======== Helpers to expose the "simple" interface
// ======================================================

template <typename T>
class SimpleFlipUpdater : public EdgeFlipPolicy {
  // General case where the pre function passes data to the post function. Both must be given.
  T data;

public:
  std::function<T(Edge)> pre;
  std::function<void(Edge, T)> post;
  void beforeEdgeFlip(const Edge& e) override { data = pre(e); }
  void afterEdgeFlip(const Edge& e) override { post(e, data); }
};
template <>
class SimpleFlipUpdater<void> : public EdgeFlipPolicy {
  // Special case where no data is passed from pre to post. Just one of the functions may be specified.
public:
  std::function<void(Edge)> pre;
  std::function<void(Edge)> post;
  void beforeEdgeFlip(const Edge& e) override {
    if (pre) pre(e);
  }
  void afterEdgeFlip(const Edge& e) override {
    if (post) post(e);
  }
};

template <typename T>
class SimpleSplitUpdater : public EdgeSplitPolicy {
  // General case where the pre function passes data to the post function. Both must be given.
  T data;
public:
  std::function<T(Edge, double)> pre;
  std::function<void(Halfedge, Halfedge, double, T)> post;
  void beforeEdgeSplit(const Edge& e, const double& tSplit) override { data = pre(e, tSplit); }
  void afterEdgeSplit(const Halfedge& he1, const Halfedge& he2, const double& tSplit) override {
    post(he1, he2, tSplit, data);
  }
};
template <>
class SimpleSplitUpdater<void> : public EdgeSplitPolicy {
  // Special case where no data is passed from pre to post. Just one of the functions may be specified.
public:
  std::function<void(Edge, double)> pre;
  std::function<void(Halfedge, Halfedge, double)> post;
  void beforeEdgeSplit(const Edge& e, const double& tSplit) override {
    if (pre) pre(e, tSplit);
  }
  void afterEdgeSplit(const Halfedge& he1, const Halfedge& he2, const double& tSplit) override {
    if (post) post(he1, he2, tSplit);
  }
};

template <typename T>
class SimpleFaceSplitUpdater : public FaceSplitPolicy {
  // General case where the pre function passes data to the post function. Both must be given.
  T data;

public:
  std::function<T(Face, std::vector<double>)> pre;
  std::function<void(Vertex, T)> post;
  void beforeFaceSplit(const Face& f, const std::vector<double>& bary) override { data = pre(f, bary); }
  void afterFaceSplit(const Vertex& newV, const std::vector<double>& bary) override { post(newV, data); }
};
template <>
class SimpleFaceSplitUpdater<void> : public FaceSplitPolicy {
  // Special case where no data is passed from pre to post. Just one of the functions may be specified.
public:
  std::function<void(Face, std::vector<double>)> pre;
  std::function<void(Vertex)> post;
  void beforeFaceSplit(const Face& f, const std::vector<double>& bary) override {
    if (pre) pre(f, bary);
  }
  void afterFaceSplit(const Vertex& newV, const std::vector<double>& bary) override {
    if (post) post(newV);
  }
};

template <typename T>
class SimpleEdgeCollapseUpdater : public EdgeCollapsePolicy {
  // General case where the pre function passes data to the post function. Both must be given.
  T data;

public:
  std::function<T(Edge, double)> pre;
  std::function<void(Vertex, double, T)> post;
  void beforeEdgeCollapse(const Edge& e, const double& tCollapse) override { data = pre(e, tCollapse); }
  void afterEdgeCollapse(const Vertex& newV, const double& tCollapse) override { post(newV, tCollapse, data); }
};
template <>
class SimpleEdgeCollapseUpdater<void> : public EdgeCollapsePolicy {
  // Special case where no data is passed from pre to post. Just one of the functions may be specified.
public:
  std::function<void(Edge, double)> pre;
  std::function<void(Vertex, double)> post;
  void beforeEdgeCollapse(const Edge& e, const double& tCollapse) override {
    if (pre) pre(e, tCollapse);
  }
  void afterEdgeCollapse(const Vertex& newV, const double& tCollapse) override {
    if (post) post(newV, tCollapse);
  }
};

// TODO add some ASSERT compile-time checks to these, so that in the likely case where a user passes a non-matching
// template function, we can give a sane diagnostic message.


// === Edge flip handler one-liners
template <typename Fpre, typename Fpost>
void MutationManager::registerEdgeFlipHandlers(Fpre pre, Fpost post) {
  using T = decltype(pre(Edge()));
  SimpleFlipUpdater<T>* updater = new SimpleFlipUpdater<T>();
  updater->pre = pre;
  updater->post = post;
  registerPolicy(updater);
}
template <typename Fpre>
void MutationManager::registerEdgeFlipPreHandler(Fpre pre) {
  SimpleFlipUpdater<void>* updater = new SimpleFlipUpdater<void>();
  updater->pre = pre;
  registerPolicy(updater);
}
template <typename Fpost>
void MutationManager::registerEdgeFlipPostHandler(Fpost post) {
  SimpleFlipUpdater<void>* updater = new SimpleFlipUpdater<void>();
  updater->post = post;
  registerPolicy(updater);
}

// === Edge split handler one-liners
template <typename Fpre, typename Fpost>
void MutationManager::registerEdgeSplitHandlers(Fpre pre, Fpost post) {
  using T = decltype(pre(Edge(), 0.5));
  SimpleSplitUpdater<T>* updater = new SimpleSplitUpdater<T>();
  updater->pre = pre;
  updater->post = post;
  registerPolicy(updater);
}
template <typename Fpre>
void MutationManager::registerEdgeSplitPreHandler(Fpre pre) {
  SimpleSplitUpdater<void>* updater = new SimpleSplitUpdater<void>();
  updater->pre = pre;
  registerPolicy(updater);
}
template <typename Fpost>
void MutationManager::registerEdgeSplitPostHandler(Fpost post) {
  SimpleSplitUpdater<void>* updater = new SimpleSplitUpdater<void>();
  updater->post = post;
  registerPolicy(updater);
}

// === Vertex insertion handler one-liners
template <typename Fpre, typename Fpost>
void MutationManager::registerFaceSplitHandlers(Fpre pre, Fpost post) {
  using T = decltype(pre(Face(), std::vector<double>{0, 0, 0}));
  SimpleFaceSplitUpdater<T>* updater = new SimpleFaceSplitUpdater<T>();
  updater->pre = pre;
  updater->post = post;
  registerPolicy(updater);
}
template <typename Fpre>
void MutationManager::registerFaceSplitPreHandler(Fpre pre) {
  SimpleFaceSplitUpdater<void>* updater = new SimpleFaceSplitUpdater<void>();
  updater->pre = pre;
  registerPolicy(updater);
}
template <typename Fpost>
void MutationManager::registerFaceSplitPostHandler(Fpost post) {
  SimpleFaceSplitUpdater<void>* updater = new SimpleFaceSplitUpdater<void>();
  updater->post = post;
  registerPolicy(updater);
}

// === Edge collapse handler one-liners
template <typename Fpre, typename Fpost>
void MutationManager::registerEdgeCollapseHandlers(Fpre pre, Fpost post) {
  using T = decltype(pre(Edge(), 0.5));
  SimpleEdgeCollapseUpdater<T>* updater = new SimpleEdgeCollapseUpdater<T>();
  updater->pre = pre;
  updater->post = post;
  registerPolicy(updater);
}
template <typename Fpre>
void MutationManager::registerEdgeCollapsePreHandler(Fpre pre) {
  SimpleEdgeCollapseUpdater<void>* updater = new SimpleEdgeCollapseUpdater<void>();
  updater->pre = pre;
  registerPolicy(updater);
}
template <typename Fpost>
void MutationManager::registerEdgeCollapsePostHandler(Fpost post) {
  SimpleEdgeCollapseUpdater<void>* updater = new SimpleEdgeCollapseUpdater<void>();
  updater->post = post;
  registerPolicy(updater);
}


} // namespace surface
} // namespace geometrycentral
