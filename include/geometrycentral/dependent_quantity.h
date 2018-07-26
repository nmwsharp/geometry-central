#pragma once

#include <functional>
#include <vector>

namespace geometrycentral {

// Helper class which manages a dependency graph of quantities
class DependentQuantity {

public:
  DependentQuantity(){};

  DependentQuantity(std::vector<DependentQuantity*> dependencies_, std::function<void()> evaluateFunc_)
      : dependencies(dependencies_), evaluateFunc(evaluateFunc_) {}

  std::vector<DependentQuantity*> dependencies;

  // Compute the quantity, if we don't have it already
  void ensureHave();

  // Compute the quantity if we need it and don't have it already
  void ensureHaveIfRequired();

  // Note that something will reqiure this quantity (increments a count of such requirements),
  // and ensure that we have this quantity
  void require();

  // Decrement the count of requirements of this quantity
  void unrequire();

  bool computed = false;
  int requireCount = 0;

  std::function<void()> evaluateFunc;
};
} // namespace geometrycentral
