#include "geometrycentral/dependent_quantity.h"

namespace geometrycentral {

void DependentQuantity::ensureHaveIfRequired() {
  if (requireCount > 0) {
    ensureHave();
  }
}

void DependentQuantity::ensureHave() {

  // If the quantity is already populated, early out
  if (computed) {
    return;
  }

  // Resolve all of the dependencies
  for (auto x : dependencies) {
    x->ensureHave();
  }

  // Compute this quantity
  evaluateFunc();

  computed = true;
};

void DependentQuantity::require() {
  requireCount++;
  ensureHave();
}

void DependentQuantity::unrequire() {
  requireCount--;

  if (requireCount < 0) {
    throw std::logic_error("Quantity was unrequire()'d more than than it was require()'d");
    requireCount = 0;
  }
}

} // namespace geometrycentral
