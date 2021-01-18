#pragma once

// Would be great not to include Eigen here, but seems necessary to specialize clearBuffer() in the .ipp. Suggestions
// for an easy workaround are welcome.
#include <Eigen/SparseCore>

#include <functional>
#include <iostream>
#include <vector>
#include <array>


namespace geometrycentral {

class DependentQuantity {

public:
  DependentQuantity(std::function<void()> evaluateFunc_, std::vector<DependentQuantity*>& listToJoin)
      : evaluateFunc(evaluateFunc_) {
    listToJoin.push_back(this);
  }

  virtual ~DependentQuantity(){};

  std::function<void()> evaluateFunc;
  bool computed = false;
  int requireCount = 0;
  bool clearable = true; // if false, clearing does nothing

  // Compute the quantity, if we don't have it already
  void ensureHave();

  // Compute the quantity if we need it and don't have it already
  void ensureHaveIfRequired();

  // Note that something will reqiure this quantity (increments a count of such requirements),
  // and ensure that we have this quantity
  void require();

  // Decrement the count of requirements of this quantity
  void unrequire();

  // Clear out the underlying quantity to reduce memory usage
  virtual void clearIfNotRequired() = 0;
};

// Wrapper class which manages a dependency graph of quantities. Templated on the underlying type of the data.
template <typename D>
class DependentQuantityD : public DependentQuantity {

public:
  DependentQuantityD(){};
  virtual ~DependentQuantityD(){};

  DependentQuantityD(D* dataBuffer_, std::function<void()> evaluateFunc_, std::vector<DependentQuantity*>& listToJoin)
      : DependentQuantity(evaluateFunc_, listToJoin), dataBuffer(dataBuffer_) {}

  D* dataBuffer = nullptr;

  // Clear out the underlying quantity to reduce memory usage
  virtual void clearIfNotRequired() override;
};

} // namespace geometrycentral

#include "geometrycentral/utilities/dependent_quantity.ipp"
