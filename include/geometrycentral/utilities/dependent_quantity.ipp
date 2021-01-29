namespace geometrycentral {

inline void DependentQuantity::ensureHaveIfRequired() {
  if (requireCount > 0) {
    ensureHave();
  }
}

inline void DependentQuantity::ensureHave() {

  // If the quantity is already populated, early out
  if (computed) {
    return;
  }

  // Compute this quantity
  evaluateFunc();

  computed = true;
};

inline void DependentQuantity::require() {
  requireCount++;
  ensureHave();
}

inline void DependentQuantity::unrequire() {
  requireCount--;

  if (requireCount < 0) {
    throw std::logic_error("Quantity was unrequire()'d more than than it was require()'d");
    requireCount = 0;
  }
}

// Helper functions to clear data
// Note: if/when we start using more types in these quantities, we might need to generalize this mechanism. But for the
// current set of uses (scalars, MeshData<>, Eigen types), this works just fine.
namespace {

// General method: call a clear function
template <typename T>
void clearBuffer(T* buffer) {
  buffer->clear();
}

// Scalars
void clearBuffer(double* buffer) {}
void clearBuffer(size_t* buffer) {}
void clearBuffer(int* buffer) {}


// Eigen sparse matrices
template <typename F>
void clearBuffer(Eigen::SparseMatrix<F>* buffer) {
  *buffer = Eigen::SparseMatrix<F>();
}

// any unique_ptr<> type
template <typename P>
void clearBuffer(std::unique_ptr<P>* buffer) {
  buffer->reset();
}

// Array of any otherwise clearable type
template <typename A, size_t N>
void clearBuffer(std::array<A*, N>* buffer) {
  for (size_t i = 0; i < N; i++) {
    // Recurse to an approriate version of this template
    A* elem = (*buffer)[i];
    clearBuffer(elem);
  }
}

// Pair of clearable types
template <typename A, typename B>
void clearBuffer(std::pair<A, B>* buffer) {
  A elemA = buffer->first;
  clearBuffer(elemA);
  B elemB = buffer->second;
  clearBuffer(elemB);
}

} // namespace

template <typename D>
void DependentQuantityD<D>::clearIfNotRequired() {
  if (clearable && requireCount <= 0 && dataBuffer != nullptr && computed) {
    clearBuffer(dataBuffer);
    computed = false;
  }
}

} // namespace geometrycentral
