#pragma once

#include <algorithm> // std::remove
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <string>
#include <typeinfo>

// Error checking macro. CONDITION should be true if life is good (like in assert(CONDITION))
#ifdef NGC_SAFETY_CHECKS
#define GC_SAFETY_ASSERT(CONDITION, MSG)
#else
#define GC_SAFETY_ASSERT(CONDITION, MSG)                                                                               \
  if (!(CONDITION)) {                                                                                                  \
    throw std::runtime_error("GC_SAFETY_ASSERT FAILURE from " + std::string(__FILE__) + ":" +                          \
                             std::to_string(__LINE__) + " - " + (MSG));                                                \
  }
#endif

namespace geometrycentral {

// === Constants

const size_t INVALID_IND = std::numeric_limits<size_t>::max();
const double PI = 3.1415926535897932384;

// === Memory management

template <typename T>
void safeDelete(T*& x) {
  if (x != nullptr) {
    delete x;
    x = nullptr;
  }
}

template <typename T>
void safeDeleteArray(T*& x) {
  if (x != nullptr) {
    delete[] x;
    x = nullptr;
  }
}

// === Type names

template <typename T>
std::string typeNameString(T& x) {
  return std::string(typeid(x).name());
}

template <typename T>
std::string typeNameString(T* x) {
  return std::string(typeid(x).name());
}

// === Arithmetic

// Clamp
template <typename T>
T clamp(T val, T low, T high);

double regularizeAngle(double theta); // Map theta in to [0,2pi)

// Missing isfinite function
inline bool isfinite(const std::complex<double>& c) { return std::isfinite(c.real()) && std::isfinite(c.imag()); }

// Conjugate which does nothing to non-complex values
template <typename T>
inline T conj(T val) {
  return val;
}
template <>
inline std::complex<double> conj(std::complex<double> val) {
  return std::conj(val);
}

// === Inline implementations
template <typename T>
inline T clamp(T val, T low, T high) {
  if (val > high) return high;
  if (val < low) return low;
  return val;
}

inline double regularizeAngle(double theta) { return theta - 2 * PI * ::std::floor(theta / (2 * PI)); }

// Applies a permutation such that d_new[i] = d_old[p[i]].
// Permutation should be injection to [0,sourceData.size()). Return vector has length permNewToOld.size().
// Any permutation indices with value INVALID_IND are skipped
template <typename T, typename A1, typename A2>
std::vector<T, A1> applyPermutation(const std::vector<T, A1>& sourceData, const std::vector<size_t, A2>& permNewToOld) {
  std::vector<T, A1> retVal(permNewToOld.size());
  for (size_t i = 0; i < permNewToOld.size(); i++) {
    if (permNewToOld[i] != INVALID_IND) {
      retVal[i] = sourceData[permNewToOld[i]];
    }
  }
  return retVal;
}

// erase-remove idiom
// modifies vector in-place removing all occurences of `obj` according to == (if any)
template <typename T, typename O>
void removeFromVector(std::vector<T>& vec, O& obj) {
  vec.erase(std::remove(vec.begin(), vec.end(), obj), vec.end());
}

// === Random number generation ===

extern std::random_device util_random_device;
extern std::mt19937 util_mersenne_twister;

inline double unitRand() {
  std::uniform_real_distribution<double> dist(0., 1.);
  return dist(util_mersenne_twister);
}

inline double randomReal(double minVal, double maxVal) {
  std::uniform_real_distribution<double> dist(minVal, maxVal);
  return dist(util_mersenne_twister);
}

// Generate a random int in the INCLUSIVE range [lower,upper]
inline int randomInt(int lower, int upper) {
  std::uniform_int_distribution<int> dist(lower, upper);
  return dist(util_mersenne_twister);
}
// Generate a random size_t in the range [0, N)
inline size_t randomIndex(size_t size) {
  std::uniform_int_distribution<size_t> dist(0, size - 1);
  return dist(util_mersenne_twister);
}

inline double randomNormal(double mean = 0.0, double stddev = 1.0) {
  std::normal_distribution<double> dist{mean, stddev};
  return dist(util_mersenne_twister);
}

// === Printing things to strings ===
template <typename T>
inline std::string to_string(std::vector<T> const& v) {
  std::stringstream ss;
  ss << "[";
  for (size_t i = 0; i < v.size(); i++) {
    if (i > 0) {
      ss << ",";
    }
    ss << v[i];
  }
  ss << "]";

  return ss.str();
}

// Printf to a std::string
template <typename... Args>
std::string str_printf(const std::string& format, Args... args) {
  size_t size = std::snprintf(nullptr, 0, format.c_str(), args...) + 1;
  std::unique_ptr<char[]> buf(new char[size]);
  std::snprintf(buf.get(), size, format.c_str(), args...);
  return std::string(buf.get(), buf.get() + size - 1);
}

} // namespace geometrycentral
