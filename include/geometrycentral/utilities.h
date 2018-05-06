#pragma once

#include <cmath>
#include <complex>
#include <complex>
#include <exception>
#include <limits>
#include <random>
#include <string>
#include <typeinfo>

#include "geometrycentral/vector2.h"
#include "geometrycentral/vector3.h"

namespace geometrycentral {

const double PI = 3.1415926535897932384;

// Lightweight class for holding 3-tuples of unsigned ints
union uint3 {
  struct {
    unsigned int first;
    unsigned int second;
    unsigned int third;
  };
  unsigned int v[3];
};

// Complex numbers are useful
using Complex = ::std::complex<double>;
const Complex IM_I(0.0, 1.0);
inline double dot(Complex x, Complex y) {
  return x.real() * y.real() + x.imag() * y.imag();
}

inline Complex inv(Complex c) { return ::std::conj(c) / ::std::norm(c); }
inline Complex unit(Complex c) { return c / ::std::abs(c); }

// Various functions
template <typename T>
T clamp(T val, T low, T high);
Vector3 clamp(Vector3 val, Vector3 low, Vector3 high);
template <typename T>
bool approxEqualsAbsolute(T a, T b, double eps = 1e-6);
double regularizeAngle(double theta);  // Map theta in to [0,2pi)

template <typename T>
T sqr(T x) {
  return x * x;
}

// === Inline implementations
template <typename T>
inline T clamp(T val, T low, T high) {
  if (val > high) return high;
  if (val < low) return low;
  return val;
}
inline Vector3 clamp(Vector3 val, Vector3 low, Vector3 high) {
  double x = clamp(val.x, low.x, high.x);
  double y = clamp(val.y, low.y, high.y);
  double z = clamp(val.z, low.z, high.z);
  return Vector3{x, y, z};
}

template <typename T>
inline bool approxEqualsAbsolute(T a, T b, double eps) {
  double absA = ::std::abs(a);
  double absB = ::std::abs(b);
  double absDiff = ::std::abs(a - b);

  if (a == b) {
    return true;
  } else {
    return absDiff < eps;
  }
}

inline double regularizeAngle(double theta) {
  return theta - 2 * PI * ::std::floor(theta / (2 * PI));
}

template <typename T>
std::string typeNameString(T& x) {
  return std::string(typeid(x).name());
}

template <typename T>
std::string typeNameString(T* x) {
  return std::string(typeid(x).name());
}

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

// Random number generation -----------------------------------------
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
  std::uniform_int_distribution<size_t> dist(0, size-1);
  return dist(util_mersenne_twister);
}

inline double randomNormal(double mean=0.0, double stddev=1.0) {
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

// === Custom error types
class FunctionalityException : public std::runtime_error {
 public:
  FunctionalityException(std::string msg)
      : std::runtime_error("Missing functionaliy: " + msg){};
};

}  // namespace geometrycentral

namespace std {
// NOTE: Technically, the lines below are illegal, because specializing standard
// library functions is ONLY allowed for user-defined types. That being said, I
// don't think this will crash any planes.
inline bool isfinite(const ::std::complex<double> c) {
  return isfinite(c.real()) && isfinite(c.imag());
}
}
