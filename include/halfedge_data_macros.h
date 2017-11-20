#pragma once

// These operators implement element-wise scalar arithmetic for VertexData<>
// (etc) types.
// Since they are exactly identical between different containers, the macros
// reduce code
// duplication and keeps things a bit tidier.

#define GC_INTERNAL_GENERATE_DATATYPE_OPERATOR_DECLARATIONS(XXXDATA)           \
  /* Operators */                                                              \
  /* Note: Technically these friend declarations are more broad than they need \
   * to be, as U need not match T. */                                          \
  /*       Doesn't seem to be any better method, however. */                   \
  XXXDATA<T>& operator*=(const XXXDATA<T>& rhs);                               \
  XXXDATA<T>& operator*=(const T& scalar);                                     \
  template <typename U>                                                        \
  friend XXXDATA<U> operator*(const XXXDATA<U>& lhs, const XXXDATA<U>& rhs);   \
  template <typename U>                                                        \
  friend XXXDATA<U> operator*(const XXXDATA<U>& lhs, const U& scalar);         \
  template <typename U>                                                        \
  friend XXXDATA<U> operator*(const U& scalar, const XXXDATA<U>& rhs);         \
                                                                               \
  XXXDATA<T>& operator/=(const XXXDATA<T>& rhs);                               \
  XXXDATA<T>& operator/=(const T& scalar);                                     \
  template <typename U>                                                        \
  friend XXXDATA<U> operator/(const XXXDATA<U>& lhs, const XXXDATA<U>& rhs);   \
  template <typename U>                                                        \
  friend XXXDATA<U> operator/(const XXXDATA<U>& lhs, const U& scalar);         \
  template <typename U>                                                        \
  friend XXXDATA<U> operator/(const U& scalar, const XXXDATA<U>& rhs);         \
                                                                               \
  XXXDATA<T>& operator+=(const XXXDATA<T>& rhs);                               \
  XXXDATA<T>& operator+=(const T& scalar);                                     \
  template <typename U>                                                        \
  friend XXXDATA<U> operator+(const XXXDATA<U>& lhs, const XXXDATA<U>& rhs);   \
  template <typename U>                                                        \
  friend XXXDATA<U> operator+(const XXXDATA<U>& lhs, const U& scalar);         \
  template <typename U>                                                        \
  friend XXXDATA<U> operator+(const U& scalar, const XXXDATA<U>& rhs);         \
                                                                               \
  XXXDATA<T>& operator-=(const XXXDATA<T>& rhs);                               \
  XXXDATA<T>& operator-=(const T& scalar);                                     \
  template <typename U>                                                        \
  friend XXXDATA<U> operator-(const XXXDATA<U>& lhs, const XXXDATA<U>& rhs);   \
  template <typename U>                                                        \
  friend XXXDATA<U> operator-(const XXXDATA<U>& lhs, const U& scalar);         \
  template <typename U>                                                        \
  friend XXXDATA<U> operator-(const U& scalar, const XXXDATA<U>& rhs);         \
                                                                               \
  XXXDATA<T>& operator%=(const XXXDATA<T>& rhs);                               \
  XXXDATA<T>& operator%=(const T& scalar);                                     \
  template <typename U>                                                        \
  friend XXXDATA<U> operator%(const XXXDATA<U>& lhs, const XXXDATA<U>& rhs);   \
  template <typename U>                                                        \
  friend XXXDATA<U> operator%(const XXXDATA<U>& lhs, const U& scalar);         \
  template <typename U>                                                        \
  friend XXXDATA<U> operator%(const U& scalar, const XXXDATA<U>& rhs);         \
                                                                               \
  XXXDATA<T>& operator-();                                                     \
                                                                               \
  template <typename U>                                                        \
  friend XXXDATA<U> operator&&(const XXXDATA<U>& lhs, const XXXDATA<U>& rhs);  \
  template <typename U>                                                        \
  friend XXXDATA<U> operator&&(const XXXDATA<U>& lhs, const U& scalar);        \
  template <typename U>                                                        \
  friend XXXDATA<U> operator&&(const U& scalar, const XXXDATA<U>& rhs);        \
                                                                               \
  template <typename U>                                                        \
  friend XXXDATA<U> operator||(const XXXDATA<U>& lhs, const XXXDATA<U>& rhs);  \
  template <typename U>                                                        \
  friend XXXDATA<U> operator||(const XXXDATA<U>& lhs, const U& scalar);        \
  template <typename U>                                                        \
  friend XXXDATA<U> operator||(const U& scalar, const XXXDATA<U>& rhs);        \
                                                                               \
  XXXDATA<T>& operator!();

#define GC_INTERNAL_GENERATE_DATATYPE_OPERATOR_DEFINITIONS(XXXDATA, MESHVAR)   \
  /* Multiplication (*) */                                                     \
  template <typename T>                                                        \
  XXXDATA<T>& XXXDATA<T>::operator*=(const XXXDATA<T>& rhs) {                  \
    if (MESHVAR != rhs.MESHVAR)                                                \
      throw std::runtime_error(                                                \
          "Meshes must match for arithmetic operations.");                     \
    for (size_t i = 0; i < data.size(); i++) data[i] *= rhs.data[i];           \
    return *this;                                                              \
  }                                                                            \
  template <typename T>                                                        \
  inline XXXDATA<T> operator*(const XXXDATA<T>& lhs, const XXXDATA<T>& rhs) {  \
    if (lhs.MESHVAR != rhs.MESHVAR)                                            \
      throw std::runtime_error(                                                \
          "Meshes must match for arithmetic operations.");                     \
    XXXDATA<T> result(lhs.MESHVAR);                                            \
    for (size_t i = 0; i < result.data.size(); i++)                            \
      result.data[i] = lhs.data[i] * rhs.data[i];                              \
    return result;                                                             \
  }                                                                            \
  template <typename T>                                                        \
  XXXDATA<T>& XXXDATA<T>::operator*=(const T& scalar) {                        \
    for (size_t i = 0; i < data.size(); i++) data[i] *= scalar;                \
    return *this;                                                              \
  }                                                                            \
  template <typename T>                                                        \
  inline XXXDATA<T> operator*(const T& scalar, const XXXDATA<T>& rhs) {        \
    XXXDATA<T> result(rhs.MESHVAR);                                            \
    for (size_t i = 0; i < result.data.size(); i++)                            \
      result.data[i] = scalar * rhs.data[i];                                   \
    return result;                                                             \
  }                                                                            \
  template <typename T>                                                        \
  inline XXXDATA<T> operator*(const XXXDATA<T>& lhs, const T& scalar) {        \
    XXXDATA<T> result(lhs.MESHVAR);                                            \
    for (size_t i = 0; i < result.data.size(); i++)                            \
      result.data[i] = lhs.data[i] * scalar;                                   \
    return result;                                                             \
  }                                                                            \
                                                                               \
  /* Division (/) */                                                           \
  template <typename T>                                                        \
  XXXDATA<T>& XXXDATA<T>::operator/=(const XXXDATA<T>& rhs) {                  \
    if (MESHVAR != rhs.MESHVAR)                                                \
      throw std::runtime_error(                                                \
          "Meshes must match for arithmetic operations.");                     \
    for (size_t i = 0; i < data.size(); i++) data[i] /= rhs.data[i];           \
    return *this;                                                              \
  }                                                                            \
  template <typename T>                                                        \
  inline XXXDATA<T> operator/(const XXXDATA<T>& lhs, const XXXDATA<T>& rhs) {  \
    if (lhs.MESHVAR != rhs.MESHVAR)                                            \
      throw std::runtime_error(                                                \
          "Meshes must match for arithmetic operations.");                     \
    XXXDATA<T> result(lhs.MESHVAR);                                            \
    for (size_t i = 0; i < result.data.size(); i++)                            \
      result.data[i] = lhs.data[i] / rhs.data[i];                              \
    return result;                                                             \
  }                                                                            \
  template <typename T>                                                        \
  XXXDATA<T>& XXXDATA<T>::operator/=(const T& scalar) {                        \
    for (size_t i = 0; i < data.size(); i++) data[i] /= scalar;                \
    return *this;                                                              \
  }                                                                            \
  template <typename T>                                                        \
  inline XXXDATA<T> operator/(const T& scalar, const XXXDATA<T>& rhs) {        \
    XXXDATA<T> result(rhs.MESHVAR);                                            \
    for (size_t i = 0; i < result.data.size(); i++)                            \
      result.data[i] = scalar / rhs.data[i];                                   \
    return result;                                                             \
  }                                                                            \
  template <typename T>                                                        \
  inline XXXDATA<T> operator/(const XXXDATA<T>& lhs, const T& scalar) {        \
    XXXDATA<T> result(lhs.MESHVAR);                                            \
    for (size_t i = 0; i < result.data.size(); i++)                            \
      result.data[i] = lhs.data[i] / scalar;                                   \
    return result;                                                             \
  }                                                                            \
                                                                               \
  /* Addition (+) */                                                           \
  template <typename T>                                                        \
  XXXDATA<T>& XXXDATA<T>::operator+=(const XXXDATA<T>& rhs) {                  \
    if (MESHVAR != rhs.MESHVAR)                                                \
      throw std::runtime_error(                                                \
          "Meshes must match for arithmetic operations.");                     \
    for (size_t i = 0; i < data.size(); i++) data[i] += rhs.data[i];           \
    return *this;                                                              \
  }                                                                            \
  template <typename T>                                                        \
  inline XXXDATA<T> operator+(const XXXDATA<T>& lhs, const XXXDATA<T>& rhs) {  \
    if (lhs.MESHVAR != rhs.MESHVAR)                                            \
      throw std::runtime_error(                                                \
          "Meshes must match for arithmetic operations.");                     \
    XXXDATA<T> result(lhs.MESHVAR);                                            \
    for (size_t i = 0; i < result.data.size(); i++)                            \
      result.data[i] = lhs.data[i] + rhs.data[i];                              \
    return result;                                                             \
  }                                                                            \
  template <typename T>                                                        \
  XXXDATA<T>& XXXDATA<T>::operator+=(const T& scalar) {                        \
    for (size_t i = 0; i < data.size(); i++) data[i] += scalar;                \
    return *this;                                                              \
  }                                                                            \
  template <typename T>                                                        \
  inline XXXDATA<T> operator+(const T& scalar, const XXXDATA<T>& rhs) {        \
    XXXDATA<T> result(rhs.MESHVAR);                                            \
    for (size_t i = 0; i < result.data.size(); i++)                            \
      result.data[i] = scalar + rhs.data[i];                                   \
    return result;                                                             \
  }                                                                            \
  template <typename T>                                                        \
  inline XXXDATA<T> operator+(const XXXDATA<T>& lhs, const T& scalar) {        \
    XXXDATA<T> result(lhs.MESHVAR);                                            \
    for (size_t i = 0; i < result.data.size(); i++)                            \
      result.data[i] = lhs.data[i] + scalar;                                   \
    return result;                                                             \
  }                                                                            \
                                                                               \
  /* Subtraction (-) */                                                        \
  template <typename T>                                                        \
  XXXDATA<T>& XXXDATA<T>::operator-=(const XXXDATA<T>& rhs) {                  \
    if (MESHVAR != rhs.MESHVAR)                                                \
      throw std::runtime_error(                                                \
          "Meshes must match for arithmetic operations.");                     \
    for (size_t i = 0; i < data.size(); i++) data[i] -= rhs.data[i];           \
    return *this;                                                              \
  }                                                                            \
  template <typename T>                                                        \
  inline XXXDATA<T> operator-(const XXXDATA<T>& lhs, const XXXDATA<T>& rhs) {  \
    if (lhs.MESHVAR != rhs.MESHVAR)                                            \
      throw std::runtime_error(                                                \
          "Meshes must match for arithmetic operations.");                     \
    XXXDATA<T> result(lhs.MESHVAR);                                            \
    for (size_t i = 0; i < result.data.size(); i++)                            \
      result.data[i] = lhs.data[i] - rhs.data[i];                              \
    return result;                                                             \
  }                                                                            \
  template <typename T>                                                        \
  XXXDATA<T>& XXXDATA<T>::operator-=(const T& scalar) {                        \
    for (size_t i = 0; i < data.size(); i++) data[i] -= scalar;                \
    return *this;                                                              \
  }                                                                            \
  template <typename T>                                                        \
  inline XXXDATA<T> operator-(const T& scalar, const XXXDATA<T>& rhs) {        \
    XXXDATA<T> result(rhs.MESHVAR);                                            \
    for (size_t i = 0; i < result.data.size(); i++)                            \
      result.data[i] = scalar - rhs.data[i];                                   \
    return result;                                                             \
  }                                                                            \
  template <typename T>                                                        \
  inline XXXDATA<T> operator-(const XXXDATA<T>& lhs, const T& scalar) {        \
    XXXDATA<T> result(lhs.MESHVAR);                                            \
    for (size_t i = 0; i < result.data.size(); i++)                            \
      result.data[i] = lhs.data[i] - scalar;                                   \
    return result;                                                             \
  }                                                                            \
                                                                               \
  /* Remainder (%) */                                                          \
  template <typename T>                                                        \
  XXXDATA<T>& XXXDATA<T>::operator%=(const XXXDATA<T>& rhs) {                  \
    if (MESHVAR != rhs.MESHVAR)                                                \
      throw std::runtime_error(                                                \
          "Meshes must match for arithmetic operations.");                     \
    for (size_t i = 0; i < data.size(); i++) data[i] %= rhs.data[i];           \
    return *this;                                                              \
  }                                                                            \
  template <typename T>                                                        \
  inline XXXDATA<T> operator%(const XXXDATA<T>& lhs, const XXXDATA<T>& rhs) {  \
    if (lhs.MESHVAR != rhs.MESHVAR)                                            \
      throw std::runtime_error(                                                \
          "Meshes must match for arithmetic operations.");                     \
    XXXDATA<T> result(lhs.MESHVAR);                                            \
    for (size_t i = 0; i < result.data.size(); i++)                            \
      result.data[i] = lhs.data[i] % rhs.data[i];                              \
    return result;                                                             \
  }                                                                            \
  template <typename T>                                                        \
  XXXDATA<T>& XXXDATA<T>::operator%=(const T& scalar) {                        \
    for (size_t i = 0; i < data.size(); i++) data[i] %= scalar;                \
    return *this;                                                              \
  }                                                                            \
  template <typename T>                                                        \
  inline XXXDATA<T> operator%(const T& scalar, const XXXDATA<T>& rhs) {        \
    XXXDATA<T> result(rhs.MESHVAR);                                            \
    for (size_t i = 0; i < result.data.size(); i++)                            \
      result.data[i] = scalar % rhs.data[i];                                   \
    return result;                                                             \
  }                                                                            \
  template <typename T>                                                        \
  inline XXXDATA<T> operator%(const XXXDATA<T>& lhs, const T& scalar) {        \
    XXXDATA<T> result(lhs.MESHVAR);                                            \
    for (size_t i = 0; i < result.data.size(); i++)                            \
      result.data[i] = lhs.data[i] % scalar;                                   \
    return result;                                                             \
  }                                                                            \
                                                                               \
  /* Unary negation */                                                         \
  template <typename T>                                                        \
  XXXDATA<T>& XXXDATA<T>::operator-() {                                        \
    for (size_t i = 0; i < data.size(); i++) data[i] = -data[i];               \
    return *this;                                                              \
  }                                                                            \
                                                                               \
  /* Logical And (&&) */                                                       \
  template <typename T>                                                        \
  inline XXXDATA<T> operator&&(const XXXDATA<T>& lhs, const XXXDATA<T>& rhs) { \
    if (lhs.MESHVAR != rhs.MESHVAR)                                            \
      throw std::runtime_error(                                                \
          "Meshes must match for arithmetic operations.");                     \
    XXXDATA<T> result(lhs.MESHVAR);                                            \
    for (size_t i = 0; i < result.data.size(); i++)                            \
      result.data[i] = lhs.data[i] && rhs.data[i];                             \
    return result;                                                             \
  }                                                                            \
  template <typename T>                                                        \
  inline XXXDATA<T> operator&&(const T& scalar, const XXXDATA<T>& rhs) {       \
    XXXDATA<T> result(rhs.MESHVAR);                                            \
    for (size_t i = 0; i < result.data.size(); i++)                            \
      result.data[i] = scalar && rhs.data[i];                                  \
    return result;                                                             \
  }                                                                            \
  template <typename T>                                                        \
  inline XXXDATA<T> operator&&(const XXXDATA<T>& lhs, const T& scalar) {       \
    XXXDATA<T> result(lhs.MESHVAR);                                            \
    for (size_t i = 0; i < result.data.size(); i++)                            \
      result.data[i] = lhs.data[i] && scalar;                                  \
    return result;                                                             \
  }                                                                            \
                                                                               \
  /* Logical Or (||) */                                                        \
  template <typename T>                                                        \
  inline XXXDATA<T> operator||(const XXXDATA<T>& lhs, const XXXDATA<T>& rhs) { \
    if (lhs.MESHVAR != rhs.MESHVAR)                                            \
      throw std::runtime_error(                                                \
          "Meshes must match for arithmetic operations.");                     \
    XXXDATA<T> result(lhs.MESHVAR);                                            \
    for (size_t i = 0; i < result.data.size(); i++)                            \
      result.data[i] = lhs.data[i] || rhs.data[i];                             \
    return result;                                                             \
  }                                                                            \
  template <typename T>                                                        \
  inline XXXDATA<T> operator||(const T& scalar, const XXXDATA<T>& rhs) {       \
    XXXDATA<T> result(rhs.MESHVAR);                                            \
    for (size_t i = 0; i < result.data.size(); i++)                            \
      result.data[i] = scalar || rhs.data[i];                                  \
    return result;                                                             \
  }                                                                            \
  template <typename T>                                                        \
  inline XXXDATA<T> operator||(const XXXDATA<T>& lhs, const T& scalar) {       \
    XXXDATA<T> result(lhs.MESHVAR);                                            \
    for (size_t i = 0; i < result.data.size(); i++)                            \
      result.data[i] = lhs.data[i] || scalar;                                  \
    return result;                                                             \
  }                                                                            \
                                                                               \
  /* Unary not */                                                              \
  template <typename T>                                                        \
  XXXDATA<T>& XXXDATA<T>::operator!() {                                        \
    for (size_t i = 0; i < data.size(); i++) data[i] = !data[i];               \
    return *this;                                                              \
  }