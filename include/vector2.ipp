
namespace geometrycentral {

inline Vector2 Vector2::operator+(const Vector2& v) const {
  return Vector2{x + v.x, y + v.y};
}

inline Vector2 Vector2::operator-(const Vector2& v) const {
  return Vector2{x - v.x, y - v.y};
}

inline Vector2 Vector2::operator*(double s) const {
  return Vector2{x * s, y * s};
}

inline Vector2 Vector2::operator/(double s) const {
  const double r = 1. / s;
  return Vector2{x * r, y * r};
}

inline const Vector2 Vector2::operator-() const { return Vector2{-x, -y}; }

inline Vector2 operator*(const double s, const Vector2& v) {
  return Vector2{s * v.x, s * v.y};
}

inline void Vector2::normalize(void) {
  double r = 1. / sqrt(x * x + y * y);
  x *= r;
  y *= r;
}

inline Vector2& Vector2::operator+=(const Vector2& other) {
  x += other.x;
  y += other.y;
  return *this;
}

inline Vector2& Vector2::operator-=(const Vector2& other) {
  x -= other.x;
  y -= other.y;
  return *this;
}

inline Vector2& Vector2::operator*=(const double& s) {
  x *= s;
  y *= s;
  return *this;
}

inline Vector2& Vector2::operator/=(const double& s) {
  x /= s;
  y /= s;
  return *this;
}

inline bool Vector2::operator==(const Vector2& other) const {
  return x == other.x && y == other.y;
}

inline bool Vector2::operator!=(const Vector2& other) const {
  return !(*this == other);
}

inline double norm(const Vector2& v) { return sqrt(v.x * v.x + v.y * v.y); }

inline double norm2(const Vector2& v) { return v.x * v.x + v.y * v.y; }

inline Vector2 unit(const Vector2& v) {
  double n = norm(v);
  return Vector2{v.x / n, v.y / n};
}

inline double dot(const Vector2& u, const Vector2& v) {
  return u.x * v.x + u.y * v.y;
}

inline double angle(const Vector2& u, const Vector2& v) {
  return acos(fmax(-1., fmin(1., dot(unit(u), unit(v)))));
}

inline Vector3 cross(const Vector2& u, const Vector2& v) {
  return Vector3{0., 0., u.x * v.y - u.y * v.x};
}

inline bool Vector2::isFinite() const {
  return std::isfinite(x) && std::isfinite(y);
}

inline bool Vector2::isDefined() const {
  return (!std::isnan(x)) && (!std::isnan(y));
}

inline Vector2 componentwiseMin(const Vector2& u, const Vector2& v) {
  return Vector2{fmin(u.x, v.x), fmin(u.y, v.y)};
}

inline Vector2 componentwiseMax(const Vector2& u, const Vector2& v) {
  return Vector2{fmax(u.x, v.x), fmax(u.y, v.y)};
}

}  // namespace geometrycentral