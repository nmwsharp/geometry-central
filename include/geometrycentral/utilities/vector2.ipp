
namespace geometrycentral {

inline Vector2 Vector2::operator+(const Vector2& v) const { return Vector2{x + v.x, y + v.y}; }

inline Vector2 Vector2::operator-(const Vector2& v) const { return Vector2{x - v.x, y - v.y}; }

inline Vector2 Vector2::operator*(const Vector2& v) const { return Vector2{x * v.x - y * v.y, x * v.y + y * v.x}; }

inline Vector2 Vector2::operator/(const Vector2& v) const {
  double denom = v.x * v.x + v.y * v.y;
  return Vector2{x * v.x + y * v.y, y * v.x - x * v.y} / denom;
}

inline Vector2 Vector2::operator*(double s) const { return Vector2{x * s, y * s}; }

inline Vector2 Vector2::operator/(double s) const {
  const double r = 1. / s;
  return Vector2{x * r, y * r};
}

inline const Vector2 Vector2::operator-() const { return Vector2{-x, -y}; }

template <typename T>
inline Vector2 operator*(const T s, const Vector2& v) {
  return Vector2{s * v.x, s * v.y};
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

inline Vector2& Vector2::operator*=(const Vector2& other) {
  Vector2 tmp = *this * other;
  *this = tmp;
  return *this;
}

inline Vector2& Vector2::operator/=(const Vector2& other) {
  Vector2 tmp = *this / other;
  *this = tmp;
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

inline bool Vector2::operator==(const Vector2& other) const { return x == other.x && y == other.y; }

inline bool Vector2::operator!=(const Vector2& other) const { return !(*this == other); }


inline Vector2::operator std::complex<double>() const { return std::complex<double>{x, y}; }

inline Vector2 Vector2::normalize() const {
  double r = 1. / std::sqrt(x * x + y * y);
  return *this * r;
}

inline Vector2 Vector2::unit() const { return normalize(); }


inline Vector2 Vector2::normalizeCutoff(double mag) const {
  double len = std::sqrt(x * x + y * y);
  if (len <= mag) len = 1.;
  double r = 1. / len;
  return *this * r;
}

inline Vector2 normalize(const Vector2& v) { return v.normalize(); }

inline Vector2 unit(const Vector2& v) { return normalize(v); }

inline Vector2 normalizeCutoff(const Vector2& v, double mag) { return v.normalizeCutoff(mag); }

inline Vector2 Vector2::rotate(double theta) const {
  double cosTh = std::cos(theta);
  double sinTh = std::sin(theta);
  return Vector2{cosTh * x + sinTh * y, -sinTh * x + cosTh * y};
}

inline Vector2 Vector2::rotate90() const { return Vector2{-y, x}; }

inline Vector2 Vector2::pow(double p) const {
  std::complex<double> c{x, y};
  c = std::pow(c, p);
  return Vector2{c.real(), c.imag()};
}

inline Vector2 Vector2::pow(Vector2 p) const {
  std::complex<double> c{x, y};
  std::complex<double> pc{p.x, p.y};
  c = std::pow(c, pc);
  return Vector2{c.real(), c.imag()};
}

inline Vector2 Vector2::conj() const { return Vector2{x, -y}; }

inline Vector2 Vector2::inv() const { return Vector2{1., 0.} / *this; }

inline double Vector2::arg() const { return std::atan2(y, x); }
inline double arg(const Vector2& v) { return std::atan2(v.y, v.x); }

inline double Vector2::norm() const { return std::sqrt(x * x + y * y); }
inline double norm(const Vector2& v) { return std::sqrt(v.x * v.x + v.y * v.y); }

inline double Vector2::norm2() const { return x * x + y * y; }
inline double norm2(const Vector2& v) { return v.x * v.x + v.y * v.y; }


inline double dot(const Vector2& u, const Vector2& v) { return u.x * v.x + u.y * v.y; }

inline double angle(const Vector2& u, const Vector2& v) {
  return std::acos(std::fmax(-1., std::fmin(1., dot(unit(u), unit(v)))));
}
inline double orientedAngle(const Vector2& u, const Vector2& v) {
  Vector2 uHat = unit(u);
  Vector2 vHat = unit(v);
  return std::atan2(cross(uHat, vHat), dot(uHat, vHat));
}

inline double cross(const Vector2& u, const Vector2& v) { return u.x * v.y - u.y * v.x; }
inline Vector3 cross3(const Vector2& u, const Vector2& v) { return Vector3{0., 0., u.x * v.y - u.y * v.x}; }
inline Vector2 clamp(const Vector2& val, const Vector2& low, const Vector2& high) {
  Vector2 rVal;
  for (int i = 0; i < 2; i++) {
    rVal[i] = clamp(val[i], low[i], high[i]);
  }
  return rVal;
}

inline bool Vector2::isFinite() const { return std::isfinite(x) && std::isfinite(y); }
inline bool isfinite(const Vector2& v) { return v.isFinite(); }
inline bool Vector2::isDefined() const { return (!std::isnan(x)) && (!std::isnan(y)); }
inline bool isDefined(const Vector2& v) { return v.isDefined(); }

inline Vector2 componentwiseMin(const Vector2& u, const Vector2& v) {
  return Vector2{std::fmin(u.x, v.x), std::fmin(u.y, v.y)};
}
inline Vector2 componentwiseMax(const Vector2& u, const Vector2& v) {
  return Vector2{std::fmax(u.x, v.x), std::fmax(u.y, v.y)};
}
inline double sum(const Vector2& u) { return u.x + u.y; }

inline std::ostream& operator<<(std::ostream& output, const Vector2& v) {
  output << "<" << v.x << ", " << v.y << ">";
  return output;
}

inline std::istream& operator>>(std::istream& input, Vector2& v) {
  double x, y;
  input >> x >> y;
  v = Vector2{x, y};
  return input;
}

} // namespace geometrycentral

namespace std {
inline std::size_t std::hash<geometrycentral::Vector2>::operator()(const geometrycentral::Vector2& v) const {
  return std::hash<double>{}(v.x) ^ (std::hash<double>{}(v.y) + (std::hash<double>{}(v.y) << 2)) ^
         (std::hash<double>{}(v.x) + (std::hash<double>{}(v.x) << 4));
}

inline std::string to_string(geometrycentral::Vector2 vec) {
  ostringstream output;
  output << vec;
  return output.str();
}


} // namespace std
