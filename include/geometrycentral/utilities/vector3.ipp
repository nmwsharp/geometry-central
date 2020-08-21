namespace geometrycentral {

inline Vector3 Vector3::operator+(const Vector3& v) const { return Vector3{x + v.x, y + v.y, z + v.z}; }

inline Vector3 Vector3::operator-(const Vector3& v) const { return Vector3{x - v.x, y - v.y, z - v.z}; }

inline Vector3 Vector3::operator*(double s) const { return Vector3{x * s, y * s, z * s}; }

inline Vector3 Vector3::operator/(double s) const {
  const double r = 1. / s;
  return Vector3{x * r, y * r, z * r};
}

inline const Vector3 Vector3::operator-() const { return Vector3{-x, -y, -z}; }

template <typename T>
inline Vector3 operator*(const T s, const Vector3& v) {
  return Vector3{s * v.x, s * v.y, s * v.z};
}

inline Vector3& Vector3::operator+=(const Vector3& other) {
  x += other.x;
  y += other.y;
  z += other.z;
  return *this;
}

inline Vector3& Vector3::operator-=(const Vector3& other) {
  x -= other.x;
  y -= other.y;
  z -= other.z;
  return *this;
}

inline Vector3& Vector3::operator*=(const double& s) {
  x *= s;
  y *= s;
  z *= s;
  return *this;
}

inline Vector3& Vector3::operator/=(const double& s) {
  x /= s;
  y /= s;
  z /= s;
  return *this;
}

inline bool Vector3::operator==(const Vector3& other) const { return x == other.x && y == other.y && z == other.z; }

inline bool Vector3::operator!=(const Vector3& other) const { return !(*this == other); }

inline double Vector3::norm() const { return std::sqrt(x * x + y * y + z * z); }
inline double norm(const Vector3& v) { return v.norm(); }

inline double Vector3::norm2() const { return x * x + y * y + z * z; }
inline double norm2(const Vector3& v) { return v.norm2(); }

inline Vector3 normalize(const Vector3& v) { return v.normalize(); }

inline Vector3 normalizeCutoff(const Vector3& v, double mag) { return v.normalizeCutoff(mag); }

inline Vector3 unit(const Vector3& v) { return normalize(v); }

inline Vector3 cross(const Vector3& u, const Vector3& v) {
  double x = u.y * v.z - u.z * v.y;
  double y = u.z * v.x - u.x * v.z;
  double z = u.x * v.y - u.y * v.x;
  return Vector3{x, y, z};
}

inline double dot(const Vector3& u, const Vector3& v) { return u.x * v.x + u.y * v.y + u.z * v.z; }
inline double sum(const Vector3& u) { return u.x + u.y + u.z; }

inline double angle(const Vector3& u, const Vector3& v) {
  return std::acos(std::fmax(-1., std::fmin(1., dot(unit(u), unit(v)))));
}

inline double angleInPlane(const Vector3& u, const Vector3& v, const Vector3& normal) {
  // Put u in plane with the normal
  Vector3 N = unit(normal);
  Vector3 uPlane = unit(u - dot(u, N) * N);
  Vector3 basisY = unit(cross(normal, uPlane));

  double xComp = dot(v, uPlane);
  double yComp = dot(v, basisY);

  return ::std::atan2(yComp, xComp);
}

inline bool Vector3::isFinite() const { return ::std::isfinite(x) && ::std::isfinite(y) && ::std::isfinite(z); }
inline bool isfinite(const Vector3& v) { return v.isFinite(); }

inline bool Vector3::isDefined() const { return (!::std::isnan(x)) && (!::std::isnan(y)) && (!::std::isnan(z)); }
inline bool isDefined(const Vector3& v) { return v.isDefined(); }

inline Vector3 clamp(const Vector3& val, const Vector3& low, const Vector3& high) {
  Vector3 rVal;
  for (int i = 0; i < 3; i++) {
    rVal[i] = clamp(val[i], low[i], high[i]);
  }
  return rVal;
}

inline Vector3 componentwiseMin(const Vector3& u, const Vector3& v) {
  return Vector3{std::fmin(u.x, v.x), std::fmin(u.y, v.y), std::fmin(u.z, v.z)};
}

inline Vector3 componentwiseMax(const Vector3& u, const Vector3& v) {
  return Vector3{std::fmax(u.x, v.x), std::fmax(u.y, v.y), std::fmax(u.z, v.z)};
}

inline Vector3 Vector3::rotateAround(Vector3 axis, double theta) const {
  Vector3 thisV = {x, y, z};
  Vector3 axisN = axis.normalize();
  Vector3 parallelComp = axisN * dot(thisV, axisN);
  Vector3 tangentComp = thisV - parallelComp;

  if (tangentComp.norm2() > 0.0) {
    Vector3 basisX = tangentComp.normalize();
    Vector3 basisY = cross(axisN, basisX);

    double tangentMag = tangentComp.norm();

    Vector3 rotatedV = tangentMag * (std::cos(theta) * basisX + std::sin(theta) * basisY);
    return rotatedV + parallelComp;
  } else {
    return parallelComp;
  }
}

inline Vector3 Vector3::normalize() const {
  double r = 1. / std::sqrt(x * x + y * y + z * z);
  return *this * r;
}

inline Vector3 Vector3::normalizeCutoff(double mag) const {
  double len = std::sqrt(x * x + y * y + z * z);
  if (len <= mag) len = 1.;
  double r = 1. / len;
  return *this * r;
}

inline Vector3 Vector3::unit() const { return normalize(); }

inline Vector3 Vector3::removeComponent(const Vector3& unitDir) const { return *this - unitDir * dot(unitDir, *this); }

inline std::array<Vector3, 2> Vector3::buildTangentBasis() const {
  Vector3 unitDir = normalize();
  Vector3 testVec{1., 0., 0.};
  if (std::fabs(dot(testVec, unitDir)) > 0.9) {
    testVec = Vector3{0., 1., 0.};
  }

  Vector3 basisX = cross(testVec, unitDir).normalize();
  Vector3 basisY = cross(unitDir, basisX).normalize();

  return std::array<Vector3, 2>{basisX, basisY};
};

inline std::ostream& operator<<(std::ostream& output, const Vector3& v) {
  output << "<" << v.x << ", " << v.y << ", " << v.z << ">";
  return output;
}

inline std::istream& operator>>(std::istream& input, Vector3& v) {
  double x, y, z;
  input >> x >> y >> z;
  v = Vector3{x, y, z};
  return input;
}

} // namespace geometrycentral

namespace std {
inline std::size_t std::hash<geometrycentral::Vector3>::operator()(const geometrycentral::Vector3& v) const {
  return std::hash<double>{}(v.x) ^ (std::hash<double>{}(v.y) + (std::hash<double>{}(v.y) << 2)) ^
         (std::hash<double>{}(v.z) + (std::hash<double>{}(v.z) << 4));
}

inline std::string to_string(geometrycentral::Vector3 vec) {
  ostringstream output;
  output << vec;
  return output.str();
}


} // namespace std
