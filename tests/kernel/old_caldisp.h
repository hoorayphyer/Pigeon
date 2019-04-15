using std::sin;
using std::cos;
using std::exp;
using std::acos;
using std::atan;

template <typename T>
struct Vec3 {
  T x, y, z;

  typedef Vec3<T> self_type;
  typedef T elem_type;

  Vec3() : x(static_cast<T>(0)), y(static_cast<T>(0)), z(static_cast<T>(0)) {}
  Vec3(T xi, T yi, T zi) : x(xi), y(yi), z(zi) {}

  Vec3(const self_type& other) = default;
  Vec3(self_type&& other) = default;

  inline T& operator[](int idx) {
    if (idx == 0)
      return x;
    else if (idx == 1)
      return y;
    // else if (idx == 2)
    else
      return z;
    // else
    //   throw std::out_of_range("Index out of bound!");
  }

  inline const T& operator[](int idx) const {
    if (idx == 0)
      return x;
    else if (idx == 1)
      return y;
    // else if (idx == 2)
    else
      return z;
    // else
      // throw std::out_of_range("Index out of bound!");
  }

  inline self_type& operator= (const self_type& other) {
    x = other.x;
    y = other.y;
    z = other.z;
    return *this;
  }

  inline bool operator== (const self_type& other) const {
    return (x == other.x && y == other.y && z == other.z);
  }

  inline bool operator!= (const self_type& other) const {
    return (x != other.x || y != other.y || z != other.z);
  }

  inline self_type& operator+=(const self_type& other) {
    x += other.x;
    y += other.y;
    z += other.z;
    return (*this);
  }

  inline self_type& operator-=(const self_type& other) {
    x -= other.x;
    y -= other.y;
    z -= other.z;
    return (*this);
  }

  inline self_type& operator*=( T n) {
    x *= n;
    y *= n;
    z *= n;
    return (*this);
  }

  inline self_type& operator/=( T n) {
    x /= n;
    y /= n;
    z /= n;
    return (*this);
  }

  inline self_type operator+(const self_type& other) const {
    self_type tmp = *this;
    tmp += other;
    return tmp;
  }

  inline self_type operator-(const self_type& other) const {
    self_type tmp = *this;
    tmp -= other;
    return tmp;
  }

  inline self_type operator* (T n) const {
    self_type tmp { x * n, y * n, z * n };
    return tmp;
  }

  inline self_type operator/ (T n) const {
    self_type tmp { x / n, y / n, z / n };
    return tmp;
  }

  inline T dot(const self_type& other) const {
    return (x * other.x + y * other.y + z * other.z);
  }

  inline self_type cross(const self_type& other) const {
    return Vec3<T>(y * other.z - z * other.y, z * other.x - x * other.z,
                   x * other.y - y * other.x);
  }
};

template < typename T >
struct ScalesLogSpherical {
  const T PI = 3.141592653589793238462643383279502884197169399375105820974944592307816406286L;

  // Transform a vector to Cartesian coordinates
  inline void VectorToCartesian(Vec3<T>& v,
                                const Vec3<T>& pos) const {
    T v1n = v.x, v2n = v.y, v3n = v.z;
    // Now these become vx, vy, vz
    v.x = v1n * sin(pos.y) * cos(pos.z) + v2n * cos(pos.y) * cos(pos.z) -
          v3n * sin(pos.z);
    v.y = v1n * sin(pos.y) * sin(pos.z) + v2n * cos(pos.y) * sin(pos.z) +
          v3n * cos(pos.z);
    v.z = v1n * cos(pos.y) - v2n * sin(pos.y);
  }

  // Transform a vector from Cartesian coordinates
  inline void VectorFromCartesian(Vec3<T>& v,
                                  const Vec3<T>& pos) const {
    T v1n = v.x, v2n = v.y, v3n = v.z;
    // These become vr, vtheta, vphi
    v.x = v1n * sin(pos.y) * cos(pos.z) + v2n * sin(pos.y) * sin(pos.z) +
          v3n * cos(pos.y);
    v.y = v1n * cos(pos.y) * cos(pos.z) + v2n * cos(pos.y) * sin(pos.z) -
          v3n * sin(pos.y);
    v.z = -v1n * sin(pos.z) + v2n * cos(pos.z);
  }

  inline void PosToCartesian(Vec3<T>& pos) const {
    T x1n = pos.x, x2n = pos.y, x3n = pos.z;
    // These become x, y, z
    pos.x = exp(x1n) * sin(x2n) * cos(x3n);
    pos.y = exp(x1n) * sin(x2n) * sin(x3n);
    pos.z = exp(x1n) * cos(x2n);
  }

  inline void PosFromCartesian(Vec3<T>& pos) const {
    T x1n = pos.x, x2n = pos.y, x3n = pos.z;
    // These become r, theta, phi
    pos.x = sqrt(x1n * x1n + x2n * x2n + x3n * x3n);
    pos.y = acos(x3n / pos.x);  // acos(z / r)
    // TODO: Check correctness of this statement!
    pos.z = atan(x2n / x1n) + (x1n < 0) * PI;  // atan(y / x)
    pos.x = log(pos.x);
  }
};


template < typename T >
inline void logsph_geodesic_move_old ( Vec3<T>& x, Vec3<T>& v, T dt ) {
  const T PI = 3.141592653589793238462643383279502884197169399375105820974944592307816406286L;
  static ScalesLogSpherical<T> _scales;
  Vec3<T> x_new(x);
  _scales . PosToCartesian(x_new);

  //FIXME: should use a half-step position for p conversion. How?
  _scales . VectorToCartesian( v, x);
  x_new += v * dt;

  _scales . PosFromCartesian(x_new);

  //FIXME: This only works in axisymmetric coordinates, how to generalize?
  if (x_new.z - x.z > PI) x_new.z -= 2.0 * PI;
  if (x_new.z - x.z < -PI) x_new.z += 2.0 * PI;

  for ( int i = 0; i < 3; ++i ) {
    x[i] = x_new[i];
  }

  // update momentum components under the new local coordinates.
  _scales . VectorFromCartesian( v, x_new );
}

inline void OLD_geodesic_move(Vec<double,3>& x, Vec<double,3>& v, const double& dt ) noexcept {
  Vec3<double> x_o, v_o;
  // convert from the new Vec<T,3> to old Vec3<T>
  for ( int i = 0; i < 3; ++i ) {
    x_o[i] = x[i];
    v_o[i] = v[i];
  }

  logsph_geodesic_move_old( x_o, v_o, dt );

  // convert back
  for ( int i = 0; i < 3; ++i ) {
    x[i] = x_o[i];
    v[i] = v_o[i];
  }

}
