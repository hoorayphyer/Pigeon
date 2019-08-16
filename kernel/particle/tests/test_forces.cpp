#include "testfw/testfw.hpp"
#include "particle/forces_impl.hpp"
#include "particle/particle.hpp"

using namespace particle;

namespace old_vay {
  template <typename T>
  struct Vec3 {
    T x, y, z;

    typedef Vec3<T> self_type;
    typedef T elem_type;

    Vec3() : x(static_cast<T>(0)), y(static_cast<T>(0)), z(static_cast<T>(0)) {}
    Vec3(T xi, T yi, T zi) : x(xi), y(yi), z(zi) {}

    // conversion from apt::Vec
    Vec3(const apt::Vec<T,3>& vec)
      : Vec3(vec[0], vec[1], vec[2]) {}

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

  // lambda = charge_to_mass * dt;
  inline Vec3<double> lorentz( double lambda, const Vec3<double>& p, const Vec3<double>& EVector, const Vec3<double>& BVector ) {
    lambda /=  2.0; // measured in units of (e/m) * R_* / c;
    Vec3<double> u_halfstep = p + EVector * lambda + p.cross( BVector ) * ( lambda / std::sqrt( 1.0 + p.dot(p) ) );
    Vec3<double> upr = u_halfstep + EVector * lambda;
    Vec3<double> tau = BVector * lambda;
    // store some repeatedly used intermediate results
    double tt = tau.dot( tau );
    double ut = upr.dot( tau );

    double sigma = 1.0 + upr.dot(upr) - tt;
    // inv_gamma2 means ( 1 / gamma^(i+1) ) ^2
    double inv_gamma2 =  2.0 / ( sigma + std::sqrt( sigma * sigma + 4.0 * ( tt + ut * ut ) ) );
    double s = 1.0 / ( 1.0 + inv_gamma2 * tt );
    return ( upr + tau * ( ut * inv_gamma2 ) + upr.cross( tau ) * std::sqrt(inv_gamma2) ) * s;
  }
}

SCENARIO("Test against old Vay Pusher", "[particle]") {
  aio::unif_real<double> unif(-1.0 , 1.0);
  const int N = 1000000;
  for ( int i = 0; i < N; ++i ) {
    apt::Vec<double, 3> E( unif(), unif(), unif() );
    apt::Vec<double, 3> B( unif(), unif(), unif() );
    apt::Vec<double, 3> p( unif(), unif(), unif() );
    p *= 10.0;

    double dt = 0.01;
    double w_gyro_unitB = 10000.0;
    double q_over_m = 1.0;

    auto p_old = old_vay::lorentz( q_over_m * dt * w_gyro_unitB, p, E, B );

    Particle<double,aio::Specs> ptc;
    ptc.p() = p;
    force::lorentz( ptc, dt, E, B, q_over_m * w_gyro_unitB );
    for ( int i = 0; i < 3; ++i ) {
      REQUIRE( p_old[i] == Approx(ptc.p()[i]) );
    }
  }

}

SCENARIO("Test gyration resolved Lorentz force", "[particle]") {
  SECTION("E = 0, uniform B in the z direction") {
    aio::unif_real<double> unif(-1.0 , 1.0);
    // unif.seed(0);
    const int N = 100000;
    apt::Vec<double, 3> E;
    for ( int i = 0; i < N; ++i ) {
      apt::Vec<double, 3> B( 0.0, 0.0, unif() );
      apt::Vec<double, 3> p( unif(), unif(), unif() );

      // B = {0.0,0.0,1.0};
      // p = {1.0,0.0,0.0};

      const double gamma = 10.0 * std::abs(unif()) + 1.0;
      B /= apt::abs(B);
      p *= std::sqrt( (gamma * gamma - 1.0)  / apt::sqabs(p) );

      double dt = 0.0001;
      double w_B = 10000.0; // qB / m_q c

      double angle = w_B * dt / gamma;
      if ( B[2] < 0 ) angle *= -1;
      double p_thy_0 = p[0] * std::cos(angle) - p[1] * std::sin(angle);
      double p_thy_1 = p[0] * std::sin(angle) + p[1] * std::cos(angle);
      double p_thy_2 = p[2];

      Particle<double,aio::Specs> ptc;
      ptc.p() = p;
      force::lorentz_exact( ptc, dt, E, B, w_B );
      REQUIRE( ptc.p()[0] == Approx(p_thy_0) );
      REQUIRE( ptc.p()[1] == Approx(p_thy_1) );
      REQUIRE( ptc.p()[2] == Approx(p_thy_2) );
    }
  }

  SECTION("Test energy conservation when E // B") {
    aio::unif_real<double> unif(-1.0 , 1.0);
    // unif.seed(0);
    const int N = 100000;
    for ( int i = 0; i < N; ++i ) {
      apt::Vec<double, 3> B( unif(), unif(), unif() );
      apt::Vec<double, 3> p( unif(), unif(), unif() );

      B /= apt::abs(B);

      apt::Vec<double, 3> E = 2 * unif() * B; // E can be larger than B

      const double gamma = 10.0 * std::abs(unif()) + 1;
      p *= std::sqrt( (gamma * gamma - 1.0)  / apt::sqabs(p) );

      double dt = 0.0002;
      double w_B = 720000.0;

      // NOTE this formula works only when E // B
      const double gamma_new = sqrt( gamma * gamma + 2.0 * apt::dot(p,E) * dt * w_B + apt::sqabs(E) * dt * dt * w_B * w_B );

      // auto p_old = old_vay::lorentz( w_B * dt, p, E, B );

      Particle<double,aio::Specs> ptc;
      ptc.p() = p;
      force::lorentz_exact( ptc, dt, E, B, w_B );

      // CHECK( 1.0 + p_old[0]*p_old[0] + p_old[1]*p_old[1] + p_old[2]*p_old[2] == Approx( gamma_new * gamma_new ) );
      REQUIRE( 1.0 + apt::sqabs(ptc.p()) == Approx(gamma_new * gamma_new) );
    }
  }

  SECTION("Test energy conservation when E \perp B and E < B, and p // ExB") {
    aio::unif_real<double> unif(-1.0 , 1.0);
    unif.seed(0);
    const int N = 1;
    for ( int i = 0; i < N; ++i ) {
      apt::Vec<double, 3> E( unif(), unif(), unif() );
      apt::Vec<double, 3> B( unif(), unif(), unif() );
      apt::Vec<double, 3> p( unif(), unif(), unif() );

      B /= apt::abs(B);
      E = E - B * apt::dot(E,B);
      E /= apt::abs(E);
      p = p - E * apt::dot(E,p);
      p = p - B * apt::dot(B,p);
      E *= 0.5;

      const double gamma = 2.0;
      p *= std::sqrt( (gamma * gamma - 1.0)  / apt::sqabs(p) );

      double dt = 0.0002;
      double w_B = 7200.0;

      // auto p_old = old_vay::lorentz( w_B * dt, p, E, B );

      Particle<double,aio::Specs> ptc;
      ptc.p() = p;
      force::lorentz_exact( ptc, dt, E, B, w_B );

      // CHECK( 1.0 + p_old[0]*p_old[0] + p_old[1]*p_old[1] + p_old[2]*p_old[2] == Approx( gamma * gamma ) );
      // TODO test wrong, there can still be a little acceleration
      // REQUIRE( 1.0 + apt::sqabs(ptc.p()) == Approx(gamma * gamma) );
    }
  }

  // SECTION("Random") {
  //   aio::unif_real<double> unif(-1.0 , 1.0);
  //   unif.seed(0);
  //   const int N = 1;
  //   for ( int i = 0; i < N; ++i ) {
  //     apt::Vec<double, 3> E( unif(), unif(), unif() );
  //     apt::Vec<double, 3> B( unif(), unif(), unif() );
  //     apt::Vec<double, 3> p( unif(), unif(), unif() );
  //     E = { 0.00760478, 0.00280827, -0.00196783 };
  //     B = { -0.00251577, 0.0078043, 0.00141511 };
  //     p = { 2.35757, -3.76661, 5.40618 };

  //     double dt = 0.003;
  //     double w_B = 3750.0;

  //     auto p_old = old_vay::lorentz( w_B * dt, p, E, B );

  //     Particle<double,aio::Specs> ptc;
  //     ptc.p() = p;
  //     force::lorentz_exact( ptc, dt, E, B, w_B );
  //     CHECK( ptc.p()[0] == 0 );
  //     CHECK( ptc.p()[1] == 0 );
  //     CHECK( ptc.p()[2] == 0 );
  //   }
  // }
}
