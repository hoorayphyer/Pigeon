#include "all_in_one.hpp"
#include "pic.hpp"
#include "particle/forces.cpp"
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

SCENARIO("Test against old Vay Pusher") {
  aio::unif_real<double> unif(-1.0 , 1.0);
  const int N = 1000000;
  for ( int i = 0; i < N; ++i ) {
    apt::Vec<double, 3> E( unif(), unif(), unif() );
    apt::Vec<double, 3> B( unif(), unif(), unif() );
    apt::Vec<double, 3> p( unif(), unif(), unif() );
    E *= 1000.0;
    B *= 1000.0;
    p *= 10.0;

    double dt = 0.01;
    double q_over_m = 1.0;

    auto p_old = old_vay::lorentz( q_over_m * dt, p, E, B );

    Particle<double,Specs> ptc;
    ptc.p() = p;
    force::lorentz( ptc, dt, E, B, q_over_m );
    for ( int i = 0; i < 3; ++i ) {
      REQUIRE( p_old[i] == Approx(ptc.p()[i]) );
    }
  }

}

