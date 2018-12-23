#ifndef  _PHYS_VECTOR_HPP_
#define  _PHYS_VECTOR_HPP_

#include <cmath>
#include <tuple>
#include <utility>

#include "traits.hpp"

// TODOL index_sequence and fold expression may help the generalization to arbitrary Dim
template < typename T >
struct Vec3 : public std::tuple<T,T,T> {
  static constexpr int Dim = 3;

  template < typename... Args >
  Vec3( Args&&... args ) : std::tuple<T, T, T>( std::forward<Args>(args)... ) {}

  constexpr Vec3& operator+= ( const Vec3<const T&>& b ) {
    std::get<0>(*this) += std::get<0>(b);
    std::get<1>(*this) += std::get<1>(b);
    std::get<2>(*this) += std::get<2>(b);
    return *this;
  }

  constexpr Vec3& operator-= ( const Vec3<const T&>& b ) {
    std::get<0>(*this) -= std::get<0>(b);
    std::get<1>(*this) -= std::get<1>(b);
    std::get<2>(*this) -= std::get<2>(b);
    return *this;
  }

  constexpr Vec3& operator*= ( const typename std::remove_reference<T>::type b) {
    std::get<0>(*this) *= b;
    std::get<1>(*this) *= b;
    std::get<2>(*this) *= b;
    return *this;
  }

  constexpr Vec3& operator/= ( const typename std::remove_reference<T>::type b) {
    std::get<0>(*this) /= b;
    std::get<1>(*this) /= b;
    std::get<2>(*this) /= b;
    return *this;
  }

  // constexpr decltype(auto) operator[] ( const int i ) {
  // }
};

template < typename T >
constexpr auto operator+ ( const Vec3<const T&>& a, const Vec3<const T&>& b ) {
  return Vec3<T>( std::get<0>(a) + std::get<0>(b),
                  std::get<1>(a) + std::get<1>(b),
                  std::get<2>(a) + std::get<2>(b) );
}

template < typename T >
constexpr auto operator- ( const Vec3<const T&>& a, const Vec3<const T&>& b ) {
  return Vec3<T>( std::get<0>(a) - std::get<0>(b),
                  std::get<1>(a) - std::get<1>(b),
                  std::get<2>(a) - std::get<2>(b) );
}

template < typename T >
constexpr auto operator* ( const Vec3<const T&>& a, const T b ) {
  return Vec3<T>( std::get<0>(a) * b,
                  std::get<1>(a) * b,
                  std::get<2>(a) * b );
}

template < typename T >
constexpr auto operator* ( const T b, const Vec3<const T&>& a ) {
  return a * b;
}

template < typename T >
constexpr auto operator/ ( const Vec3<const T&>& a, const T b ) {
  return Vec3<T>( std::get<0>(a) / b,
                  std::get<1>(a) / b,
                  std::get<2>(a) / b );
}


struct vecop {
private:
  template < int I = 0, typename Vector >
  static constexpr auto sum_square ( const Vector& vec ) {
    static_assert( I < Vector::Dim, "Invalid I" );
    if constexpr ( I == Vector::Dim - 1 ) return std::get<I>(vec) * std::get<I>(vec);
    else return std::get<I>(vec) * std::get<I>(vec) + sum_square< I+1 >(vec) ;
  };

  public:

  template < typename Vector >
  static constexpr decltype(auto) abs_sq( const Vector& vec ) {
    return sum_square<>(vec);
  }

  template < typename Vector >
  static constexpr decltype(auto) abs( const Vector& vec ) {
    return std::sqrt( abs_sq( vec ) );
  }
};

#endif
