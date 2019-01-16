#ifndef _PARTICLE_HPP_
#define _PARTICLE_HPP_

#include "types.hpp"
#include <array>
#include <vector>
#include "vector.hpp"

template < class PtcArray >
struct Particle {
  static constexpr auto Dim = PtcArray::Dptc;

  PtcArray& arr;
  const int i;
};

template < std::size_t Dim_Ptc >
struct ParticleArray {
  std::array<std::vector<T>, Dim_Ptc> q;
  std::array<std::vector<T>, Dim_Ptc> p;
  std::vector<encoded_bits_t> state;

  static constexpr auto DPtc = Dim_Ptc;

  inline auto operator[] ( int i ) noexcept {
    return Particle<ParticleArray>{ *this, i };
  }

  inline auto operator[] ( int i ) const noexcept {
    return Particle<const ParticleArray>{ *this, i };
  }
};

namespace mem {

  template <typename Ptc >
  inline auto q( Ptc& ptc ) noexcept {
    return vec::per_dim::tie<Ptc::Dim>
      ( [i=ptc.i]( auto&& x ) { return x[i]; }, ptc.arr.q );
  }

  template <typename Ptc >
  inline auto p( Ptc& ptc ) noexcept {
    return vec::per_dim::tie<Ptc::Dim>
      ( [i=ptc.i]( auto&& x ) { return x[i]; }, ptc.arr.p );
  }

  template < std::size_t D, typename Ptc >
  inline auto q( Ptc& ptc ) noexcept {
    static_assert( D <= Ptc::Dim );
    return vec::per_dim::tie<D>
      ( [i=ptc.i]( auto&& x ) { return x[i]; }, ptc.arr.q );
  }
}

#endif
