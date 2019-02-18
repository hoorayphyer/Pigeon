#ifndef _C_PARTICLE_HPP_
#define _C_PARTICLE_HPP_

#include "particle/particle_expression.hpp"
#include "apt/foreach.hpp"

namespace particle {
  // cParticle uses C-native types for the sake of interfacing with MPI
  template < typename T, int DPtc, typename state_t >
  class cParticle : public PtcExpression<cParticle<T,DPtc,state_t>, T[DPtc], state_t> { // TODO check T[DPtc]
  private:
    T _q [DPtc];
    T _p [DPtc];
    state_t _s;

  public:
    static constexpr int Dim = DPtc;
    using vec_type = T[DPtc];
    using state_type = state_t;
    // TODO double check this: std::get on C-style array is supported with bound checking. At least, std::begin and std::end do.
    constexpr auto& q() noexcept { return T(&_q)[DPtc]; }
    constexpr const auto& q() const noexcept { return (const T)(&_q)[DPtc]; }

    constexpr auto& p() noexcept { return T(&_p)[DPtc]; }
    constexpr const auto& p() const noexcept { return (const T)(&_p)[DPtc]; }

    constexpr auto& state() noexcept { return _s; }
    constexpr const auto& state() const noexcept { return _s; }

    template < typename E >
    cParticle( PtcExpression<E>&& ptc ) noexcept {
      apt::foreach<0,DPtc>([](auto& a, auto& b){ std::swap(a,b);}, q(), ptc.q() );
      apt::foreach<0,DPtc>([](auto& a, auto& b){ std::swap(a,b);}, p(), ptc.p() );
      std::swap( _s, ptc.state() );
      ptc.set(flag::empty);
    }

  };
}

#endif
