#ifndef _C_PARTICLE_HPP_
#define _C_PARTICLE_HPP_

#include "particle/particle_expression.hpp"
#include "apt/foreach.hpp"

namespace particle {
  // cParticle uses C-native types for the sake of interfacing with MPI
  template < typename T, int DPtc, typename state_t >
  class cParticle : public PtcExpression<cParticle<T,DPtc,state_t>, T[DPtc], state_t> {
  private:
    std::array<T,DPtc> _q;
    std::array<T,DPtc> _p;
    state_t _s;

  public:
    static constexpr int Dim = DPtc;
    using vec_type = std::array<T,DPtc>;
    using state_type = state_t;

    constexpr auto& q() noexcept { return _q; }
    constexpr const auto& q() const noexcept { return _q; }

    constexpr auto& p() noexcept { return _p; }
    constexpr const auto& p() const noexcept { return _p; }

    constexpr auto& state() noexcept { return _s; }
    constexpr const auto& state() const noexcept { return _s; }

    cParticle() = default;

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
