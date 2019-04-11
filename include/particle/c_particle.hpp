#ifndef _C_PARTICLE_HPP_
#define _C_PARTICLE_HPP_

#include "particle/particle_expression.hpp"
#include "apt/foreach.hpp"
#include "apt/array.hpp"

namespace particle {
  // for communication
  template < typename T, int DPtc, typename state_t >
  class cParticle : public PtcExpression<cParticle<T,DPtc,state_t>, apt::array<T,DPtc>, state_t> {
  private:
    apt::array<T,DPtc> _q {};
    apt::array<T,DPtc> _p {};
    state_t _s {};

  public:
    static constexpr int NDim = DPtc;
    using vec_type = apt::array<T,DPtc>;
    using state_type = state_t;

    constexpr vec_type& q() noexcept { return _q; }
    constexpr const vec_type& q() const noexcept { return _q; }

    constexpr vec_type& p() noexcept { return _p; }
    constexpr const vec_type& p() const noexcept { return _p; }

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

    template < typename E >
    cParticle& operator= ( const PtcExpression<E>& ptc ) noexcept {
      apt::foreach<0,DPtc>([](auto& a, auto& b){ a = b;}, q(), ptc.q() );
      apt::foreach<0,DPtc>([](auto& a, auto& b){ a = b;}, p(), ptc.p() );
      _s = ptc.state();

      return *this;
    }

    template < typename E >
    cParticle& operator= ( PtcExpression<E>&& ptc ) noexcept {
      apt::foreach<0,DPtc>([](auto& a, auto& b){ std::swap(a,b);}, q(), ptc.q() );
      apt::foreach<0,DPtc>([](auto& a, auto& b){ std::swap(a,b);}, p(), ptc.p() );
      std::swap( _s, ptc.state() );
      ptc.set(flag::empty);

      return *this;
    }

    constexpr void swap( cParticle& other ) noexcept {
      apt::foreach<0,DPtc>([](auto& a, auto& b){ std::swap(a,b);}, q(), other.q() );
      apt::foreach<0,DPtc>([](auto& a, auto& b){ std::swap(a,b);}, p(), other.p() );
      std::swap( _s, other.state() );
    }

  };
}

#endif
