#ifndef _VIRTUAL_PARTICLE_HPP_
#define _VIRTUAL_PARTICLE_HPP_

#include "particle/particle_expression.hpp"
#include "apt/virtual_vec.hpp"

namespace particle {
  template < typename T, int DPtc, typename state_t >
  struct vParticle : public PtcExpression< vParticle<T, DPtc, state_t>, apt::vVec<T,DPtc>, state_t > {
  private:
    static_assert( 8 * sizeof( state_t) >= 64 );

    apt::vVec<T, DPtc> _q;
    apt::vVec<T, DPtc> _p;
    state_t& _state;

  public:
    static constexpr int NDim = DPtc;
    using vec_type = apt::vVec<T, DPtc>;
    using state_type = state_t;

    constexpr auto& q() noexcept { return _q; }
    constexpr const auto& q() const noexcept { return _q; }

    constexpr auto& p() noexcept { return _p; }
    constexpr const auto& p() const noexcept { return _p; }

    constexpr auto& state() noexcept { return _state; }
    constexpr const auto& state() const noexcept { return _state; }

    template < typename Q, typename P >
    vParticle( Q&& q, P&& p, state_t& state ) noexcept
      : _q(std::forward<Q>(q)),
        _p(std::forward<P>(p)),
        _state(state) {}

    vParticle( vParticle&& ptc ) noexcept
      : _q( std::move(ptc._q)), _p( std::move(ptc._p) ), _state( ptc._state ) {}

    vParticle() = delete;
    vParticle( const vParticle& ) = delete;
    template < typename E >
    vParticle( const PtcExpression<E>& ) = delete;
    template < typename E >
    vParticle( PtcExpression<E>&& ) = delete;
  };

}

#endif
