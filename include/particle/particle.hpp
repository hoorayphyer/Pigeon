#ifndef _PARTICLE_HPP_
#define _PARTICLE_HPP_

#include "particle/particle_expression.hpp"
#include "apt/vec.hpp"

namespace particle {

  template < typename T, int DPtc, typename state_t >
  struct Particle : public PtcExpression< Particle<T, DPtc, state_t>, apt::Vec<T,DPtc>, state_t > {
  private:
    static_assert( 8 * sizeof( state_t ) >= 64 );

    apt::Vec<T, DPtc> _q;
    apt::Vec<T, DPtc> _p;
    state_t _state;

  public:
    static constexpr int NDim = DPtc;
    using vec_type = apt::Vec<T, DPtc>;
    using state_type = state_t;

    constexpr auto& q() noexcept { return _q; }
    constexpr const auto& q() const noexcept { return _q; }

    constexpr auto& p() noexcept { return _p; }
    constexpr const auto& p() const noexcept { return _p; }

    constexpr auto& state() noexcept { return _state; }
    constexpr const auto& state() const noexcept { return _state; }

    template < typename E1, typename E2, typename... Attrs >
    Particle( const apt::VecExpression<E1>& q, const apt::VecExpression<E2>& p, const Attrs&... attrs ) noexcept
      : _q(q), _p(p), _state(0) {
      // NOTE need this before set because of dependent base lookup
      if constexpr( sizeof...(Attrs) > 0 ) this->set(attrs...);
    }

    template < typename E >
    Particle( const PtcExpression<E>& ptc )
      : _q(ptc.q()), _p(ptc.p()), _state( ptc.state() ) {}

    template < typename E >
    Particle( PtcExpression<E>&& ptc ) noexcept
      : _q(std::move(ptc.q())), _p(std::move(ptc.p())) {
      std::swap( _state, ptc.state() );
    }
  };

}


#endif
