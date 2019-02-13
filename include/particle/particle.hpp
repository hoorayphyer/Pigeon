#ifndef _PARTICLE_HPP_
#define _PARTICLE_HPP_

#include "particle/particle_expression.hpp"
#include "apt/vec.hpp"

namespace particle {
  template < typename T, int DPtc, typename state_t >
  struct vParticle : public PtcExpression< vParticle<T, DPtc, state_t> > {
  private:
    static_assert( 8 * sizeof( state_t) >= 64 );

    apt::vVec<T, DPtc> _q;
    apt::vVec<T, DPtc> _p;
    state_t& _state;

  public:
    static constexpr auto Dim = DPtc;

    constexpr auto& q() noexcept { return _q; }
    constexpr const auto& q() const noexcept { return _q; }

    constexpr auto& p() noexcept { return _p; }
    constexpr const auto& p() const noexcept { return _p; }

    constexpr auto& state() noexcept { return _state; }
    constexpr const auto& state() const noexcept { return _state; }

    vParticle( apt::vVec<T, DPtc>&& q, apt::vVec<T, DPtc>&& p, state_t& state ) noexcept
      : _q(std::move(q)), _p(std::move(p)), _state(state) {}

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

namespace particle {

  template < typename T, int DPtc, typename state_t >
  struct Particle : public PtcExpression< Particle<T, DPtc, state_t> > {
  private:
    static_assert( 8 * sizeof( state_t ) >= 64 );

    apt::Vec<T, DPtc> _q;
    apt::Vec<T, DPtc> _p;
    state_t _state;

  public:
    static constexpr auto Dim = DPtc;

    constexpr auto& q() noexcept { return _q; }
    constexpr const auto& q() const noexcept { return _q; }

    constexpr auto& p() noexcept { return _p; }
    constexpr const auto& p() const noexcept { return _p; }

    constexpr auto& state() noexcept { return _state; }
    constexpr const auto& state() const noexcept { return _state; }

    template < typename... Attr >
    Particle( const apt::Vec<T, DPtc>& q, const apt::Vec<T, DPtc>& p, const Attr&... attrs ) noexcept
      : _q(q), _p(p), _state(0) {
      if constexpr( sizeof...(Attr) > 0 ) set(attrs...);
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
