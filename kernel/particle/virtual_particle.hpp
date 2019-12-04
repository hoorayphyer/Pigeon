#ifndef _VIRTUAL_PARTICLE_HPP_
#define _VIRTUAL_PARTICLE_HPP_

#include "particle/particle_expression.hpp"
#include "apt/virtual_vec.hpp"

namespace particle {
  template < typename T, template < typename > class Specs >
  struct vParticle : public PtcExpression< vParticle<T, Specs>, apt::vVec<T,Specs<T>::Dim>, T, typename Specs<T>::state_type > {
  private:
    using state_t = typename Specs<T>::state_type;

    apt::vVec<T, Specs<T>::Dim> _q;
    apt::vVec<T, Specs<T>::Dim> _p;
    T& _frac;
    state_t& _state;

  public:
    static constexpr int NDim = Specs<T>::Dim;
    using vec_type = apt::vVec<T, Specs<T>::Dim>;
    using state_type = state_t;

    constexpr auto& q() noexcept { return _q; }
    constexpr const auto& q() const noexcept { return _q; }
    constexpr auto& q(int i) noexcept { return _q[i]; }
    constexpr const auto& q(int i) const noexcept { return _q[i]; }

    constexpr auto& p() noexcept { return _p; }
    constexpr const auto& p() const noexcept { return _p; }
    constexpr auto& p(int i) noexcept { return _p[i]; }
    constexpr const auto& p(int i) const noexcept { return _p[i]; }

    constexpr auto& frac() noexcept { return _frac; }
    constexpr const auto& frac() const noexcept { return _frac; }

    constexpr auto& state() noexcept { return _state; }
    constexpr const auto& state() const noexcept { return _state; }

    template < typename Q, typename P >
    vParticle( Q&& q, P&& p, T& frac, state_t& state ) noexcept
      : _q(std::forward<Q>(q)),
        _p(std::forward<P>(p)),
        _frac(frac),
        _state(state) {}

    vParticle() = delete;
    vParticle( const vParticle& ) = delete;
    vParticle( vParticle&& ptc ) noexcept = default;

    template < typename E >
    constexpr vParticle& operator= ( const PtcExpression<E>& ptc ) noexcept {
      _q = ptc.q();
      _p = ptc.p();
      _frac = ptc.frac();
      _state = ptc.state();
      return *this;
    }

    template < typename E >
    constexpr vParticle& operator= ( PtcExpression<E>&& ptc ) noexcept {
      for ( int i = 0; i < NDim; ++i ) {
        std::swap( _q[i], ptc.q(i) );
        std::swap( _p[i], ptc.p(i) );
      }
      std::swap( _frac, ptc.frac() );
      std::swap( _state, ptc.state() );
      ptc.reset(flag::exist);
      return *this;
    }
  };

}

#endif
