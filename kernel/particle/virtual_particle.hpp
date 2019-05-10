#ifndef _VIRTUAL_PARTICLE_HPP_
#define _VIRTUAL_PARTICLE_HPP_

#include "particle/particle_expression.hpp"
#include "apt/virtual_vec.hpp"

namespace particle {
  template < typename T, template < typename > class Specs >
  struct vParticle : public PtcExpression< vParticle<T, Specs>, apt::vVec<T,Specs<T>::Dim>, typename Specs<T>::state_type > {
  private:
    using state_t = typename Specs<T>::state_type;

    apt::vVec<T, Specs<T>::Dim> _q;
    apt::vVec<T, Specs<T>::Dim> _p;
    state_t& _state;

  public:
    static constexpr int NDim = Specs<T>::Dim;
    using vec_type = apt::vVec<T, Specs<T>::Dim>;
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
    constexpr vParticle& operator= ( const PtcExpression<E>& ptc ) noexcept {
      _q = ptc.q();
      _p = ptc.p();
      _state = ptc.state();
      return *this;
    }

  };

}

#endif
