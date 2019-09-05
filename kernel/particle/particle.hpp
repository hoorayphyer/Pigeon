#ifndef _PARTICLE_HPP_
#define _PARTICLE_HPP_

#include "particle/particle_expression.hpp"
#include "apt/vec.hpp"

namespace particle {

  template < typename T, template < typename > class Specs >
  struct Particle : public PtcExpression< Particle<T, Specs>, apt::Vec<T,Specs<T>::Dim>, T, typename Specs<T>::state_type > {
  private:
    using state_t = typename Specs<T>::state_type;

    apt::Vec<T, Specs<T>::Dim> _q;
    apt::Vec<T, Specs<T>::Dim> _p;
    T _frac {};
    state_t _state {};

  public:
    static constexpr int NDim = Specs<T>::Dim;
    using vec_type = apt::Vec<T, Specs<T>::Dim>;
    using state_type = state_t; // TODOL

    constexpr auto& q() noexcept { return _q; }
    constexpr const auto& q() const noexcept { return _q; }

    constexpr auto& p() noexcept { return _p; }
    constexpr const auto& p() const noexcept { return _p; }

    constexpr auto& frac() noexcept { return _frac; }
    constexpr const auto& frac() const noexcept { return _frac; }

    constexpr auto& state() noexcept { return _state; }
    constexpr const auto& state() const noexcept { return _state; }

    Particle() = default;

    template < typename E1, typename E2, typename... Attrs >
    Particle( const apt::VecExpression<E1>& q, const apt::VecExpression<E2>& p, T frac, const Attrs&... attrs ) noexcept
      : _q(q), _p(p), _frac(frac) {
      this->set(flag::exist);
      // NOTE need "this" before set because of dependent base lookup
      if constexpr( sizeof...(Attrs) > 0 ) this->set(attrs...);
    }

    // this version is to enable brace-init
    template < typename... Attrs >
    Particle( const vec_type& q, const vec_type& p, T frac, const Attrs&... attrs ) noexcept
      : _q(q), _p(p), _frac(frac) {
      this->set(flag::exist);
      // NOTE need "this" before set because of dependent base lookup
      if constexpr( sizeof...(Attrs) > 0 ) this->set(attrs...);
    }

    template < typename E >
    Particle( const PtcExpression<E>& ptc )
      : _q(ptc.q()), _p(ptc.p()), _frac( ptc.frac() ), _state( ptc.state() ) {}

    template < typename E >
    Particle( PtcExpression<E>&& ptc ) noexcept
      : _q(std::move(ptc.q())), _p(std::move(ptc.p())), _frac( ptc.frac() ), _state( ptc.state() ) {
      ptc.reset(flag::exist);
    }

    constexpr void swap( Particle& other ) noexcept {
      apt::foreach<0,Specs<T>::Dim>([](auto& a, auto& b){ std::swap(a,b);}, q(), other.q() );
      apt::foreach<0,Specs<T>::Dim>([](auto& a, auto& b){ std::swap(a,b);}, p(), other.p() );
      std::swap( frac(), other.frac() );
      std::swap( state(), other.state() );
    }
  };

}


#endif
