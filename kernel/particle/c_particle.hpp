#ifndef _C_PARTICLE_HPP_
#define _C_PARTICLE_HPP_

#include "particle/particle_expression.hpp"
#include "apt/foreach.hpp"
#include "apt/array.hpp"

namespace particle {
  // for communication
  template < typename T, template < typename > class Specs >
  class cParticle : public PtcExpression< cParticle<T, Specs>, apt::array<T,Specs<T>::Dim>, typename Specs<T>::state_type > {
  private:
    using state_t = typename Specs<T>::state_type;

    apt::array<T,Specs<T>::Dim> _q {};
    apt::array<T,Specs<T>::Dim> _p {};
    state_t _s {};
    char _extra{}; // holds extra information on communication for example

  public:
    static constexpr int NDim = Specs<T>::Dim;
    using vec_type = apt::array<T,Specs<T>::Dim>;
    using state_type = state_t;

    constexpr vec_type& q() noexcept { return _q; }
    constexpr const vec_type& q() const noexcept { return _q; }

    constexpr vec_type& p() noexcept { return _p; }
    constexpr const vec_type& p() const noexcept { return _p; }

    constexpr auto& state() noexcept { return _s; }
    constexpr const auto& state() const noexcept { return _s; }

    constexpr auto& extra() noexcept { return _extra; }
    constexpr const auto& extra() const noexcept { return _extra; }

    cParticle() = default;

    cParticle( const cParticle& ptc ) noexcept
      : _s(ptc.state()), _extra(ptc.extra()) {
      apt::foreach<0,Specs<T>::Dim>([](auto& a, const auto& b){ a = b;}, q(), ptc.q() );
      apt::foreach<0,Specs<T>::Dim>([](auto& a, const auto& b){ a = b;}, p(), ptc.p() );
    }

    cParticle( cParticle&& ptc ) noexcept {
      apt::foreach<0,Specs<T>::Dim>([](auto& a, auto& b){ std::swap(a,b);}, q(), ptc.q() );
      apt::foreach<0,Specs<T>::Dim>([](auto& a, auto& b){ std::swap(a,b);}, p(), ptc.p() );
      std::swap( _s, ptc.state() );
      std::swap( _extra, ptc.extra() );
      ptc.reset(flag::exist);
    }

    template < typename E >
    cParticle( PtcExpression<E>&& ptc ) noexcept {
      apt::foreach<0,Specs<T>::Dim>([](auto& a, auto& b){ std::swap(a,b);}, q(), ptc.q() );
      apt::foreach<0,Specs<T>::Dim>([](auto& a, auto& b){ std::swap(a,b);}, p(), ptc.p() );
      std::swap( _s, ptc.state() );
      ptc.reset(flag::exist);
    }

    cParticle& operator= ( const cParticle& ptc ) noexcept {
      apt::foreach<0,Specs<T>::Dim>([](auto& a, const auto& b){ a = b;}, q(), ptc.q() );
      apt::foreach<0,Specs<T>::Dim>([](auto& a, const auto& b){ a = b;}, p(), ptc.p() );
      _s = ptc.state();
      _extra = ptc.extra();

      return *this;
    }

    template < typename E >
    cParticle& operator= ( const PtcExpression<E>& ptc ) noexcept {
      apt::foreach<0,Specs<T>::Dim>([](auto& a, const auto& b){ a = b;}, q(), ptc.q() );
      apt::foreach<0,Specs<T>::Dim>([](auto& a, const auto& b){ a = b;}, p(), ptc.p() );
      _s = ptc.state();

      return *this;
    }

    cParticle& operator= ( cParticle&& ptc ) noexcept {
      apt::foreach<0,Specs<T>::Dim>([](auto& a, auto& b){ std::swap(a,b);}, q(), ptc.q() );
      apt::foreach<0,Specs<T>::Dim>([](auto& a, auto& b){ std::swap(a,b);}, p(), ptc.p() );
      std::swap( _s, ptc.state() );
      std::swap( _extra, ptc.extra() );
      ptc.reset(flag::exist);

      return *this;
    }

    template < typename E >
    cParticle& operator= ( PtcExpression<E>&& ptc ) noexcept {
      apt::foreach<0,Specs<T>::Dim>([](auto& a, auto& b){ std::swap(a,b);}, q(), ptc.q() );
      apt::foreach<0,Specs<T>::Dim>([](auto& a, auto& b){ std::swap(a,b);}, p(), ptc.p() );
      std::swap( _s, ptc.state() );
      ptc.reset(flag::exist);

      return *this;
    }

    // NOTE only swap handles _s_extra
    constexpr void swap( cParticle& other ) noexcept {
      apt::foreach<0,Specs<T>::Dim>([](auto& a, auto& b){ std::swap(a,b);}, q(), other.q() );
      apt::foreach<0,Specs<T>::Dim>([](auto& a, auto& b){ std::swap(a,b);}, p(), other.p() );
      std::swap( _s, other.state() );
      std::swap( _extra, other._extra );
    }

  };
}

#endif
