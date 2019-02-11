#ifndef _PARTICLE_HPP_
#define _PARTICLE_HPP_

#include "particle/state.hpp"

namespace particle {
  template < typename T, int DPtc, template < typename, int > class Vec, typename state_t >
  struct Particle : private particle::state_expression<Particle<T, DPtc, Vec, state_t>> {
  private:
    static_assert( 8 * sizeof( std::remove_reference_t<state_t>) >= 64 );
    state_t _state;
    inline auto& _state() noexcept { return _state; }

  public:
    static constexpr auto Dim = DPtc;
    friend class particle::state_expression<Particle<T, DPtc, Vec, state_t>>;

    Vec<T, DPtc> q;
    Vec<T, DPtc> p;

    template < typename... Attr >
    Particle( const Vec<T, DPtc>& q_val, const Vec<T, DPtc>& p_val,
              state_t state_val, const Attr&... attrs ) noexcept
      : q(q_val), p(p_val), _state(state_val) {
      if constexpr( sizeof...(Attr) > 0 ) set(attrs...);
    }
  };

}


#endif
