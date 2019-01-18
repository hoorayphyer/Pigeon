#ifndef _PARTICLE_HPP_
#define _PARTICLE_HPP_

#include "particle_state_codec.hpp"
#include "vector.hpp"


// Particle is only a proxy
// TODO this limits user from passing a single particle into a function that only works on particles. Such as template <typename Ptc> void foo( Ptc ), where one cannot use Ptc&. One solution is use Ptc&&. The other is use generic lambda.
// Or let's put it this way, those functions only work with proxies. When a material particle is passed, a proxy will generated for it.
// Or, we will use Ptc&, this only disallows passing proxy&&. const Ptc& works for everything.

template < std::size_t Dim_Ptc, typename T, typename state_t >
struct Particle : private state_codec<state_t> {
  static_assert( std::is_reference_v<T> == std::is_reference_v<state_t> );
public:
  static constexpr bool is_proxy = std::is_reference_v<T>;
  static constexpr auto Dim = Dim_Ptc;
  Vec<T, Dim_Ptc> q;
  Vec<T, Dim_Ptc> p;

  Particle( Vec<T, Dim_Ptc> qq, Vec<T, Dim_Ptc> pp, state_t state ) noexcept
    : q(std::move(qq)), p(std::move(pp)), state_codec<state_t>(state) {}
};


#endif
