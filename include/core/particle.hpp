#ifndef _PARTICLE_HPP_
#define _PARTICLE_HPP_

#include "particle_state_codec.hpp"
#include "vector.hpp"

// TODO this limits user from passing a single particle into a function that only works on particles. Such as template <typename Ptc> void foo( Ptc ), where one cannot use Ptc&. One solution is use Ptc&&. The other is use generic lambda.
// Or let's put it this way, those functions only work with proxies. When a material particle is passed, a proxy will generated for it.
// Or, we will use Ptc&, this only disallows passing proxy&&. const Ptc& works for everything.


template < typename T, std::size_t DPtc, typename state_t >
struct Particle : private state_codec<apt::copy_constref_t<T, state_t>> {
  static constexpr auto Dim = DPtc;
  Vec<T, Dim_Ptc> q;
  Vec<T, Dim_Ptc> p;

  // TODO check constructor
  Particle( Vec<T, DPtc> qq, Vec<T, DPtc> pp,
            apt::copy_constref_t<T, state_t> state ) noexcept
    : q(std::move(qq)), p(std::move(pp)), state_codec<state_t>(state) {}
};


#endif
