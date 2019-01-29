#ifndef _PARTICLE_HPP_
#define _PARTICLE_HPP_

#include "state_codec.hpp"
#include "vector.hpp"

// TODO this limits user from passing a single particle into a function that only works on particles. Such as template <typename Ptc> void foo( Ptc ), where one cannot use Ptc&. One solution is use Ptc&&. The other is use generic lambda.
// Or let's put it this way, those functions only work with proxies. When a material particle is passed, a proxy will generated for it.
// Or, we will use Ptc&, this only disallows passing proxy&&. const Ptc& works for everything.
template < typename T >
using state_t = apt::copy_cvref_t<T, unsigned long long>;

template < typename T, std::size_t DPtc >
struct Particle : private state_codec<state_t<T>> {
  static constexpr auto Dim = DPtc;
  Vec<T, DPtc> q;
  Vec<T, DPtc> p;

  // TODO check constructor, allow virtual real converson
  Particle( const Vec<T, DPtc>& q_other, const Vec<T, DPtc>& p_other
            const state_t<T>& state_other = 0 ) noexcept
    : q(q_other), p(p_other), state_codec<state_t<T>>(state_other) {}
};


#endif
