#ifndef _PARTICLE_HPP_
#define _PARTICLE_HPP_

#include "particle/state.hpp"
#include "apt/type_traits.hpp"

// TODO this limits user from passing a single particle into a function that only works on particles. Such as template <typename Ptc> void foo( Ptc ), where one cannot use Ptc&. One solution is use Ptc&&. The other is use generic lambda.
// Or let's put it this way, those functions only work with proxies. When a material particle is passed, a proxy will generated for it.
// Or, we will use Ptc&, this only disallows passing proxy&&. const Ptc& works for everything.
namespace {
  template < typename T, std::size_t DPtc >
  struct Particle : private particle::state_t<apt::copy_cvref_t<T, particle::state_underlying_t >> {
    static constexpr auto Dim = DPtc;
    apt::Vec<T, DPtc> q;
    apt::Vec<T, DPtc> p;

    // TODO check constructor, allow virtual real converson
    template < typename... Attr,
               class = std::enable_if_t< std::is_same_v<T, std::remove_reference_t<T>>, int > >
    Particle( const apt::Vec<T, DPtc>& q_val, const apt::Vec<T, DPtc>& p_val, const Attr&... attrs ) noexcept
      : q(q_val), p(p_val),
      particle::state_t<apt::copy_cvref_t<T, particle::state_underlying_t >>( attrs... ) {}
  };
}


#endif
