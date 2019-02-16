#ifndef _PARTICLE_EXPRESSION_HPP_
#define _PARTICLE_EXPRESSION_HPP_

#include "particle/state.hpp"

namespace particle {
  template < typename Ptc, typename Vec = typename Ptc::vec_type, typename State = typename Ptc::state_type >
  struct PtcExpression : public particle::state_expression<Ptc> {
    static constexpr int Dim = Ptc::Dim;
    // q and p must be std::get-able
    constexpr Vec& q() noexcept { return static_cast<Ptc&>(*this).q(); }
    constexpr const Vec& q() const noexcept { return static_cast<const Ptc&>(*this).q(); }

    constexpr Vec& p() noexcept { return static_cast<Ptc&>(*this).p(); }
    constexpr const Vec& p() const noexcept { return static_cast<const Ptc&>(*this).p(); }

    constexpr State& state() noexcept { return static_cast<Ptc&>(*this).state(); }
    constexpr const State& state() const noexcept { return static_cast<const Ptc&>(*this).state(); }
  };

}



#endif
