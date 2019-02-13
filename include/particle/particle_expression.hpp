#ifndef _PARTICLE_EXPRESSION_HPP_
#define _PARTICLE_EXPRESSION_HPP_

#include "particle/state.hpp"

namespace particle {
  template < typename Ptc >
  struct PtcExpression : public particle::state_expression<Ptc> {
    static constexpr auto Dim = Ptc::Dim;
    // q and p must be std::get-able
    constexpr auto& q() noexcept { return static_cast<Ptc&>(*this).q(); }
    constexpr const auto& q() const noexcept { return static_cast<const Ptc&>(*this).q(); }

    constexpr auto& p() noexcept { return static_cast<Ptc&>(*this).p(); }
    constexpr const auto& p() const noexcept { return static_cast<const Ptc&>(*this).p(); }

    constexpr auto& state() noexcept { return static_cast<Ptc&>(*this).state(); }
    constexpr const auto& state() const noexcept { return static_cast<const Ptc&>(*this).state(); }
  };

}



#endif
