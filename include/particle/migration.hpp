#ifndef  _MIGRATION_HPP_
#define  _MIGRATION_HPP_

#include "particle/particle_expression.hpp"
#include "apt/foreach.hpp"
#include <array>
#include <vector>
#include <optional>

namespace particle {
  // cParticle uses C-native types for the sake of interfacing with MPI
  template < typename T, int DPtc, typename state_t >
  class cParticle : public PtcExpression<cParticle<T,DPtc,state_t>> {
  private:
    T _q [DPtc];
    T _p [DPtc];
    state_t _s;

  public:
    static constexpr auto Dim = DPtc;
    // TODO double check this: std::get on C-style array is supported with bound checking. At least, std::begin and std::end do.
    constexpr auto& q() noexcept { return T(&_q)[DPtc]; }
    constexpr const auto& q() const noexcept { return (const T)(&_q)[DPtc]; }

    constexpr auto& p() noexcept { return T(&_p)[DPtc]; }
    constexpr const auto& p() const noexcept { return (const T)(&_p)[DPtc]; }

    constexpr auto& state() noexcept { return _s; }
    constexpr const auto& state() const noexcept { return _s; }

    template < typename E >
    cParticle( PtcExpression<E>&& ptc ) noexcept {
      apt::foreach<0,DPtc>([](auto& a, auto& b){ std::swap(a,b);}, q(), ptc.q() );
      apt::foreach<0,DPtc>([](auto& a, auto& b){ std::swap(a,b);}, p(), ptc.p() );
      std::swap( _s, ptc.state() );
      ptc.set(flag::empty);
    }

  };
}

namespace mpi { struct InterComm; }

namespace particle {
  template < typename Vec, int DGrid, typename T >
  bool is_migrate( const Vec& q, const std::array< std::array<T, 2>, DGrid>& borders ) noexcept;

  template < typename T, int DPtc, typename state_t, int DGrid >
  void migrate ( std::vector<cParticle<T,DPtc,state_t>>& buffer,
                 const std::array< std::array<std::optional<mpi::InterComm>,2>, DGrid >& intercomms,
                 const std::array< std::array<T,2>, DGrid>& borders, unsigned int pairing_shift );
  // NOTE pairing_shift is to help achieve amortized load balancing during intercommunication. This value needs to be synchronized across at least all active processes ( preferrably all processes to avoid future sync. ) Using timestep is a good choice.
}

#endif
