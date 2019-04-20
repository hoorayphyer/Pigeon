#ifndef _DYE_DYNAMIC_BALANCE_HPP_
#define _DYE_DYNAMIC_BALANCE_HPP_

#include "particle/array.hpp"
#include "particle/map.hpp"

#include "ensemble.hpp"

namespace mpi {
  struct Comm;
  struct CartComm;
}

namespace dye {
  template < typename T, template < typename > class PtcSpecs >
  void detailed_balance ( particle::array<T, PtcSpecs>& ptcs,
                          const mpi::Comm& intra );

  // NOTE fields are not taken care of during dynamic_adjust, so data such as pair creation rate on each ensemble is simply lost. The solution is to do dynamic_ajust always afeter a data export, which is reset that kind of data anyway.
  template < typename T, template < typename > class PtcSpecs, int DGrid >
  void dynamic_load_balance( particle::map<particle::array<T,PtcSpecs>>& particles,
                             std::optional<Ensemble<DGrid>>& ens_opt,
                             const std::optional<mpi::CartComm>& cart_opt,
                             unsigned int target_load );
}

#endif

