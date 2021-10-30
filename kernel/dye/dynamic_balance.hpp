#ifndef _DYE_DYNAMIC_BALANCE_HPP_
#define _DYE_DYNAMIC_BALANCE_HPP_

#include "dye/ensemble.hpp"
#include "particle/array.hpp"
#include "particle/load_type.hpp"
#include "particle/map.hpp"

namespace mpi {
struct Comm;
struct CartComm;
}  // namespace mpi

namespace dye {
template <typename T, template <typename> class PtcSpecs>
void detailed_balance(particle::array<T, PtcSpecs>& ptcs,
                      const mpi::Comm& intra);

// Get ens_opt_new
template <int DGrid>
std::optional<Ensemble<DGrid>> deploy(
    particle::load_t my_tot_load,
    const std::optional<Ensemble<DGrid>>& ens_opt_old,
    const std::optional<mpi::CartComm>& cart_opt, unsigned int target_load);

// NOTE fields are not taken care of during dynamic_adjust, so data such as pair
// creation rate on each ensemble is simply lost. The solution is to do
// dynamic_ajust always afeter a data export, which is reset that kind of data
// anyway. NOTE TODOL current implementation assumes that particles have
// initialized arrays for species that may be relevant during detailed balance.
// This is due to confilict between touch create and mpi communication
template <typename T, template <typename> class PtcSpecs, int DGrid>
void dynamic_load_balance(
    particle::map<particle::array<T, PtcSpecs>>& particles,
    std::optional<Ensemble<DGrid>>& ens_opt,
    const std::optional<mpi::CartComm>& cart_opt, unsigned int target_load);
}  // namespace dye

#endif
