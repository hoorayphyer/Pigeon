#include "dye/dynamic_balance_impl.hpp"
#include "pic.hpp"

namespace dye {
using namespace pic;
template void detailed_balance(particle::array<real_t, particle::Specs>& ptcs,
                               const mpi::Comm& intra);

template void dynamic_load_balance(
    particle::map<particle::array<real_t, particle::Specs>>& particles,
    std::optional<Ensemble<DGrid>>& ens_opt,
    const std::optional<mpi::CartComm>& cart_opt, unsigned int target_load);
}  // namespace dye
