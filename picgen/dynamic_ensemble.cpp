#include "ensemble/dynamic_balance.cpp"
#include "pic.hpp"

namespace dye {
  using namespace pic;
  using T = real_t;
  using state_t = ptc_state_t;

  template
  void detailed_balance<T, DPtc, state_t> ( particle::array<T, DPtc, state_t>& ptcs, const mpi::Comm& intra );

  template
  void dynamic_load_balance< T, DPtc, state_t, DGrid > ( particle::map<particle::array<T, DPtc, state_t>>& particles,
                                                         std::optional<Ensemble<DGrid>>& ens_opt,
                                                         const std::optional<mpi::CartComm>& cart_opt,
                                                         unsigned int target_load );
}
