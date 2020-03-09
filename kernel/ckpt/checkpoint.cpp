#include "ckpt/checkpoint_impl.hpp"
#include "pic.hpp"

namespace ckpt {
  using namespace pic;

  template
  std::string save_checkpoint( std::string prefix, const int num_parts,
                               const std::optional<dye::Ensemble<DGrid>>& ens_opt,
                               int timestep,
                               const field::Field<real_t, 3, DGrid>& E,
                               const field::Field<real_t, 3, DGrid>& B,
                               const particle::map<particle::array<real_t,particle::Specs>>& particles,
                               const particle::map<particle::Properties>& properties
                               );

  template
  int load_checkpoint( std::string dir,
                       std::optional<dye::Ensemble<DGrid>>& ens_opt,
                       const std::optional<mpi::CartComm>& cart_opt,
                       field::Field<real_t, 3, DGrid>& E,
                       field::Field<real_t, 3, DGrid>& B,
                       particle::map<particle::array<real_t,particle::Specs>>& particles,
                       const particle::map<particle::Properties>& properties,
                       int target_load
                       );

  template
  std::string save_tracing( std::string prefix, const int num_parts,
                            const std::optional<dye::Ensemble<DGrid>>& ens_opt,
                            int timestep,
                            const particle::map<particle::array<real_t,particle::Specs>>& particles
                            );
}
